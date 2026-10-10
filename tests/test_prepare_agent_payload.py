"""Published payload metadata must describe its own sample and saved tables."""

import importlib.util
import hashlib
import json
from pathlib import Path
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_agent_payload.py"
spec = importlib.util.spec_from_file_location("prepare_agent_payload", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_payload_uses_registered_sample_and_actual_coverage(tmp_path):
    sample = tmp_path / "another_sample"
    tables = sample / "analysis" / "01_bioinformatic_analysis" / "tables"
    tables.mkdir(parents=True)
    (sample / "summary.json").write_text(json.dumps({
        "sample_id": "registered_id", "cancer_type": "BRCA", "modality": "Visium",
        "status": "partial", "scientific_acceptance": False,
        "reconstruction_summary": {
            "total_reconstructed_cells": 123, "target_lineages": ["Tumor"],
        },
    }))
    (tables / "pathway_scores_by_subtype.csv").write_text(
        "subcluster,pathway,score\nTumor_1,PROGRAM,2.0\n"
    )
    details = tables.parent / "detailed_data"
    details.mkdir()
    deg_path = details / "focal_subtypes_top_degs.csv"
    deg_path.write_text("subcluster,gene\nTumor_1,GENE_A\n")
    payload = module.prepare_agent_payload(sample)
    metadata = payload["sample_metadata"]
    assert metadata["sample_id"] == "registered_id"
    assert metadata["cancer_type"] == "BRCA"
    assert metadata["modality"] == "Visium"
    assert metadata["total_cells"] == 123
    assert metadata["lineages_covered"] == ["Tumor"]
    assert metadata["subclusters_universe"] is None
    assert metadata["evaluation_scope"]["hallmark_ora_rows"] == 1
    assert metadata["evaluation_scope"]["selected_cci_rows"] is None
    assert "significant_cci_pairs" not in metadata["evaluation_scope"]
    assert metadata["scientific_acceptance"] is False
    assert metadata["pipeline_status"] == "partial"
    assert "selected_lr_candidates" in metadata["missing_inputs"]
    assert "literature_spatial_context" in metadata["missing_inputs"]
    assert "significant_lr_pairs" not in metadata["missing_inputs"]
    assert metadata["evidence_files"]["focal_top_degs"] == {
        "path": "analysis/01_bioinformatic_analysis/detailed_data/focal_subtypes_top_degs.csv",
        "sha256": hashlib.sha256(deg_path.read_bytes()).hexdigest(), "rows": 1,
    }


def test_missing_measurement_stays_null_in_saved_payload(tmp_path):
    sample = tmp_path / "sample"
    tables = sample / "analysis" / "01_bioinformatic_analysis" / "tables"
    tables.mkdir(parents=True)
    (sample / "summary.json").write_text('{"sample_id": "sample"}')
    (tables / "pathway_scores_by_subtype.csv").write_text(
        "subcluster,pathway,score\nTumor_1,PROGRAM,\n"
    )
    module.prepare_agent_payload(sample)
    saved = json.loads((sample / "agent_input_payload.json").read_text())
    assert saved["quantitative_evidence"]["pathways_ora_hallmark"][0]["score"] is None


def test_invalid_table_does_not_overwrite_existing_payload(tmp_path):
    sample = tmp_path / "sample"
    tables = sample / "analysis" / "01_bioinformatic_analysis" / "tables"
    tables.mkdir(parents=True)
    (sample / "summary.json").write_text('{"sample_id": "sample"}')
    (tables / "pathway_scores_by_subtype.csv").write_text(
        "subcluster,pathway,score\nTumor_1,PROGRAM,inf\n"
    )
    payload_path = sample / "agent_input_payload.json"
    payload_path.write_text('{"previous": "valid"}')
    with pytest.raises(ValueError):
        module.prepare_agent_payload(sample)
    assert json.loads(payload_path.read_text()) == {"previous": "valid"}


def test_payload_separates_quantitative_and_contextual_evidence(tmp_path):
    sample = tmp_path / "sample"
    analysis = sample / "analysis"
    table_paths = {
        "historical_count_comparisons": (
            analysis / "02_reconstruction_impact" / "tables" / "raw_vs_svc_unmasking_evidence.csv",
            "comparison,fold\nT_1,2.5\n",
        ),
        "systemic_impact_5layers": (
            analysis / "02_reconstruction_impact" / "tables" / "systemic_impact_5layers.csv",
            "layer,description\n1,context\n",
        ),
        "focal_subtype_descriptions": (
            analysis / "03_biological_discoveries" / "tables" / "notable_subtypes_summary.csv",
            "subcluster,description\nT_1,desc\n",
        ),
        "pathways_ora_hallmark": (
            analysis / "01_bioinformatic_analysis" / "tables" / "pathway_scores_by_subtype.csv",
            "subcluster,score\nT_1,\n",
        ),
        "selected_lr_candidates": (
            analysis / "01_bioinformatic_analysis" / "tables" / "selected_lr_candidates.csv",
            "pair,score\nL-R,3.0\n",
        ),
        "literature_spatial_context": (
            analysis / "01_bioinformatic_analysis" / "tables" / "literature_spatial_context.csv",
            "pair,evidence\nT_1|T_2,observed\n",
        ),
        "hypothesis_candidates": (
            analysis / "03_biological_discoveries" / "tables" / "hypotheses_evidence_matrix.csv",
            "hypothesis_id,score\nH1,0.8\n",
        ),
    }
    for path, content in table_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    (sample / "summary.json").write_text(json.dumps({
        "sample_id": "sample", "scientific_acceptance": False,
        "reconstruction_summary": {},
    }))

    payload = module.prepare_agent_payload(sample)
    metadata = payload["sample_metadata"]

    assert metadata["pipeline_status"] == "complete"
    assert metadata["missing_inputs"] == []
    assert metadata["scientific_acceptance"] is False
    assert metadata["evaluation_scope"]["selected_cci_rows"] == 1
    assert set(payload["quantitative_evidence"]) == {
        "historical_count_comparisons", "pathways_ora_hallmark", "preselected_lr_rows",
    }
    assert set(payload["contextual_evidence"]) == {
        "focal_subtype_descriptions", "systemic_impact_5layers",
        "literature_spatial_colocalizations", "hypothesis_candidates",
    }
    assert payload["quantitative_evidence"]["pathways_ora_hallmark"][0]["score"] is None
    assert payload["quantitative_evidence"]["historical_count_comparisons"][0]["fold"] == 2.5
    assert payload["quantitative_evidence"]["preselected_lr_rows"][0]["score"] == 3.0
    assert payload["contextual_evidence"]["focal_subtype_descriptions"][0]["description"] == "desc"
    assert "significant_cci_pairs" not in json.dumps(payload)
    assert "target_lineage_programs" not in json.dumps(payload)

    for key, (path, _) in table_paths.items():
        evidence = metadata["evidence_files"][key]
        assert evidence["path"] == str(path.relative_to(sample))
        assert evidence["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
        assert evidence["rows"] == 1


def test_missing_scientific_acceptance_remains_null(tmp_path):
    sample = tmp_path / "sample"
    tables = sample / "analysis" / "01_bioinformatic_analysis" / "tables"
    tables.mkdir(parents=True)
    (sample / "summary.json").write_text('{"sample_id": "sample"}')
    (tables / "pathway_scores_by_subtype.csv").write_text(
        "subcluster,pathway,score\nTumor_1,PROGRAM,1.0\n"
    )

    payload = module.prepare_agent_payload(sample)

    assert payload["sample_metadata"]["scientific_acceptance"] is None
