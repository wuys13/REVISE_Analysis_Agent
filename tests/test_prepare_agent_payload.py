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
        "status": "partial", "reconstruction_summary": {
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
    assert metadata["evaluation_scope"]["significant_cci_pairs"] is None
    assert metadata["pipeline_status"] == "partial"
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
