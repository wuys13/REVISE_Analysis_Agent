#!/usr/bin/env python3
"""Prepare the standardized agent input payload from quantitative analysis results.

Specification: docs/delivery_framework/05_agent_data_contract_and_payload.md
Extracts key quantitative metrics from tables/ and detailed_data/ to supply
the synthesis Agent with a clean, bounded factual payload for generating narrative texts.
Strictly avoids synthetic formulas, mock statistics, or hardcoded fake metrics.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path
import pandas as pd


def prepare_agent_payload(sample_dir: Path, output_file: Path | None = None) -> dict:
    sample_dir = sample_dir.resolve()
    analysis_dir = sample_dir / "analysis"
    
    if not analysis_dir.exists():
        print(f"[!] Error: analysis directory not found at {analysis_dir}", file=sys.stderr)
        sys.exit(1)

    t1_tables = analysis_dir / "01_bioinformatic_analysis/tables"
    t1_details = analysis_dir / "01_bioinformatic_analysis/detailed_data"
    t2_tables = analysis_dir / "02_reconstruction_impact/tables"
    t3_tables = analysis_dir / "03_biological_discoveries/tables"

    # Define required files for full generation
    required_tables = {
        "unmasking_evidence": t2_tables / "raw_vs_svc_unmasking_evidence.csv",
        "systemic_5layers": t2_tables / "systemic_impact_5layers.csv",
        "notable_subtypes": t3_tables / "notable_subtypes_summary.csv",
        "pathway_scores": t1_tables / "pathway_scores_by_subtype.csv",
        "significant_cci": t1_tables / "significant_lr_pairs.csv",
        "colocalizations": t1_tables / "colocalization_zscores.csv",
        "hypotheses_matrix": t3_tables / "hypotheses_evidence_matrix.csv"
    }

    missing_tables = [k for k, p in required_tables.items() if not p.exists()]
    generation_status = "complete" if not missing_tables else "partial"

    summary_path = sample_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not isinstance(summary, dict) or not summary.get("sample_id"):
        raise ValueError("summary.json must register this sample's sample_id")
    reconstruction = summary.get("reconstruction_summary", {})
    frames = {
        key: pd.read_csv(path) for key, path in required_tables.items() if path.exists()
    }
    deg_path = t1_details / "focal_subtypes_top_degs.csv"
    deg_frame = pd.read_csv(deg_path) if deg_path.exists() else None
    registry_path = sample_dir / "source_registry.json"
    registry = json.loads(registry_path.read_text(encoding="utf-8")) if registry_path.exists() else None
    evidence_files = {
        key: {"path": str(path.relative_to(sample_dir)),
              "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "rows": len(frames[key])}
        for key, path in required_tables.items() if key in frames
    }
    if deg_frame is not None:
        evidence_files["focal_top_degs"] = {
            "path": str(deg_path.relative_to(sample_dir)),
            "sha256": hashlib.sha256(deg_path.read_bytes()).hexdigest(),
            "rows": len(deg_frame),
        }
    for key, path in (("sample_summary", summary_path), ("source_registry", registry_path)):
        if path.exists():
            evidence_files[key] = {
                "path": str(path.relative_to(sample_dir)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }

    def row_count(key: str) -> int | None:
        return len(frames[key]) if key in frames else None

    # 1. Metadata and provenance registration
    payload = {
        "sample_metadata": {
            "sample_id": summary["sample_id"],
            "cancer_type": summary.get("cancer_type"),
            "modality": summary.get("modality"),
            "total_cells": reconstruction.get("total_reconstructed_cells"),
            "lineages_covered": reconstruction.get("target_lineages"),
            "subclusters_universe": reconstruction.get("subclusters_universe"),
            "subclusters_per_lineage": reconstruction.get("subclusters_per_lineage"),
            "paper_reference": summary.get("paper_alignment"),
            "provenance_registry": registry,
            "evaluation_scope": {
                "focal_subtypes_count": row_count("notable_subtypes"),
                "top_degs_rows": len(deg_frame) if deg_frame is not None else None,
                "hallmark_ora_rows": row_count("pathway_scores"),
                "significant_cci_pairs": row_count("significant_cci"),
                "hypotheses_count": row_count("hypotheses_matrix")
            },
            "pipeline_status": generation_status,
            "missing_inputs": missing_tables,
            "source_status": summary.get("status", "unverified"),
            "acceptance_level": summary.get("acceptance_level", "unverified"),
            "evidence_files": evidence_files
        },
        "quantitative_evidence": {}
    }

    def records(key: str) -> list[dict]:
        frame = frames[key]
        return frame.astype(object).where(frame.notna(), None).to_dict(orient="records")

    # 2. Extract Tier 2 Unmasking evidence & 5-layer systemic impact
    payload["quantitative_evidence"]["unmasking_impact"] = {}
    if required_tables["unmasking_evidence"].exists():
        payload["quantitative_evidence"]["unmasking_impact"]["key_comparisons"] = records("unmasking_evidence")
    if required_tables["systemic_5layers"].exists():
        payload["quantitative_evidence"]["unmasking_impact"]["systemic_impact_5layers"] = records("systemic_5layers")

    # 3. Extract saved notable subclusters
    if required_tables["notable_subtypes"].exists():
        payload["quantitative_evidence"]["top_subtypes_multimetric"] = records("notable_subtypes")

    # 4. Extract Tier 1 Pathways (ORA Hallmark)
    if required_tables["pathway_scores"].exists():
        payload["quantitative_evidence"]["pathways_ora_hallmark"] = records("pathway_scores")

    # 5. Extract Tier 1 Spatial colocalization & CCI
    if required_tables["significant_cci"].exists():
        payload["quantitative_evidence"]["cell_cell_interactions"] = records("significant_cci")
    if required_tables["colocalizations"].exists():
        payload["quantitative_evidence"]["spatial_colocalizations"] = records("colocalizations")

    # 6. Extract Tier 3 Hypotheses Evidence Matrix
    if required_tables["hypotheses_matrix"].exists():
        payload["quantitative_evidence"]["hypotheses_evidence"] = records("hypotheses_matrix")

    view_paths = [analysis_dir / slot / "FINDINGS.md" for slot in (
        "01_bioinformatic_analysis", "02_reconstruction_impact", "03_biological_discoveries"
    )]
    view_paths.extend(sorted((analysis_dir / "03_biological_discoveries/hypotheses").glob("*.md")))
    view_paths.extend([sample_dir / "README.md", sample_dir / "report.html"])
    payload["generation_targets"] = [
        str(path.relative_to(sample_dir)) for path in view_paths if path.is_file()
    ]

    target_json = output_file or (sample_dir / "agent_input_payload.json")
    serialized = json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False)
    target_json.write_text(serialized, encoding="utf-8")
    
    print(f"[✓] Generated agent input payload: {target_json} (status: {generation_status})")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Agent Input Payload")
    parser.add_argument("--sample-dir", type=Path, default=Path(__file__).resolve().parents[2] / "output/CRC/P2_CRC_Xenium")
    args = parser.parse_args()

    prepare_agent_payload(args.sample_dir)
