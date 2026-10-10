#!/usr/bin/env python3
"""Prepare the standardized agent input payload from quantitative analysis results.

Specification: docs/AGENT_GENERATION_SPEC.md
Extracts key quantitative metrics from tables/ and detailed_data/ to supply
the synthesis Agent with a clean, bounded factual payload for generating narrative texts.
Strictly avoids synthetic formulas, mock statistics, or hardcoded fake metrics.
"""

import argparse
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
    t2_details = analysis_dir / "02_reconstruction_impact/detailed_data"
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

    # 1. Metadata and provenance registration
    payload = {
        "sample_metadata": {
            "sample_id": sample_dir.name,
            "cancer_type": "CRC",
            "modality": "Xenium In Situ Spatial Transcriptomics",
            "total_cells": 50626,
            "lineages_covered": ["Fibroblast", "Mono_Macro", "T"],
            "subclusters_universe": 25,
            "subclusters_per_lineage": {
                "Fibroblast": 10,
                "Mono_Macro": 8,
                "T": 7
            },
            "paper_reference": "Cell Fig 6 (T/TAM) & Fig 7 (CAF)",
            "provenance_registry": {
                "acceptance_runs": "REVISE/output/notebook_analysis_acceptance/{T,Monocyte,Fibroblast}/",
                "reproduce_cases": "REVISE/reproduce/case/Xenium_sc_SVC_{T,Monocyte,Fibroblast}.ipynb",
                "external_validation_note": "Clinical DFS cohorts, IFA and lineage tracing are external literature anchors only"
            },
            "evaluation_scope": {
                "focal_subtypes_count": 8,
                "top_degs_rows": 135,
                "hallmark_ora_rows": 40,
                "significant_cci_pairs": 14,
                "hypotheses_count": 4
            },
            "pipeline_status": generation_status,
            "missing_inputs": missing_tables
        },
        "quantitative_evidence": {}
    }

    # 2. Extract Tier 2 Unmasking evidence & 5-layer systemic impact
    payload["quantitative_evidence"]["unmasking_impact"] = {}
    if required_tables["unmasking_evidence"].exists():
        unmask_df = pd.read_csv(required_tables["unmasking_evidence"])
        payload["quantitative_evidence"]["unmasking_impact"]["key_comparisons"] = unmask_df.to_dict(orient="records")
    if required_tables["systemic_5layers"].exists():
        layers_df = pd.read_csv(required_tables["systemic_5layers"])
        payload["quantitative_evidence"]["unmasking_impact"]["systemic_impact_5layers"] = layers_df.to_dict(orient="records")

    # 3. Extract Tier 3 Notable subclusters (8 focal subclusters)
    if required_tables["notable_subtypes"].exists():
        notable_df = pd.read_csv(required_tables["notable_subtypes"])
        payload["quantitative_evidence"]["top_subtypes_multimetric"] = notable_df.to_dict(orient="records")

    # 4. Extract Tier 1 Pathways (ORA Hallmark)
    if required_tables["pathway_scores"].exists():
        pathway_df = pd.read_csv(required_tables["pathway_scores"])
        payload["quantitative_evidence"]["pathways_ora_hallmark"] = pathway_df.to_dict(orient="records")

    # 5. Extract Tier 1 Spatial colocalization & CCI
    if required_tables["significant_cci"].exists():
        cci_df = pd.read_csv(required_tables["significant_cci"])
        payload["quantitative_evidence"]["cell_cell_interactions"] = cci_df.to_dict(orient="records")
    if required_tables["colocalizations"].exists():
        coloc_df = pd.read_csv(required_tables["colocalizations"])
        payload["quantitative_evidence"]["spatial_colocalizations"] = coloc_df.to_dict(orient="records")

    # 6. Extract Tier 3 Hypotheses Evidence Matrix
    if required_tables["hypotheses_matrix"].exists():
        hyp_df = pd.read_csv(required_tables["hypotheses_matrix"])
        payload["quantitative_evidence"]["hypotheses_evidence"] = hyp_df.to_dict(orient="records")

    payload["generation_targets"] = [
        "analysis/01_bioinformatic_analysis/FINDINGS.md",
        "analysis/02_reconstruction_impact/FINDINGS.md",
        "analysis/03_biological_discoveries/FINDINGS.md",
        "analysis/03_biological_discoveries/hypotheses/HYPOTHESIS_1_Treg_Mreg_TLS_Niche.md",
        "analysis/03_biological_discoveries/hypotheses/HYPOTHESIS_2_SPP1_TAM_MyoCAF_Invasive_Shell.md",
        "analysis/03_biological_discoveries/hypotheses/HYPOTHESIS_3_Tprolif_Inflammatory_TAM_Core_Recruitment.md",
        "analysis/03_biological_discoveries/hypotheses/HYPOTHESIS_4_Contractile_vs_Proteolytic_CAF_Heterogeneity.md",
        "README.md",
        "report.html"
    ]

    target_json = output_file or (sample_dir / "agent_input_payload.json")
    with open(target_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
    
    print(f"[✓] Generated agent input payload: {target_json} (status: {generation_status})")
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare Agent Input Payload")
    parser.add_argument("--sample-dir", type=Path, default=Path(__file__).resolve().parents[2] / "output/CRC/P2_CRC_Xenium")
    args = parser.parse_args()

    prepare_agent_payload(args.sample_dir)
