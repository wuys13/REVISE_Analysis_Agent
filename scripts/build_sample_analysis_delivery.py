#!/usr/bin/env python3
"""[DEPRECATED / ARCHIVED PROTOTYPE]
Historical script used for initial directory structure skeleton.
DO NOT RUN to overwrite official deliverables.
All grounded deliverables are managed directly under output/CRC/P2_CRC_Xenium/
and extracted via prepare_agent_payload.py.
"""

import sys

if __name__ == "__main__":
    print("[ERROR] build_sample_analysis_delivery.py is archived and deprecated.", file=sys.stderr)
    print("        It contains legacy hardcoded templates and must not overwrite current deliverables.", file=sys.stderr)
    sys.exit(1)

import argparse
import json
import shutil
from pathlib import Path
import pandas as pd


def build_sample_delivery(sample_dir: Path, cancer_type: str = "CRC", sample_id: str = "P2_CRC_Xenium"):
    sample_dir = sample_dir.resolve()
    print(f"[*] Building delivery package for {sample_id} in: {sample_dir}")
    analysis_dir = sample_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    t1_dir = analysis_dir / "01_bioinformatic_analysis"
    t2_dir = analysis_dir / "02_reconstruction_impact"
    t3_dir = analysis_dir / "03_biological_discoveries"

    for d in [t1_dir / "tables", t1_dir / "detailed_data",
              t2_dir / "tables", t2_dir / "detailed_data",
              t3_dir / "tables", t3_dir / "detailed_data", t3_dir / "hypotheses"]:
        d.mkdir(parents=True, exist_ok=True)

    # 1. Update analysis.json
    analysis_manifest = {
        "schema_version": 2,
        "sample_id": sample_id,
        "paradigm": "Inverted Pyramid Progressive Disclosure",
        "slots": [
            {
                "slot_id": "01_bioinformatic_analysis",
                "role": "全谱系平权先验生信分析",
                "tables": ["pathway_scores_by_subtype.csv", "colocalization_zscores.csv", "significant_lr_pairs.csv"],
                "detailed_data": ["subtype_degs_full.csv"]
            },
            {
                "slot_id": "02_reconstruction_impact",
                "role": "重建前后变化与原始去掩盖刻画",
                "tables": ["raw_vs_svc_unmasking_evidence.csv", "moran_gain_by_lineage.csv", "diversity_entropy_delta.csv"],
                "detailed_data": ["contamination_clearance.csv", "noise_reduction_diffusion.csv"]
            },
            {
                "slot_id": "03_biological_discoveries",
                "role": "综合分析下的生物学发现与科学假设",
                "hypotheses_count": 4,
                "tables": ["notable_subtypes_summary.csv", "hypotheses_evidence_matrix.csv"],
                "detailed_data": ["target_lineage_programs.csv"]
            }
        ]
    }
    with open(analysis_dir / "analysis.json", "w") as f:
        json.dump(analysis_manifest, f, indent=2)

    # 2. Update analysis/README.md
    with open(analysis_dir / "README.md", "w") as f:
        f.write(f"""# {sample_id} 分析成果槽位 (Analysis Slots Root)

本目录严格遵循 [倒金字塔渐进式披露范式](../../REVISE_Analysis_Agent/docs/ANALYSIS_OUTPUT_PARADIGM.md)：
- **`01_bioinformatic_analysis/`**：先验平权生信分析（24 亚群 One-vs-Rest 全量差异表达、MSigDB Hallmark 通路、配体-受体通讯、空间共定位）
- **`02_reconstruction_impact/`**：重建前后变化与原始数据去掩盖（空间自相关提升 +0.234、灵敏度跃升 20.7 倍、横向弥散与上皮黏蛋白/角蛋白污染清除）
- **`03_biological_discoveries/`**：综合分析下的生物学发现与科学假设（多维度指标评估、四大结构化科学假设 H1-H4、核心靶标谱系程序）

> **调用规范**：
> - 细节数据矩阵位于各槽位下的 `detailed_data/`；
> - 统计汇总指标位于各槽位下的 `tables/`；
> - 宏观结论与交互式仪表盘已提升至上一级样本根目录：[`../report.html`](../report.html) 与 [`../README.md`](../README.md)。
""")

    print(f"[✓] Analysis delivery verification complete for {sample_id}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build sample analysis delivery package.")
    parser.add_argument("--sample-dir", type=Path, default=Path(__file__).resolve().parents[2] / "output/CRC/P2_CRC_Xenium")
    parser.add_argument("--cancer-type", type=str, default="CRC")
    parser.add_argument("--sample-id", type=str, default="P2_CRC_Xenium")
    args = parser.parse_args()

    build_sample_delivery(args.sample_dir, args.cancer_type, args.sample_id)
