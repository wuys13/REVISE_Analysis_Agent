"""Five-Layer Quantitative Evaluation Extractor for REVISE Analysis Agent.

Extracts standard metrics across ST Unit, Gene, Pathway, CCI, and Region layers.
"""
from __future__ import annotations

from typing import Any, Dict, List
import numpy as np
import pandas as pd

class FiveLayerEvaluator:
    """Computes or compiles standard 5-layer reconstruction impact metrics."""

    @staticmethod
    def summarize_layer1_st_unit(subtypes_dict: Dict[str, List[str]]) -> Dict[str, Any]:
        total_clusters = sum(len(v) for v in subtypes_dict.values())
        return {
            "resolved_subtypes_count": total_clusters,
            "lineage_breakdown": {k: len(v) for k, v in subtypes_dict.items()},
            "status": "pass",
        }

    @staticmethod
    def summarize_layer2_gene(raw_degs: int, svc_degs: int, moran_same: float, moran_diff: float) -> Dict[str, Any]:
        fold = round(svc_degs / max(raw_degs, 1), 2)
        return {
            "raw_deg_count": raw_degs,
            "svc_deg_count": svc_degs,
            "expansion_fold": fold,
            "moran_same": moran_same,
            "moran_diff": moran_diff,
            "cross_boundary_leakage_status": "acceptable" if moran_diff <= 0.02 else "warning",
        }

    @staticmethod
    def summarize_layer3_pathway(pathway_scores: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        return {
            "evaluated_pathways_count": len(pathway_scores),
            "top_activated_pathways": pathway_scores,
        }

    @staticmethod
    def summarize_layer4_cci(active_pairs: List[Dict[str, Any]]) -> Dict[str, Any]:
        return {
            "active_pairs_count": len(active_pairs),
            "top_interactions": active_pairs[:5],
        }

    @staticmethod
    def summarize_layer5_region(niches_status: Dict[str, bool]) -> Dict[str, Any]:
        return {
            "niches_detected": niches_status,
            "microenvironment_intact": all(niches_status.values()),
        }
