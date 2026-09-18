"""Reusable analysis methods migrated from REVISE.

Source provenance: selected scientific helpers from REVISE main@c83dc97 and
reconstruction-impact@e82dd13 (MIT licensed).  Optional providers remain lazy.
"""

from .partition import PartitionComparison, compare_membership, compare_partitions, compute_partition
from .programs import compute_pathway_scores
from .spatial import compute_moran

__all__ = [
    "PartitionComparison", "compare_membership", "compare_partitions",
    "compute_partition", "compute_moran", "compute_pathway_scores",
]
