"""AUCell pathway scoring, intentionally one supplied gene set at a time."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import pandas as pd

from ._shared import analysis_parameters, deterministic_subset, save_json, save_table, unavailable


def compute_pathway_scores(adata: Any, *, genes: list[str], score_name: str, auc_threshold: float, seed: int):
    """Expose the provider-backed AUCell calculation for notebooks."""
    from revise_analysis.methods.programs import compute_pathway_scores as implementation
    return implementation(adata, genes=genes, score_name=score_name, auc_threshold=auc_threshold, seed=seed)


def run(sample: Any, output_dir: Path, parameters: dict | None = None) -> dict:
    supplied = analysis_parameters(parameters, {"gene_sets", "geneset_path", "gene_set_names", "pathway_auc_threshold", "sample_n_units", "random_state"})
    gene_sets, resource_error = resolve_gene_sets(supplied)
    effective = {
        "gene_sets": gene_sets or {}, "geneset_path": supplied.get("geneset_path"), "gene_set_names": supplied.get("gene_set_names"),
        "pathway_auc_threshold": float(supplied.get("pathway_auc_threshold", 0.05)),
        "sample_n_units": supplied.get("sample_n_units"),
        "random_state": int(supplied.get("random_state", 42)),
    }
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs, missing, metadata = {}, [], {}
    if resource_error:
        unavailable("pathway_resource", resource_error, missing)
        return {"status": "skipped", "parameters": effective, "outputs": outputs, "unavailable": missing}
    for side in ("raw", "svc"):
        adata = deterministic_subset(getattr(sample, side), effective["sample_n_units"], effective["random_state"] + (side == "svc"))
        for name, genes in gene_sets.items():
            component = f"pathway:{side}:{name}"
            try:
                result = compute_pathway_scores(
                    adata, genes=genes, score_name=name,
                    auc_threshold=effective["pathway_auc_threshold"], seed=effective["random_state"],
                )
                score = result["scores"].rename(name)
                save_table(score, output_dir, f"tables/pathway_{side}_{_slug(name)}.csv", outputs)
                metadata[component] = {"coverage": result["coverage"], "metadata": result["metadata"]}
            except ImportError as exc:
                # Absence of the optional AUCell provider leaves other pathways/sides usable.
                unavailable(component, str(exc), missing)
            except ValueError as exc:
                unavailable(component, f"skipped: {exc}", missing)
    if metadata:
        save_json(metadata, output_dir, "tables/pathway_metadata.json", outputs)
    _write_figures(output_dir, outputs, missing)
    status = "partial" if outputs and missing else ("succeeded" if outputs else "skipped")
    return {"status": status, "parameters": effective, "outputs": outputs, "unavailable": missing}


def _slug(value: str) -> str:
    return "".join(char if char.isalnum() else "_" for char in value).strip("_") or "pathway"


def resolve_gene_sets(parameters: dict) -> tuple[dict | None, str | None]:
    """Resolve explicit sets or a runner-resolved GMT resource within this analysis."""
    explicit = parameters.get("gene_sets")
    if explicit is not None and parameters.get("geneset_path") is not None:
        return None, "use either gene_sets or geneset_path, not both"
    if explicit is not None and parameters.get("gene_set_names") is not None:
        return None, "gene_set_names requires geneset_path"
    if explicit is not None:
        if not isinstance(explicit, dict) or not explicit:
            return None, "gene_sets must be a non-empty mapping of pathway name to gene IDs"
        if not all(isinstance(name, str) and isinstance(genes, list) and genes for name, genes in explicit.items()):
            return None, "gene_sets values must be non-empty lists of gene IDs"
        if len({_slug(name) for name in explicit}) != len(explicit):
            raise ValueError("gene-set names collide after conversion to artifact filenames")
        return explicit, None
    path = parameters.get("geneset_path")
    if not path:
        return None, "no gene sets supplied"
    try:
        from revise_analysis.io import read_gene_sets
        selected = read_gene_sets(path, names=parameters.get("gene_set_names"))
        if not selected:
            return None, "gene set resource selected no non-empty gene sets"
        if len({_slug(name) for name in selected}) != len(selected):
            raise ValueError("gene-set names collide after conversion to artifact filenames")
        return selected, None
    except (FileNotFoundError, ValueError, OSError) as exc:
        return None, f"gene set resource unavailable: {exc}"


def _write_figures(output_dir: Path, outputs: dict, missing: list[dict]) -> None:
    try:
        from revise_analysis.plotting import plot_pathway_scores
        for table in (output_dir / "tables").glob("pathway_*.csv"):
            relative = f"figures/{table.stem}.png"
            score = pd.read_csv(table, index_col=0).iloc[:, 0]
            plot_pathway_scores(score, output_dir / relative, title=table.stem.replace("_", " "))
            outputs[relative] = relative
    except ImportError as exc:
        unavailable("figures", str(exc), missing)
