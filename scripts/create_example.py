"""Create explicitly synthetic, linear-count inputs for local smoke testing."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy import sparse
from anndata import AnnData
import yaml


def create_example(root: Path):
    root = root.resolve()
    destination = root / "data" / "example"
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("raw.h5ad", "SVC.h5ad", "sample.yaml"):
        if (destination / name).exists():
            raise FileExistsError(f"Example already exists; choose another --root: {destination / name}")
    generator = np.random.default_rng(42)
    genes = [f"G{i}" for i in range(600)]
    broad = np.array(["Fibroblast", "Mono/Macro", "T", "Tumor", "Intestinal Epithelial", "B"])
    for side, n_units in (("raw", 240), ("svc", 210)):
        ids = [f"u{i:04}" for i in range(n_units)] if side == "raw" else [f"u{i:04}" for i in range(120)] + [f"new{i:04}" for i in range(90)]
        labels = broad[np.arange(n_units) % len(broad)]
        if side == "svc":
            labels = np.where(labels == "Mono/Macro", "Mono_Macro", labels)
        means = np.full((n_units, len(genes)), 0.25)
        for i in range(n_units):
            start = (i % 6) * 60
            means[i, start:start + 60] += 3.0
        matrix = generator.poisson(means).astype(np.float32)
        if side == "svc":
            matrix = matrix * 0.8  # nonnegative linear synthetic reconstruction
        x = (np.arange(n_units) % 20) * 8.0
        y = (np.arange(n_units) // 20) * 8.0
        adata = AnnData(sparse.csr_matrix(matrix), obs=pd.DataFrame({"Level1": labels}, index=ids), var=pd.DataFrame(index=genes))
        adata.obsm["spatial"] = np.column_stack([x, y])
        adata.uns["example"] = "synthetic linear expression; not biological evidence"
        adata.write_h5ad(destination / ("raw.h5ad" if side == "raw" else "SVC.h5ad"))
    sample = {
        "schema_version": 1, "sample_id": "example",
        "files": {"raw": "raw.h5ad", "svc": "SVC.h5ad"},
        "expression": {"scale": "untransformed_nonnegative"},
        "columns": {"broad": "Level1", "subtype": "Level2"},
        "label_aliases": {"Mono_Macro": "Mono/Macro"},
        "spatial": {"key": "spatial", "unit": "micron", "microns_per_coordinate": 1.0},
    }
    (destination / "sample.yaml").write_text(yaml.safe_dump(sample, sort_keys=False), encoding="utf-8")
    configs = root / "configs"
    configs.mkdir(exist_ok=True)
    project = {
        "schema_version": 1, "output_dir": "../output", "samples": ["../data/example/sample.yaml"],
        "analyses": {
            "reconstruction_impact": {"membership_same_units": True, "sample_n_units": 180, "n_top_genes": 100, "window_side_microns": 32.0,
                                       "gene_sets": {"DEMO_PROGRAM": genes[:50]}, "pathway_auc_threshold": 0.05},
            "spatial_autocorrelation": {"sample_n_units": 180, "moran_n_neighbors": 6},
            "pathway_activity": {"gene_sets": {"DEMO_PROGRAM": genes[:50]}, "pathway_auc_threshold": 0.05, "sample_n_units": 180},
        },
    }
    project_path = configs / "example_project.yaml"
    if project_path.exists():
        raise FileExistsError(project_path)
    project_path.write_text(yaml.safe_dump(project, sort_keys=False), encoding="utf-8")
    print(f"Created synthetic example: {destination}")
    print(f"Run: revise-analysis batch --config {project_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    create_example(parser.parse_args().root)
