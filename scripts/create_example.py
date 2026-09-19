"""Create synthetic reconstruction inputs for local smoke testing."""
from pathlib import Path
import argparse
import tempfile
import numpy as np
import pandas as pd
from scipy import sparse
from anndata import AnnData
import yaml


def _safe_sample_id(value: str) -> str:
    if (not isinstance(value, str) or not value.strip() or value in {".", ".."}
            or value.startswith(".") or any(token in value for token in ("/", "\\", "\0"))):
        raise ValueError(f"Unsafe sample ID: {value!r}")
    return value


def _project(sample_id: str, genes: list[str]) -> dict:
    return {
        "schema_version": 1,
        "output_dir": "../output",
        "samples": [f"../data/{sample_id}/sample.yaml"],
        "analyses": {
            "reconstruction_impact": {
                "membership_same_units": True,
                "sample_n_units": 180,
                "n_top_genes": 100,
                "window_side_microns": 32.0,
                "anatomy_window_side_microns": 40.0,
                "min_window_units": 4,
                "n_window_draws": 50,
                "region_n_bootstrap": 100,
                "gene_sets": {"DEMO_PROGRAM": genes[:50]},
                "pathway_auc_threshold": 0.05,
            },
        },
    }


def create_example(root: Path, *, sample_id: str = "example_reconstruction"):
    root = root.resolve()
    sample_id = _safe_sample_id(sample_id)
    data_root = root / "data"
    configs = root / "configs"
    destination = data_root / sample_id
    project_path = configs / "example_project.yaml"
    # Check every final target before producing either H5AD. An existing
    # project config is user/manual content: keep it and generate only the new
    # sample directory. An existing sample directory is never merged into or
    # overwritten.
    if destination.exists():
        raise FileExistsError(f"Example already exists; choose another --sample-id: {destination}")
    if project_path.exists() and not project_path.is_file():
        raise FileExistsError(f"Project config path is not a regular file: {project_path}")
    data_root.mkdir(parents=True, exist_ok=True)
    configs.mkdir(parents=True, exist_ok=True)
    generator = np.random.default_rng(42)
    genes = [f"G{i}" for i in range(600)]
    broad = np.array(["Fibroblast", "Mono/Macro", "T", "Tumor", "Intestinal Epithelial", "B"])
    with tempfile.TemporaryDirectory(prefix=f".{sample_id}-", dir=data_root) as temporary:
        staging = Path(temporary)
        for side, n_units in (("raw", 240), ("svc", 210)):
            ids = ([f"u{i:04}" for i in range(n_units)] if side == "raw"
                   else [f"u{i:04}" for i in range(120)] + [f"new{i:04}" for i in range(90)])
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
            adata = AnnData(
                sparse.csr_matrix(matrix),
                obs=pd.DataFrame({"Level1": labels}, index=ids),
                var=pd.DataFrame(index=genes),
            )
            if side == "svc":
                adata.obs["SVC_cluster"] = [f"synthetic_{i % 4}" for i in range(n_units)]
            adata.obsm["spatial"] = np.column_stack([x, y])
            adata.uns["example"] = "synthetic linear expression; not biological evidence"
            adata.write_h5ad(staging / ("raw.h5ad" if side == "raw" else "SVC.h5ad"))
        sample = {
            "schema_version": 1,
            "sample_id": sample_id,
            "files": {"raw": "raw.h5ad", "svc": "SVC.h5ad"},
            "expression": {
                side: {"matrix": "X", "identity": "synthetic_expression", "scale": "untransformed_nonnegative"}
                for side in ("raw", "svc")
            },
            "columns": {"broad": "Level1", "subtype": "Level2", "reconstruction": "SVC_cluster"},
            "label_aliases": {"Mono_Macro": "Mono/Macro"},
            "spatial": {"key": "spatial", "unit": "micron", "microns_per_coordinate": 1.0},
        }
        (staging / "sample.yaml").write_text(yaml.safe_dump(sample, sort_keys=False), encoding="utf-8")

        created_project = False
        if not project_path.exists():
            # Exclusive creation protects a config that appears after preflight.
            with project_path.open("x", encoding="utf-8") as stream:
                yaml.safe_dump(_project(sample_id, genes), stream, sort_keys=False)
            created_project = True
        try:
            staging.rename(destination)
        except Exception:
            if created_project:
                project_path.unlink(missing_ok=True)
            raise
    print(f"Created synthetic example: {destination}")
    if created_project:
        print(f"Created project config: {project_path}")
    else:
        print(f"Preserved existing project config: {project_path}")
    print(f"Run this sample: revise-analysis run --sample {destination / 'sample.yaml'} --analysis reconstruction_impact --output {root / 'output'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--sample-id", default="example_reconstruction")
    arguments = parser.parse_args()
    create_example(arguments.root, sample_id=arguments.sample_id)
