from __future__ import annotations

import ast
import inspect
import logging
import sys
from types import ModuleType, SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from revise_analysis.methods import biological_metrics
from revise_analysis.methods import cell_communication
from revise_analysis.methods import differential_expression
from revise_analysis.methods import enrichment
from revise_analysis.methods import gene_set_scoring
from revise_analysis.methods import metrics
from revise_analysis.methods import spatial_autocorrelation
from revise_analysis.methods import trajectory
from revise_analysis.methods import unsupervised


METHODS_DIR = (
    __import__("pathlib").Path(__file__).resolve().parents[1]
    / "revise_analysis"
    / "methods"
)


def _adata(
    x: np.ndarray,
    *,
    labels: list[str] | None = None,
    genes: list[str] | None = None,
) -> ad.AnnData:
    n_obs, n_vars = x.shape
    return ad.AnnData(
        x,
        obs=pd.DataFrame(
            {"group": labels or ["A"] * n_obs},
            index=[f"cell-{i}" for i in range(n_obs)],
        ),
        var=pd.DataFrame(index=genes or [f"G{i}" for i in range(n_vars)]),
    )


def test_generic_method_signatures_keep_source_contracts():
    assert list(inspect.signature(differential_expression.get_degs).parameters) == [
        "adata",
        "groupby",
        "method",
        "fc_threshold",
        "reference",
    ]
    assert list(inspect.signature(gene_set_scoring.score_genes).parameters) == [
        "adata",
        "genes",
        "score_name",
    ]
    assert list(inspect.signature(unsupervised.compute_leiden_sweep).parameters) == [
        "adata",
        "resolutions",
        "reference_col",
        "random_state",
    ]
    assert list(
        inspect.signature(
            spatial_autocorrelation.compute_groupwise_spatial_autocorrelation
        ).parameters
    ) == ["adata", "groupby", "mode", "min_obs"]
    assert list(inspect.signature(enrichment.get_enrichment).parameters) == [
        "deg_genes",
        "geneset_file",
        "cutoff",
    ]
    assert list(inspect.signature(enrichment.get_enrichment_local).parameters) == [
        "deg_genes",
        "geneset_file",
        "cutoff",
        "background",
    ]
    assert list(inspect.signature(cell_communication.run_cellphonedb_v5).parameters) == [
        "adata",
        "database_path",
        "celltype_key",
        "min_cell_fraction",
        "min_genes",
        "min_cells",
        "iterations",
        "threshold",
        "pvalue",
        "threads",
        "output_dir",
        "cleanup_temp",
    ]
    assert list(inspect.signature(trajectory.infer_palantir).parameters) == [
        "adata",
        "groupby",
        "basis",
        "use_rep",
        "n_comps",
        "origin_cells",
        "terminal_cells",
        "num_waypoints",
    ]


def test_optional_providers_are_not_imported_at_module_scope():
    optional_names = {"squidpy", "omicverse", "gseapy", "cellphonedb"}
    for path in METHODS_DIR.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.Import):
                imported = {alias.name.split(".")[0] for alias in node.names}
                assert not optional_names.intersection(imported), path.name
            if isinstance(node, ast.ImportFrom):
                assert (node.module or "").split(".")[0] not in optional_names, path.name


def test_differential_expression_preserves_tidy_result_and_filters(monkeypatch):
    calls = []

    def normalize_total(work, **kwargs):
        calls.append(("normalize_total", work, kwargs))

    def log1p(work):
        calls.append(("log1p", work))

    def rank_genes_groups(work, **kwargs):
        calls.append(("rank_genes_groups", work, kwargs))
        work.uns["rank_genes_groups"] = {
            "names": {"A": np.array(["G1", "G2"]), "reference": np.array(["G2", "G1"])},
            "logfoldchanges": {"A": np.array([2.0, -2.0]), "reference": np.array([1.0, 1.0])},
            "pvals": {"A": np.array([0.01, 0.02]), "reference": np.array([0.1, 0.1])},
            "pvals_adj": {"A": np.array([0.02, 0.04]), "reference": np.array([0.2, 0.2])},
        }

    fake_scanpy = SimpleNamespace(
        pp=SimpleNamespace(normalize_total=normalize_total, log1p=log1p),
        tl=SimpleNamespace(rank_genes_groups=rank_genes_groups),
    )
    monkeypatch.setattr(differential_expression, "sc", fake_scanpy)
    data = _adata(
        np.ones((3, 2)),
        labels=["A", "reference", "A"],
        genes=["G1", "G2"],
    )
    result = differential_expression.get_degs(
        data, "group", fc_threshold=1, reference="reference"
    )

    assert [call[0] for call in calls] == [
        "normalize_total",
        "log1p",
        "rank_genes_groups",
    ]
    assert list(result.columns) == [
        "group",
        "gene",
        "logfoldchanges",
        "pvals",
        "pvals_adj",
        "log_q",
    ]
    assert result["gene"].tolist() == ["G1"]
    assert result.index.tolist() == [0]


def test_gmt_reader_and_gene_scoring_preserve_overlap_and_copy(monkeypatch, tmp_path):
    path = tmp_path / "sets.gmt"
    path.write_text("SET_B\tdescription\tG2\tG1\n", encoding="utf-8")
    assert gene_set_scoring.read_gmt(path) == {"SET_B": ["G2", "G1"]}

    malformed = tmp_path / "malformed.gmt"
    malformed.write_text("SET\tdescription\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Malformed GMT row 1"):
        gene_set_scoring.read_gmt(malformed)

    def score_genes(work, *, gene_list, score_name, use_raw):
        assert gene_list == ["G1"]
        assert use_raw is False
        work.obs[score_name] = np.arange(work.n_obs, dtype=float)

    monkeypatch.setattr(
        gene_set_scoring,
        "sc",
        SimpleNamespace(tl=SimpleNamespace(score_genes=score_genes)),
    )
    data = _adata(np.ones((3, 2)), genes=["G1", "G2"])
    scored, key = gene_set_scoring.score_genes(
        data, ["MISSING", "G1"], score_name="SIGNATURE"
    )
    assert key == "SIGNATURE"
    assert key in scored.obs
    assert key not in data.obs


def test_leiden_sweep_preserves_pipeline_schema_and_metric_authority(monkeypatch):
    def filter_cells(_work, **_kwargs):
        return None

    def filter_genes(_work, **_kwargs):
        return None

    def normalize_total(_work, **_kwargs):
        return None

    def log1p(_work):
        return None

    def pca(work, **_kwargs):
        work.obsm["X_pca"] = np.ones((work.n_obs, 1))

    def neighbors(work, **_kwargs):
        work.obsp["connectivities"] = sparse.eye(work.n_obs, format="csr")
        work.obsp["distances"] = sparse.eye(work.n_obs, format="csr")

    def leiden(work, *, key_added, **_kwargs):
        work.obs[key_added] = pd.Categorical(["0", "1", "0", "1"])

    monkeypatch.setattr(
        unsupervised,
        "sc",
        SimpleNamespace(
            pp=SimpleNamespace(
                filter_cells=filter_cells,
                filter_genes=filter_genes,
                normalize_total=normalize_total,
                log1p=log1p,
                neighbors=neighbors,
            ),
            tl=SimpleNamespace(pca=pca, leiden=leiden),
        ),
    )
    authority_calls = []

    def fake_metrics(work, pred_key, true_key):
        authority_calls.append((work, pred_key, true_key))
        return 0.25, 0.75

    monkeypatch.setattr(unsupervised, "compute_clustering_metrics", fake_metrics)
    data = _adata(
        np.ones((4, 3)),
        labels=["A", "A", "B", "B"],
        genes=["MT-DROP", "G1", "G2"],
    )
    processed, result = unsupervised.compute_leiden_sweep(
        data, resolutions=[0.3], reference_col="group"
    )

    assert "MT-DROP" in data.var_names
    assert "MT-DROP" not in processed.var_names
    assert "X_pca" in processed.obsm
    assert {"connectivities", "distances"} <= set(processed.obsp)
    assert list(result.columns) == ["resolution", "ARI", "NMI", "cluster_num"]
    assert result.loc[0, ["ARI", "NMI"]].tolist() == [0.25, 0.75]
    assert authority_calls == [(processed, "leiden_res_0.3", "group")]


def test_groupwise_spatial_autocorrelation_loads_squidpy_lazily(monkeypatch):
    def normalize_total(_work, **_kwargs):
        return None

    def log1p(_work):
        return None

    monkeypatch.setattr(
        spatial_autocorrelation,
        "sc",
        SimpleNamespace(
            pp=SimpleNamespace(normalize_total=normalize_total, log1p=log1p)
        ),
    )
    fake_squidpy = ModuleType("squidpy")
    graph_calls = []

    def spatial_neighbors(work):
        graph_calls.append(tuple(work.obs_names))

    def spatial_autocorr(work, *, mode, genes, n_perms, n_jobs):
        statistic = "I" if mode == "moran" else "C"
        key = "moranI" if mode == "moran" else "gearyC"
        work.uns[key] = pd.DataFrame(
            {statistic: np.arange(work.n_vars, dtype=float) + work.n_obs},
            index=work.var_names,
        )

    fake_squidpy.gr = SimpleNamespace(
        spatial_neighbors=spatial_neighbors,
        spatial_autocorr=spatial_autocorr,
    )
    monkeypatch.setitem(sys.modules, "squidpy", fake_squidpy)
    data = _adata(
        np.arange(16, dtype=float).reshape(4, 4) + 1,
        labels=["A", "A", "B", "B"],
    )
    data.obsm["spatial"] = np.arange(8, dtype=float).reshape(4, 2)
    result = spatial_autocorrelation.compute_groupwise_spatial_autocorrelation(
        data, groupby="group", min_obs=2
    )

    assert list(result.columns) == ["All", "A", "B"]
    assert graph_calls == [
        tuple(data.obs_names),
        tuple(data.obs_names[:2]),
        tuple(data.obs_names[2:]),
    ]


def test_enrichment_keeps_positional_cutoff_and_provider_schema(monkeypatch):
    expected = pd.DataFrame({"Term": ["P1"], "Adjusted P-value": [0.01]})
    calls = []

    def enrichr(**kwargs):
        calls.append(kwargs)
        return SimpleNamespace(results=expected)

    monkeypatch.setattr(
        enrichment,
        "_require_gseapy",
        lambda: SimpleNamespace(enrichr=enrichr),
    )
    result = enrichment.get_enrichment(["G1"], "genes.gmt", 0.1)
    assert result is expected
    assert calls == [
        {
            "gene_list": ["G1"],
            "gene_sets": "genes.gmt",
            "organism": "human",
            "cutoff": 0.1,
        }
    ]

    empty = enrichment.get_enrichment([], "genes.gmt")
    assert empty.empty
    assert empty.columns.tolist() == list(enrichment.ENRICHMENT_RESULT_COLUMNS)


def test_cell_communication_restores_table_validation_and_forwarding(monkeypatch, tmp_path):
    table = pd.DataFrame(
        {
            "gene_a": ["L1"],
            "gene_b": ["R1"],
            "A|A": [1.0],
            "A|B": [2.0],
        }
    )
    results = {
        "pvalues": table.copy(),
        "significant_means": table.copy(),
        "interaction_scores": table.copy(),
    }
    normalized = cell_communication.normalize_cellphonedb_tables(results, ["A", "B"])
    assert normalized["pvalues"].index.tolist() == ["L1_R1"]
    assert normalized["pvalues"].index.name == "LR_pair"

    calls = []

    def run(_adata, **kwargs):
        calls.append(kwargs)
        return results, "returned-adata"

    monkeypatch.setattr(
        cell_communication,
        "_require_cellphonedb_dependencies",
        lambda: SimpleNamespace(single=SimpleNamespace(run_cellphonedb_v5=run)),
    )
    database = tmp_path / "cellphonedb.zip"
    database.write_bytes(b"db")
    data = _adata(np.ones((2, 2)), labels=["A", "B"])
    actual_results, actual_adata = cell_communication.run_cellphonedb_v5(
        data,
        database_path=database,
        celltype_key="group",
        min_cell_fraction=0.005,
        min_genes=200,
        min_cells=3,
        iterations=1000,
        threshold=0.1,
        pvalue=0.05,
        threads=2,
        output_dir=tmp_path / "out",
        cleanup_temp=True,
    )
    assert actual_results is results
    assert actual_adata == "returned-adata"
    assert calls[0]["cpdb_file_path"] == str(database)
    assert calls[0]["celltype_key"] == "group"


def test_trajectory_restores_validation_and_pseudotime_contract(monkeypatch):
    calls = []

    class FakeTrajectory:
        def __init__(self, data):
            self.data = data

        def set_origin_cells(self, value):
            calls.append(("origin", value))

        def set_terminal_cells(self, value):
            calls.append(("terminal", value))

        def inference(self, **kwargs):
            calls.append(("inference", kwargs))
            self.data.obs["palantir_pseudotime"] = [0.0, 1.0]

    def infer(data, **kwargs):
        calls.append(("TrajInfer", kwargs))
        return FakeTrajectory(data)

    fake_omicverse = ModuleType("omicverse")
    fake_omicverse.single = SimpleNamespace(TrajInfer=infer)
    monkeypatch.setitem(sys.modules, "omicverse", fake_omicverse)
    data = _adata(
        np.ones((2, 2)),
        labels=["SMC", "Fibroblast"],
    )
    data.obs.rename(columns={"group": "Level1"}, inplace=True)
    result = trajectory.infer_palantir(
        data,
        groupby="Level1",
        basis="X_umap",
        use_rep="X_pca",
        n_comps=2,
        origin_cells="SMC",
        terminal_cells=["Fibroblast"],
        num_waypoints=20,
    )

    assert calls == [
        (
            "TrajInfer",
            {
                "basis": "X_umap",
                "groupby": "Level1",
                "use_rep": "X_pca",
                "n_comps": 2,
            },
        ),
        ("origin", "SMC"),
        ("terminal", ["Fibroblast"]),
        ("inference", {"method": "palantir", "num_waypoints": 20}),
    ]
    assert result.tolist() == [0.0, 1.0]


def test_biological_metrics_keep_moran_marker_and_identity_contracts():
    data = _adata(
        np.array(
            [
                [4.0, 1.0, 0.0],
                [3.0, 1.0, 0.0],
                [0.0, 4.0, 1.0],
                [0.0, 3.0, 1.0],
            ]
        ),
        labels=["A", "A", "B", "B"],
        genes=["G1", "G2", "G3"],
    )
    data.obsp["spatial_connectivities"] = sparse.csr_matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]],
        dtype=float,
    )
    moran = biological_metrics.compute_conditional_moran_i(
        data,
        cell_type_col="group",
        normalize=False,
        log1p=False,
    )
    assert list(moran.columns) == [
        "Gene",
        "MISC",
        "MIDC",
        "n_same_edges",
        "n_diff_edges",
    ]

    tmp, summary = biological_metrics.compute_tmp_mer(
        data,
        "group",
        marker_map={"A": ["G1"], "B": ["G2"]},
    )
    assert set(tmp["metric_status"]) == {"ok"}
    assert summary.loc[0, "n_targets_ok"] == 2

    data.obsm["X_pca"] = np.array([[0, 0], [0, 1], [4, 4], [4, 5]], dtype=float)
    identity = biological_metrics.compute_identity_metrics(data, "group")
    assert identity.loc[0, "metric"] == "ASW"
    assert np.isfinite(identity.loc[0, "value"])


def test_metrics_keep_compatibility_columns_and_label_scores(monkeypatch):
    data = _adata(
        np.arange(16, dtype=float).reshape(8, 2) + 1,
        labels=["A", "A", "A", "A", "B", "B", "B", "B"],
        genes=["G1", "G2"],
    )
    predicted = _adata(
        np.arange(16, dtype=float).reshape(8, 2) + 1.5,
        labels=["A", "A", "A", "A", "B", "B", "B", "B"],
        genes=["G1", "G2"],
    )
    monkeypatch.setattr(
        metrics,
        "sc",
        SimpleNamespace(
            pp=SimpleNamespace(
                normalize_total=lambda *_args, **_kwargs: None,
                log1p=lambda *_args, **_kwargs: None,
            )
        ),
    )
    result = metrics.compute_metric(
        data,
        predicted,
        logging.getLogger("generic-parity"),
    )
    assert list(result.columns) == ["Gene", "PCC", "SSIM", "MSE", "NRMSE"]
    assert result["Gene"].tolist() == ["G1", "G2"]
    assert metrics.compute_clustering_metrics(
        data, "group", "group"
    ) == (1.0, 1.0)

