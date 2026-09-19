from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from anndata import AnnData

from revise_analysis.io import Sample, load_sample


def _objects(values=None):
    matrix = np.asarray(values if values is not None else [[0.125, 2.5], [3.75, 0.0]])
    obs = pd.DataFrame(index=["u1", "u2"])
    var = pd.DataFrame(index=["g1", "g2"])
    raw = AnnData(matrix.copy(), obs=obs.copy(), var=var.copy())
    return raw, raw.copy()


def _write_sample(tmp_path, *, expression, values=None):
    raw, svc = _objects(values)
    raw.write_h5ad(tmp_path / "raw.h5ad")
    svc.write_h5ad(tmp_path / "SVC.h5ad")
    source = tmp_path / "sample.yaml"
    source.write_text(yaml.safe_dump({
        "schema_version": 1,
        "sample_id": "linear",
        "files": {"raw": "raw.h5ad", "svc": "SVC.h5ad"},
        "expression": expression,
    }))
    return source


def test_known_identity_without_scale_uses_fixed_linear_contract_and_preserves_decimals(tmp_path):
    source = _write_sample(tmp_path, expression={
        side: {"matrix": "X", "identity": "measured_expression"}
        for side in ("raw", "svc")
    })
    sample = load_sample(source)
    assert sample.expression("raw") == {
        "identity": "measured_expression",
        "scale": "untransformed_nonnegative",
        "matrix": "X",
    }
    assert sample.expression_unavailable("raw") is None
    np.testing.assert_array_equal(sample.raw.X, [[0.125, 2.5], [3.75, 0.0]])


def test_identity_and_legacy_unknown_scale_remain_unavailable(tmp_path):
    source = _write_sample(tmp_path, expression={
        "raw": {"identity": "unknown"},
        "svc": {"identity": "reconstructed_expression", "scale": "unknown"},
    })
    sample = load_sample(source)
    assert "identity is unknown" in sample.expression_unavailable("raw")
    assert "explicitly unknown" in sample.expression_unavailable("svc")
    assert sample.expression("svc")["scale"] == "untransformed_nonnegative"


@pytest.mark.parametrize("scale", ["log", "log1p", "log1p_nonnegative"])
def test_legacy_log_declarations_are_rejected_by_loader_and_direct_sample(tmp_path, scale):
    expression = {side: {"identity": "known", "scale": scale} for side in ("raw", "svc")}
    source = _write_sample(tmp_path, expression=expression)
    with pytest.raises(ValueError, match="linear X contract"):
        load_sample(source)
    raw, svc = _objects()
    with pytest.raises(ValueError, match="linear X contract"):
        Sample("direct", raw, svc, {"expression": expression}, Path("sample.yaml"))


def test_unknown_historical_matrix_remains_loadable_for_label_space(tmp_path):
    source = _write_sample(
        tmp_path,
        expression={side: {"identity": "unknown"} for side in ("raw", "svc")},
        values=[[1.0, -0.01], [np.nan, 3.0]],
    )
    sample = load_sample(source)
    assert sample.raw.n_obs == 2
    assert "identity is unknown" in sample.expression_unavailable("raw")


def test_legacy_untransformed_alias_remains_compatible(tmp_path):
    source = _write_sample(tmp_path, expression={
        side: {"identity": "known", "scale": "untransformed"}
        for side in ("raw", "svc")
    })
    sample = load_sample(source)
    assert sample.expression_unavailable("raw") is None
    assert sample.expression("raw")["scale"] == "untransformed_nonnegative"


@pytest.mark.parametrize("shared, per_side", [
    ("log1p", "untransformed_nonnegative"),
    ("untransformed_nonnegative", "log1p"),
])
def test_per_side_scale_cannot_hide_an_incompatible_shared_or_side_declaration(tmp_path, shared, per_side):
    source = _write_sample(tmp_path, expression={
        "scale": shared,
        "raw": {"identity": "known", "scale": per_side},
        "svc": {"identity": "known"},
    })
    with pytest.raises(ValueError, match="linear X contract"):
        load_sample(source)


def test_per_side_linear_scale_cannot_hide_shared_unknown(tmp_path):
    source = _write_sample(tmp_path, expression={
        "scale": "unknown",
        "raw": {"identity": "known", "scale": "untransformed_nonnegative"},
        "svc": {"identity": "known"},
    })
    sample = load_sample(source)
    assert "explicitly unknown" in sample.expression_unavailable("raw")
