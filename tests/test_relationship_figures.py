from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from revise_analysis.plotting.relationships import render_relationship_figures


def _save(root: Path, relative: str, frame: pd.DataFrame) -> str:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    return relative


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_relationship_figures_read_only_registered_contained_artifacts(tmp_path, monkeypatch):
    root = tmp_path / "result"
    registered = _save(
        root,
        "tables/integrated_changed_units_All.csv",
        pd.DataFrame({
            "state_region": [True, False],
            "anatomy_region": ["Tumor", "Normal"],
            "n_compared": [10, 20],
            "n_changed": [3, 4],
            "changed_fraction": [.3, .2],
            "scope": ["All", "All"],
        }),
    )
    stale = _save(
        root,
        "tables/integrated_changed_units_stale.csv",
        pd.DataFrame({"n_compared": [5], "n_changed": [1], "changed_fraction": [.2]}),
    )
    outside = _save(
        tmp_path,
        "outside.csv",
        pd.DataFrame({"n_compared": [5], "n_changed": [1], "changed_fraction": [.2]}),
    )
    before = {path: _digest(root / path) for path in (registered, stale)}

    from revise_analysis.plotting import relationships
    reads = []
    original = relationships.pd.read_csv

    def recording_read(path, *args, **kwargs):
        reads.append(Path(path).resolve())
        return original(path, *args, **kwargs)

    monkeypatch.setattr(relationships.pd, "read_csv", recording_read)
    result = render_relationship_figures(
        root,
        {
            registered: registered,
            "missing": "tables/integrated_changed_units_missing.csv",
            "escape": "../outside.csv",
        },
    )

    expected = "figures/relationships_integrated_changed_units_All.png"
    assert result == {expected: expected}
    assert reads == [(root / registered).resolve()]
    assert not (root / "figures/relationships_integrated_changed_units_stale.png").exists()
    assert {path: _digest(root / path) for path in (registered, stale)} == before


def test_missing_unscored_and_unavailable_relationships_do_not_make_empty_figures(tmp_path):
    root = tmp_path / "result"
    files = [
        _save(
            root,
            "tables/integrated_program_svc_All_P_units.csv",
            pd.DataFrame({
                "x": [0.0, 1.0], "y": [0.0, 1.0], "score": [np.nan, np.nan],
                "side": ["svc", "svc"], "scope": ["All", "All"], "program": ["P", "P"],
            }),
        ),
        _save(
            root,
            "tables/integrated_program_svc_All_P_anatomy.csv",
            pd.DataFrame({
                "side": ["svc"], "scope": ["All"], "program": ["P"],
                "anatomy_region": ["Tumor"], "n_units": [2], "n_scored": [0],
                "score_mean": [np.nan], "score_median": [np.nan], "missing_status": ["unavailable"],
            }),
        ),
        _save(
            root,
            "tables/integrated_region_anatomy_summary_All.csv",
            pd.DataFrame({
                "scope": ["All", "All"], "region_kind": ["gain", "gain"],
                "in_region": [np.nan, np.nan], "anatomy_region": ["Tumor", "Normal"],
                "n_units": [3, 7], "denominator_units": [10, 10],
                "fraction_of_parent_units": [.3, .7],
                "region_definition_status": ["unavailable", "unavailable"],
            }),
        ),
    ]

    result = render_relationship_figures(root, {path: path for path in files})

    assert result == {}
    assert not (root / "figures").exists()


def test_program_summary_keeps_unscored_group_as_missing_when_another_group_is_scored(tmp_path, monkeypatch):
    root = tmp_path / "result"
    relative = _save(
        root,
        "tables/integrated_program_svc_All_P_anatomy.csv",
        pd.DataFrame({
            "side": ["svc", "svc"], "scope": ["All", "All"], "program": ["P", "P"],
            "anatomy_region": ["Tumor", "Normal"], "n_units": [6, 4], "n_scored": [5, 0],
            "score_mean": [.4, np.nan], "score_median": [.3, np.nan],
            "missing_status": ["partial", "unavailable"],
        }),
    )
    from revise_analysis.plotting import relationships
    labels = []
    original_save = relationships._save

    def recording_save(figure, destination):
        labels.extend(text.get_text() for axis in figure.axes for text in axis.texts)
        return original_save(figure, destination)

    monkeypatch.setattr(relationships, "_save", recording_save)
    result = render_relationship_figures(root, {relative: relative})

    expected = "figures/relationships_integrated_program_svc_All_P_anatomy.png"
    assert result == {expected: expected}
    assert "未评分 0/4" in labels


def test_relationship_figures_render_each_saved_table_role(tmp_path):
    root = tmp_path / "result"
    files = []
    files.append(_save(
        root,
        "tables/integrated_region_anatomy_summary_Fibroblast.csv",
        pd.DataFrame({
            "scope": ["Fibroblast"] * 5,
            "region_kind": ["state", "state", "state", "gain", "gain"],
            "in_region": [True, False, np.nan, np.nan, np.nan],
            "anatomy_region": ["Tumor", "Tumor", "Normal", "Tumor", "Normal"],
            "n_units": [2, 3, 5, 4, 6],
            "denominator_units": [10] * 5,
            "fraction_of_parent_units": [.2, .3, .5, .4, .6],
            "region_definition_status": ["available"] * 3 + ["unavailable"] * 2,
        }),
    ))
    files.append(_save(
        root,
        "tables/integrated_label_composition_Fibroblast.csv",
        pd.DataFrame({
            "state_region": [True, True, False, False],
            "anatomy_region": ["Tumor"] * 4,
            "reconstruction_label": ["A", "B", "A", "B"],
            "n_units": [3, 2, 1, 4], "fraction": [.6, .4, .2, .8],
            "denominator_units": [5] * 4, "scope": ["Fibroblast"] * 4,
        }),
    ))
    files.append(_save(
        root,
        "tables/integrated_program_svc_Fibroblast_P_units.csv",
        pd.DataFrame({
            "x": [0.0, 1.0, 2.0], "y": [0.0, 1.0, 2.0], "score": [.1, np.nan, .8],
            "side": ["svc"] * 3, "scope": ["Fibroblast"] * 3, "program": ["P"] * 3,
        }),
    ))
    files.append(_save(
        root,
        "tables/integrated_program_svc_Fibroblast_P_anatomy.csv",
        pd.DataFrame({
            "side": ["svc", "svc"], "scope": ["Fibroblast", "Fibroblast"], "program": ["P", "P"],
            "anatomy_region": ["Tumor", "Normal"], "n_units": [8, 4], "n_scored": [6, 2],
            "score_mean": [.4, .6], "score_median": [.3, .5], "missing_status": ["partial", "partial"],
        }),
    ))
    files.append(_save(
        root,
        "tables/integrated_program_svc_Fibroblast_P_region.csv",
        pd.DataFrame({
            "side": ["svc", "svc"], "scope": ["Fibroblast", "Fibroblast"], "program": ["P", "P"],
            "state_region": [True, False], "gain_region": [np.nan, np.nan],
            "anatomy_region": ["Tumor", "Normal"], "n_units": [8, 4], "n_scored": [6, 2],
            "score_mean": [.4, .6], "score_median": [.3, .5], "missing_status": ["partial", "partial"],
        }),
    ))
    files.append(_save(
        root,
        "tables/integrated_changed_units_Fibroblast.csv",
        pd.DataFrame({
            "state_region": [True, False], "anatomy_region": ["Tumor", "Normal"],
            "n_compared": [8, 4], "n_changed": [2, 1], "changed_fraction": [.25, .25],
            "scope": ["Fibroblast", "Fibroblast"],
        }),
    ))
    files.append(_save(
        root,
        "tables/moran_shared_genes_Fibroblast.csv",
        pd.DataFrame({
            "gene_id": ["a", "b", "c"], "raw_moran": [.1, .2, np.nan], "svc_moran": [.2, .1, .4],
            "comparison_available": [True, True, False],
            "raw_n_units": [40] * 3, "svc_n_units": [35] * 3,
            "raw_n_native_genes": [100] * 3, "svc_n_native_genes": [90] * 3,
            "n_common_genes": [3] * 3, "scope": ["Fibroblast"] * 3,
        }),
    ))
    metadata = root / "tables/pathway_svc_P_metadata.json"
    metadata.write_text(json.dumps({"coverage": {"fraction_present": .75}}), encoding="utf-8")
    files.append("tables/pathway_svc_P_metadata.json")
    before = {path: _digest(root / path) for path in files}

    result = render_relationship_figures(root, {path: path for path in files})

    expected = {
        f"figures/relationships_{Path(path).stem}.png"
        for path in files if path.endswith(".csv")
    }
    assert set(result) == expected
    assert all((root / path).stat().st_size > 0 for path in expected)
    assert {path: _digest(root / path) for path in files} == before
