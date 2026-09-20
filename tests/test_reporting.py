import json
from pathlib import Path

import pandas as pd

from revise_analysis.reporting import render_report
from revise_analysis.reporting.impact_reading import reading_summary
from revise_analysis.reporting import report as report_module


def test_static_report_reads_only_saved_artifacts(tmp_path: Path):
    (tmp_path / "tables").mkdir()
    (tmp_path / "tables/overview.csv").write_text("side,n_units\nraw,7\nsvc,5\n", encoding="utf-8")
    (tmp_path / "result.json").write_text(json.dumps({
        "status": "partial", "parameters": {"random_state": 42},
        "outputs": {"overview": "tables/overview.csv"},
        "unavailable": [{"component": "pathway:raw:EMT", "reason": "OmicVerse unavailable"}],
    }), encoding="utf-8")

    report = render_report(tmp_path)
    text = report.read_text(encoding="utf-8")
    assert report.name == "report.html"
    assert "raw" in text and "OmicVerse unavailable" in text
    assert "AnnData" not in text


def test_rawbaseline_report_appends_level2_coverage_to_partition_fact(tmp_path: Path):
    tables = tmp_path / "tables"
    tables.mkdir()
    partition = tables / "partitions_raw_All.csv"
    partition.write_text("partition\n0\n1\n", encoding="utf-8")
    coverage = tables / "raw_level2_baseline_coverage.csv"
    coverage.write_text(
        "scope,n_total_units,n_valid_units,n_missing_units,status\nAll,12,9,3,partial\n",
        encoding="utf-8",
    )
    artifacts = [
        {"relative": "tables/partitions_raw_All.csv", "path": partition, "suffix": ".csv"},
        {"relative": "tables/raw_level2_baseline_coverage.csv", "path": coverage, "suffix": ".csv"},
    ]
    fact = report_module._scope_fact("rawbaseline", "All", "All", artifacts, {})
    assert "Raw baseline" in fact
    assert "有效 9/12" in fact and "排除 3" in fact


def test_static_report_uses_explicit_sections_and_keeps_failed_stage_evidence(tmp_path: Path):
    (tmp_path / "tables").mkdir()
    (tmp_path / "tables/support_curves.csv").write_text(
        "curve,window_side_microns,retained_fraction\n"
        "svc,16,0.9\nsvc,32,0.9\nraw,16,0.8\nraw,32,0.8\n",
        encoding="utf-8",
    )
    (tmp_path / "tables/state_threshold.json").write_text(
        '{"status": "no_stable_threshold", "threshold": null}\n',
        encoding="utf-8",
    )
    (tmp_path / "result.json").write_text(json.dumps({
        "analysis": "reconstruction_impact",
        "sample_id": "P2CRC_Xenium",
        "status": "partial",
        "parameters": {"random_state": 42, "window_side_microns": 32},
        "outputs": {
            "support_curves": "tables/support_curves.csv",
            "state_threshold": "tables/state_threshold.json",
        },
        "sections": [
            {
                "id": "support",
                "title": "Support and scale",
                "description": "Raw, SVC, and common-valid support remain visible.",
                "artifacts": ["support_curves", "missing_curve.csv"],
            },
            {
                "id": "regions",
                "title": "State and Gain regions",
                "artifacts": ["state_threshold"],
            },
        ],
        "stages": [
            {"stage": "state:Fibroblast", "status": "unavailable", "reason": "no_stable_threshold"},
        ],
        "stage_errors": [
            {"stage": "molecular:AUCell", "reason": "provider unavailable"},
        ],
        "unavailable": [
            {"component": "state:Fibroblast", "reason": "threshold remains unavailable"},
        ],
    }), encoding="utf-8")

    report = render_report(tmp_path)
    text = report.read_text(encoding="utf-8")
    assert "Support and scale" in text
    assert "Raw, SVC, and common-valid support remain visible." in text
    assert "missing_curve.csv" in text
    assert "Saved file index" in text and "tables/support_curves.csv" in text
    assert "no_stable_threshold" in text
    assert "provider unavailable" in text
    assert "reconstructed-state labels" in text
    assert "Gain candidate" in text


def test_support_and_threshold_figures_are_saved_from_tables(tmp_path: Path):
    pytest = __import__("pytest")
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from revise_analysis.plotting import (
        plot_support_curves,
        plot_threshold_bootstrap,
        plot_window_field,
    )

    support = pd.DataFrame({
        "curve": ["svc", "svc", "raw", "raw", "common-valid", "common-valid"],
        "window_side_microns": [16, 32, 16, 32, 16, 32],
        "n_valid_windows": [8, 9, 7, 8, 6, 7],
        "retained_parent_unit_fraction": [0.8, 0.8, 0.7, 0.7, None, None],
    })
    thresholds = pd.DataFrame({"bootstrap_threshold": [1.0, 1.1, 0.9]})
    support_path = plot_support_curves(support, tmp_path / "support.png")
    threshold_path = plot_threshold_bootstrap(thresholds, tmp_path / "threshold.png")
    assert support_path.exists() and threshold_path.exists()

    failed_threshold = pd.DataFrame({
        "window_x": [0.0, 1.0],
        "window_y": [0.0, 1.0],
        "neff": [1.5, 2.0],
        "valid_window": [True, True],
        "in_region": pd.Series([pd.NA, pd.NA], dtype="boolean"),
    })
    field_path = plot_window_field(
        failed_threshold,
        tmp_path / "field.png",
        value_column="neff",
        title="state_Fibroblast_region_windows",
    )
    assert field_path.exists()


def test_impact_reading_summary_is_bounded_deterministic_and_evidence_linked(tmp_path: Path):
    (tmp_path / "tables").mkdir()
    (tmp_path / "tables/input_overview.csv").write_text(
        "side,n_units,n_genes\nraw,340837,18000\nsvc,17455,17000\n", encoding="utf-8"
    )
    (tmp_path / "tables/reconstruction_label_summary.csv").write_text(
        "reconstruction_label,n_units\na,10000\nb,7455\n", encoding="utf-8"
    )
    (tmp_path / "tables/membership_summary_All.json").write_text(json.dumps({
        "summary": [{"n_units": 100, "st_unit_change_fraction": 0.25}], "mapping": [],
    }), encoding="utf-8")
    (tmp_path / "tables/moran_shared_genes_All.csv").write_text(
        "gene_id,delta_moran,comparison_available\na,0.2,true\nb,-0.1,true\n", encoding="utf-8"
    )
    (tmp_path / "tables/integrated_region_anatomy_summary_All.csv").write_text(
        "region_kind,in_region,n_units,denominator_units,anatomy_region\nstate,true,20,80,Tumor\nstate,false,60,80,Normal\n",
        encoding="utf-8",
    )
    outputs = {path.name: path.relative_to(tmp_path).as_posix() for path in (tmp_path / "tables").iterdir()}
    result = {"analysis": "reconstruction_impact", "parameters": {"scopes": ["All"]}, "outputs": outputs}
    facts = reading_summary(tmp_path, result)
    assert len(facts) == 5
    assert all(set(fact) == {"text", "anchor"} for fact in facts)
    assert "RAW 340837" in facts[0]["text"] and "SVC 17455" in facts[0]["text"]
    assert "对象范围独立" in facts[0]["text"] and facts[0]["anchor"].startswith("artifact-")
    assert "20/80" in facts[-1]["text"]


def test_impact_report_uses_fixed_question_tree_scope_order_and_single_image_embedding(tmp_path: Path):
    tables, figures = tmp_path / "tables", tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()
    (tables / "input_overview.csv").write_text(
        "side,n_units,n_genes\nraw,340837,18000\nsvc,17455,17000\n", encoding="utf-8"
    )
    (tables / "reconstruction_label_summary.csv").write_text(
        "reconstruction_label,n_units\na,17455\n", encoding="utf-8"
    )
    (tables / "integrated_label_composition_All.csv").write_text(
        "state_region,anatomy_region,reconstruction_label,n_units,denominator_units,fraction\n"
        "true,Tumor,a,20,80,0.25\n", encoding="utf-8"
    )
    (tables / "integrated_label_units_All.csv").write_text(
        "window_id,state_region,reconstruction_label,scope\nw1,true,a,All\n", encoding="utf-8"
    )
    (tables / "misleading_moran_notes.csv").write_text("value\n1\n", encoding="utf-8")
    image = figures / "relationships_integrated_label_composition_All.png"
    image.write_bytes(b"not-decoded-by-static-renderer")
    paths = [*tables.iterdir(), image]
    outputs = {path.name: path.relative_to(tmp_path).as_posix() for path in paths}
    result = {
        "analysis": "reconstruction_impact", "sample_id": "P2", "status": "partial",
        "parameters": {"scopes": ["All", "Fibroblast"], "parent_window_side_microns": {"All": 40, "Fibroblast": 40}},
        "outputs": outputs,
        "unavailable": [{"component": "membership:Fibroblast", "reason": "no common unit IDs"}],
        "sections": [{"id": "support", "title": "Support and scale", "description": "legacy description"}],
    }
    (tmp_path / "result.json").write_text(json.dumps(result), encoding="utf-8")

    text = render_report(tmp_path).read_text(encoding="utf-8")
    assert all(title in text for title in (
        "重建后特征与两侧差异", "状态与差异的空间位置", "空间位置上的分子与成员关联",
    ))
    assert text.index("scope-label-composition-All") < text.index("scope-label-composition-Fibroblast")
    assert 'id="scope-label-composition-All"' in text
    assert "denominator_units=80" in text
    assert text.count('src="figures/relationships_integrated_label_composition_All.png"') == 1
    assert "misleading_moran_notes.csv" in text  # file index only, not semantically guessed into Moran
    genes_start = text.index('id="topic-genes"')
    genes_end = text.index('id="topic-program-overall"')
    assert "misleading_moran_notes.csv" not in text[genes_start:genes_end]
    assert "不可用：membership:Fibroblast" in text and "此状态不表示没有变化" in text
    assert "Support and scale" in text and "legacy description" in text
    assert "hashchange" in text and "popstate" in text and "IntersectionObserver" in text


def test_impact_figures_render_one_section_from_registered_tables(tmp_path: Path):
    pytest = __import__("pytest")
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from revise_analysis.plotting import render_impact_figures

    tables = tmp_path / "tables"
    tables.mkdir()
    (tables / "input_manifest.json").write_text(json.dumps({
        "spatial": {"microns_per_coordinate": 2},
    }), encoding="utf-8")
    pd.DataFrame({
        "window_id": ["0_0", "1_0"], "window_x": [5.0, 15.0], "window_y": [5.0, 5.0],
        "n_units": [5, 2], "valid_window": [True, False],
        "window_side_length": [10.0, 10.0], "window_side_microns": [20.0, 20.0],
    }).to_csv(tables / "window_support_grid_All.csv", index=False)
    anatomy = pd.DataFrame({
        "window_id": ["0_0", "1_0"], "window_x": [5.0, 15.0], "window_y": [5.0, 5.0],
        "level1_region": ["Tumor", "Normal"],
    })
    anatomy.to_csv(tables / "raw_anatomy_windows.csv", index=False)
    state = pd.DataFrame({
        "window_id": ["0_0", "1_0"], "window_x": [5.0, 15.0], "window_y": [5.0, 5.0],
        "n_units": [5, 2], "valid_window": [True, False],
        "k_obs": [2.0, None], "entropy": [.6, None], "neff": [1.8, None], "evenness": [.9, None],
        "in_region": [True, None],
    })
    state.to_csv(tables / "state_All_region_windows.csv", index=False)
    svc = state.drop(columns=["in_region"])
    svc.to_csv(tables / "window_diversity_svc_All.csv", index=False)
    raw = svc.copy()
    raw.loc[0, ["k_obs", "entropy", "neff", "evenness"]] = [1.0, .2, 1.2, .8]
    raw.to_csv(tables / "window_diversity_raw_All.csv", index=False)
    outputs = {path.name: path.relative_to(tmp_path).as_posix() for path in tables.iterdir()}
    parameters = {
        "anatomy_window_side_microns": 20,
        "parent_window_side_microns": {"All": 20},
        "min_window_units": 4, "n_window_draws": 200,
    }

    diversity = render_impact_figures(tmp_path, outputs, parameters, section="diversity")
    assert diversity["figures/window_diversity_svc_All.png"] == "diversity"
    assert diversity["figures/baseline_window_metrics_All.png"] == "diversity"
    assert not (tmp_path / "figures/anatomy_neff_state_All.png").exists()

    support = render_impact_figures(tmp_path, outputs, parameters, section="support")
    anatomy_figures = render_impact_figures(tmp_path, outputs, parameters, section="anatomy")
    assert support["figures/window_support_grid_All.png"] == "support"
    assert set(support.values()) == {"support"}
    assert anatomy_figures["figures/anatomy_neff_state_All.png"] == "anatomy"
    assert all((tmp_path / relative).stat().st_size > 0 for relative in {*diversity, *support, *anatomy_figures})


def test_impact_report_embeds_new_exact_paths_once_and_explains_neff_parameters(tmp_path: Path):
    tables, figures = tmp_path / "tables", tmp_path / "figures"
    tables.mkdir()
    figures.mkdir()
    pd.DataFrame({
        "window_id": ["0_0"], "valid_window": [True], "neff": [2.0],
    }).to_csv(tables / "window_diversity_svc_All.csv", index=False)
    pd.DataFrame({
        "scope": ["All"], "region_kind": ["state"], "anatomy_region": ["Tumor"],
        "n_units": [5], "n_inside": [3], "n_outside": [1], "n_unknown": [1], "n_known": [4],
        "fraction_inside_known": [.75], "fraction_unknown": [.2],
        "fraction_inside_parent": [.3], "denominator_parent_units": [10],
    }).to_csv(tables / "integrated_anatomy_conditionals_All.csv", index=False)
    for name in ("baseline_window_metrics_All.png", "anatomy_neff_state_All.png",
                 "relationships_integrated_anatomy_conditionals_All.png"):
        (figures / name).write_bytes(b"saved-figure")
    paths = [*tables.iterdir(), *figures.iterdir()]
    result = {
        "analysis": "reconstruction_impact", "sample_id": "P2", "status": "succeeded",
        "parameters": {"scopes": ["All"], "parent_window_side_microns": {"All": 40},
                       "min_window_units": 4, "n_window_draws": 200},
        "outputs": {path.name: path.relative_to(tmp_path).as_posix() for path in paths},
    }
    (tmp_path / "result.json").write_text(json.dumps(result), encoding="utf-8")

    text = render_report(tmp_path).read_text(encoding="utf-8")
    assert "每次抽取 4 个单位，重复 200 次后取指标均值" in text
    for name in ("baseline_window_metrics_All.png", "anatomy_neff_state_All.png",
                 "relationships_integrated_anatomy_conditionals_All.png"):
        assert text.count(f'src="figures/{name}"') == 1
