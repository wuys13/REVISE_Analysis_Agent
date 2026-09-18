"""Composable reconstruction-impact workflow.

Each stage writes a small saved artifact.  The functions are deliberately usable
from a notebook so scientific exploration does not have to call one monolith.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from ._shared import (DEFAULT_SCOPES, analysis_parameters, coordinates, deterministic_subset, save_json, save_table, unavailable)


def prepare_side(sample: Any, side: str, *, sample_n_units: int | None, random_state: int):
    """Return one independently sampled side and its canonical broad labels."""
    adata = deterministic_subset(getattr(sample, side), sample_n_units, random_state)
    labels = sample.labels(side).reindex(adata.obs_names)
    return adata, labels


def compute_impact(sample: Any, output_dir: Path, parameters: dict | None = None) -> dict:
    """Alias for :func:`run`, retained as the notebook-facing workflow name."""
    return run(sample, output_dir, parameters)


def compute_partitions(adata: Any, *, resolution: float, n_top_genes: int, random_state: int):
    from revise_analysis.methods.partition import compute_partition
    return compute_partition(
        adata, resolution=resolution, n_top_genes=n_top_genes,
        random_state=random_state, min_genes=0, min_cells=0,
    )


def partition_scope(adata: Any, labels: pd.Series, scope: str, *, resolution: float, n_top_genes: int, random_state: int):
    """Partition one side/scope without requiring units from the other side."""
    if scope != "All":
        selected = labels.index[labels == scope]
        if len(selected) == 0:
            raise ValueError(f"no units labelled {scope!r}")
        adata = adata[selected]
    return compute_partitions(adata, resolution=resolution, n_top_genes=n_top_genes, random_state=random_state)


def _partition_support_reason(adata: Any) -> str | None:
    if adata.n_obs < 3:
        return "Leiden partitioning requires at least three units"
    if adata.n_vars < 2:
        return "Leiden partitioning requires at least two genes"
    return None


def raw_anatomy_context(sample: Any, *, tumor_label: str = "Tumor",
                        normal_label: str = "Intestinal Epithelial") -> pd.DataFrame:
    """Full Raw broad anatomy context, never an inferred Level2 mapping."""
    raw = sample.raw
    labels = sample.labels("raw", sample.broad_key)
    frame = coordinates(raw, sample.spatial_key)
    frame["broad_label"] = labels.reindex(frame.index)
    frame["anatomy_class"] = frame["broad_label"].map(
        lambda value: _anatomy_class(value, tumor_label=tumor_label, normal_label=normal_label)
    )
    return frame


def native_window_diversity(sample: Any, side: str, labels: pd.Series, *, adata: Any | None = None,
                            window_side_length: float, origin: tuple[float, float], min_units: int,
                            n_draws: int, random_state: int) -> pd.DataFrame:
    """Calculate diversity in one side's native units and coordinates."""
    from revise_analysis.methods.regions import compute_window_diversity
    adata = getattr(sample, side) if adata is None else adata
    return compute_window_diversity(coordinates(adata, sample.spatial_key).reindex(labels.index), labels,
                                    window_side_length=window_side_length, origin=origin,
                                    min_units=min_units, n_draws=n_draws, random_state=random_state)


def run(sample: Any, output_dir: Path, parameters: dict | None = None) -> dict:
    """Run independent impact components and save their inspectable artifacts."""
    supplied = analysis_parameters(parameters, {
        "scopes", "resolution", "n_top_genes", "sample_n_units", "random_state", "membership_same_units",
        "window_side_microns", "min_window_units", "n_window_draws", "region_n_bootstrap",
        "matched_k_resolutions",
        "gene_sets", "geneset_path", "gene_set_names", "pathway_auc_threshold", "moran_n_neighbors",
        "anatomy_tumor_label", "anatomy_normal_label",
    })
    scopes = tuple(supplied.get("scopes", DEFAULT_SCOPES))
    if (not scopes or any(not isinstance(scope, str) or not scope.strip() for scope in scopes)
            or len(set(scopes)) != len(scopes)):
        raise ValueError("scopes must be unique non-empty broad-label names")
    if not isinstance(supplied.get("membership_same_units", False), bool):
        raise ValueError("membership_same_units must be a boolean")
    if len({_slug(scope) for scope in scopes}) != len(scopes):
        raise ValueError("scopes collide after conversion to artifact filenames")
    effective = {
        "scopes": list(scopes),
        "resolution": float(supplied.get("resolution", 0.5)),
        "n_top_genes": int(supplied.get("n_top_genes", 2000)),
        "sample_n_units": supplied.get("sample_n_units"),
        "random_state": int(supplied.get("random_state", 42)),
        "membership_same_units": bool(supplied.get("membership_same_units", False)),
        "window_side_microns": float(supplied.get("window_side_microns", 100.0)),
        "min_window_units": int(supplied.get("min_window_units", 4)),
        "n_window_draws": int(supplied.get("n_window_draws", 200)),
        "region_n_bootstrap": int(supplied.get("region_n_bootstrap", 500)),
        "matched_k_resolutions": list(supplied.get("matched_k_resolutions", (0.2, 0.4, 0.6, 0.8, 1.0))),
        "gene_sets": supplied.get("gene_sets"), "geneset_path": supplied.get("geneset_path"), "gene_set_names": supplied.get("gene_set_names"),
        "pathway_auc_threshold": float(supplied.get("pathway_auc_threshold", 0.05)),
        "moran_n_neighbors": int(supplied.get("moran_n_neighbors", 6)),
        "anatomy_tumor_label": str(supplied.get("anatomy_tumor_label", "Tumor")),
        "anatomy_normal_label": str(supplied.get("anatomy_normal_label", "Intestinal Epithelial")),
    }
    if (effective["resolution"] <= 0 or effective["n_top_genes"] < 1
            or effective["window_side_microns"] <= 0 or effective["min_window_units"] < 1
            or effective["n_window_draws"] < 1 or effective["region_n_bootstrap"] < 1):
        raise ValueError("partition, window, draw, and threshold parameters must be positive")
    if (not effective["matched_k_resolutions"]
            or any(not isinstance(value, (int, float)) or float(value) <= 0
                   for value in effective["matched_k_resolutions"])):
        raise ValueError("matched_k_resolutions must be a nonempty sequence of positive resolutions")
    effective["matched_k_resolutions"] = sorted({float(value) for value in effective["matched_k_resolutions"]})
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs, missing = {}, []

    side_data: dict[str, tuple[Any, pd.Series]] = {}
    overview = []
    for side in ("raw", "svc"):
        adata, labels = prepare_side(sample, side, sample_n_units=effective["sample_n_units"], random_state=effective["random_state"] + (side == "svc"))
        side_data[side] = (adata, labels)
        overview.append({"side": side, "n_units": int(adata.n_obs), "n_genes": int(adata.n_vars)})
    save_table(pd.DataFrame(overview), output_dir, "tables/overview.csv", outputs)

    partition_labels: dict[tuple[str, str], pd.Series] = {}
    partition_rows = []
    for side, (adata, labels) in side_data.items():
        for scope in scopes:
            component = f"partition:{side}:{scope}"
            if scope != "All" and not (labels == scope).any():
                unavailable(component, f"no {scope!r} labels on {side}", missing)
                continue
            candidate = adata if scope == "All" else adata[labels.index[labels == scope]]
            if reason := _partition_support_reason(candidate):
                unavailable(component, reason, missing)
                continue
            try:
                partition = partition_scope(adata, labels, scope, resolution=effective["resolution"], n_top_genes=effective["n_top_genes"], random_state=effective["random_state"])
                memberships = partition["labels"].rename("partition")
                partition_labels[side, scope] = memberships
                table = memberships.to_frame()
                table["scope"] = scope
                table["side"] = side
                save_table(table, output_dir, f"tables/partitions_{side}_{_slug(scope)}.csv", outputs)
                partition_rows.append({"side": side, "scope": scope, **partition["summary"]})
            except (ImportError, KeyError) as exc:
                unavailable(component, str(exc), missing)
    if partition_rows:
        save_table(pd.DataFrame(partition_rows), output_dir, "tables/partition_summary.csv", outputs)

    has_raw_level2 = _write_raw_level2_baseline(sample, output_dir, outputs, missing)
    _membership_and_change_map(sample, output_dir, effective, side_data, partition_labels, outputs, missing)
    try:
        geometry = _spatial_geometry(sample, effective)
    except (KeyError, ValueError) as exc:
        geometry = None
        unavailable("spatial_geometry", str(exc), missing)
    if geometry is not None:
        _write_raw_anatomy(sample, output_dir, geometry, effective, outputs, missing)
        _write_state_and_gain(sample, output_dir, effective, geometry, side_data, partition_labels,
                              has_raw_level2, outputs, missing)
    _write_moran_and_pathways(sample, output_dir, effective, side_data, outputs, missing)
    _write_figures(output_dir, outputs, missing)
    status = "partial" if outputs and missing else ("succeeded" if outputs else "skipped")
    return {"status": status, "parameters": effective, "outputs": outputs, "unavailable": missing}


def _membership_and_change_map(sample, output_dir, effective, side_data, partition_labels, outputs, missing) -> None:
    if not effective["membership_same_units"]:
        return
    from revise_analysis.methods.partition import compare_membership
    for scope in effective["scopes"]:
        raw = partition_labels.get(("raw", scope))
        # Membership is the sole paired calculation.  Its Raw-defined cohort
        # comes from the Raw partition sample, but it must intersect the full
        # SVC object rather than an independently sampled SVC analysis side.
        svc_adata = sample.svc
        component = f"membership:{scope}"
        if raw is None:
            unavailable(component, "Raw parent partition unavailable", missing)
            continue
        # The Raw partition defines the parent population.  SVC labels are
        # recomputed on those common IDs; SVC's own broad annotation never
        # silently redefines the membership population.
        common = raw.index[raw.index.isin(svc_adata.obs_names)]
        if common.empty:
            unavailable(component, "no common unit IDs", missing)
            continue
        if reason := _partition_support_reason(svc_adata[common]):
            unavailable(component, f"SVC partition for Raw-defined parent unavailable: {reason}", missing)
            continue
        try:
            svc_membership = compute_partitions(
                svc_adata[common], resolution=effective["resolution"], n_top_genes=effective["n_top_genes"],
                random_state=effective["random_state"],
            )["labels"]
        except (ImportError, KeyError) as exc:
            unavailable(component, f"SVC partition for Raw-defined parent unavailable: {exc}", missing)
            continue
        comparison = compare_membership(raw.loc[common], svc_membership)
        save_table(comparison.assignments, output_dir, f"tables/membership_{_slug(scope)}.csv", outputs)
        save_table(comparison.contingency, output_dir, f"tables/membership_contingency_{_slug(scope)}.csv", outputs)
        save_json({"summary": comparison.summary, "mapping": comparison.mapping}, output_dir, f"tables/membership_summary_{_slug(scope)}.json", outputs)
        # Raw coordinates are a location reference only for the explicitly common IDs.
        try:
            changed = coordinates(sample.raw, sample.spatial_key).reindex(comparison.assignments.index).join(comparison.assignments)
            save_table(changed, output_dir, f"tables/spatial_changed_map_{_slug(scope)}.csv", outputs)
        except (KeyError, ValueError) as exc:
            unavailable(f"spatial_changed_map:{scope}", str(exc), missing)
        _write_matched_k_membership(
            raw.loc[common], svc_adata[common], scope, effective, output_dir, outputs, missing,
        )


def _write_matched_k_membership(raw: pd.Series, svc_cohort: Any, scope: str, effective: dict,
                                output_dir: Path, outputs: dict, missing: list[dict]) -> None:
    """Run the optional controlled-K membership diagnostic on the Raw-defined cohort.

    This is deliberately saved beside, never merged into, the same-resolution
    membership comparison.  It changes a clustering control, not unit scope.
    """
    from revise_analysis.methods.partition import compare_membership

    target = int(raw.nunique())
    if reason := _partition_support_reason(svc_cohort):
        save_table(pd.DataFrame([{"resolution": pd.NA, "target_raw_clusters": target,
                                  "svc_clusters": pd.NA, "status": "unavailable", "reason": reason}]),
                   output_dir, f"tables/membership_matched_k_sweep_{_slug(scope)}.csv", outputs)
        unavailable(f"membership_matched_k:{scope}", reason, missing)
        return
    sweep = []
    selected = None
    for resolution in effective["matched_k_resolutions"]:
        try:
            partition = compute_partitions(
                svc_cohort, resolution=resolution, n_top_genes=effective["n_top_genes"],
                random_state=effective["random_state"],
            )
        except (ImportError, KeyError) as exc:
            sweep.append({"resolution": resolution, "target_raw_clusters": target,
                          "svc_clusters": None, "status": "unavailable", "reason": str(exc)})
            continue
        labels = partition["labels"]
        observed = int(labels.nunique())
        row = {"resolution": resolution, "target_raw_clusters": target, "svc_clusters": observed,
               "status": "exact" if observed == target else "not_exact", "reason": ""}
        sweep.append(row)
        if observed == target and selected is None:
            selected = (resolution, labels)
    save_table(pd.DataFrame(sweep), output_dir, f"tables/membership_matched_k_sweep_{_slug(scope)}.csv", outputs)
    if selected is None:
        unavailable(f"membership_matched_k:{scope}", f"no candidate produced Raw target K={target}", missing)
        return
    resolution, labels = selected
    comparison = compare_membership(raw, labels)
    save_table(comparison.assignments, output_dir, f"tables/membership_matched_k_{_slug(scope)}.csv", outputs)
    save_table(comparison.contingency, output_dir, f"tables/membership_matched_k_contingency_{_slug(scope)}.csv", outputs)
    save_json({"target_raw_clusters": target, "selected_svc_resolution": resolution,
               "summary": comparison.summary, "mapping": comparison.mapping}, output_dir,
              f"tables/membership_matched_k_summary_{_slug(scope)}.json", outputs)


def _write_raw_level2_baseline(sample, output_dir, outputs, missing) -> bool:
    """Persist an upstream Raw subtype baseline without deriving replacement labels."""
    if sample.subtype_key not in sample.raw.obs:
        unavailable("raw_level2_baseline", f"Raw column {sample.subtype_key!r} is unavailable", missing)
        return False
    labels = sample.raw.obs[sample.subtype_key].astype("string")
    if labels.isna().any():
        unavailable("raw_level2_baseline", f"Raw column {sample.subtype_key!r} contains missing values", missing)
        return False
    baseline = pd.DataFrame({"unit_id": sample.raw.obs_names.astype(str), "raw_level2": labels.to_numpy()})
    counts = baseline.groupby("raw_level2", dropna=False).size().rename("n_units").reset_index()
    save_table(baseline, output_dir, "tables/raw_level2_baseline.csv", outputs)
    save_table(counts, output_dir, "tables/raw_level2_baseline_summary.csv", outputs)
    return True


def _spatial_geometry(sample, effective) -> dict:
    scale = sample.config.get("spatial", {}).get("microns_per_coordinate")
    if scale is None or float(scale) <= 0:
        raise ValueError("spatial.microns_per_coordinate is required for micrometer windows")
    full_raw = coordinates(sample.raw, sample.spatial_key)
    return {"origin": (float(full_raw.x.min()), float(full_raw.y.min())),
            "window_side_length": effective["window_side_microns"] / float(scale)}


def _write_raw_anatomy(sample, output_dir, geometry, effective, outputs, missing) -> None:
    try:
        from revise_analysis.methods.regions import assign_anatomy_candidates, assign_square_windows
        context = raw_anatomy_context(
            sample, tumor_label=effective["anatomy_tumor_label"],
            normal_label=effective["anatomy_normal_label"],
        )
        assignments = assign_square_windows(context[["x", "y"]], window_side_length=geometry["window_side_length"], origin=geometry["origin"])
        candidates = assign_anatomy_candidates(
            assignments, context["broad_label"], tumor_label=effective["anatomy_tumor_label"],
            normal_source_label=effective["anatomy_normal_label"],
        )
        context = context.join(assignments[["window_id"]]).join(candidates.set_index("window_id")[["level1_region"]], on="window_id")
        save_table(context, output_dir, "tables/raw_anatomy_context.csv", outputs)
        save_table(candidates, output_dir, "tables/raw_anatomy_windows.csv", outputs)
    except (ImportError, KeyError, ValueError) as exc:
        unavailable("raw_anatomy_context", str(exc), missing)


def _write_state_and_gain(sample, output_dir, effective, geometry, side_data, partition_labels,
                          has_raw_level2, outputs, missing) -> None:
    windows = {}
    for (side, scope), labels in partition_labels.items():
        try:
            # ``side_data`` is independently sampled upstream.  Neither the
            # coordinates nor the rarefaction generator is reused across sides.
            adata, _ = side_data[side]
            seed = effective["random_state"] + (0 if side == "raw" else 1)
            table = native_window_diversity(
                sample, side, labels, adata=adata,
                window_side_length=geometry["window_side_length"], origin=geometry["origin"],
                min_units=effective["min_window_units"], n_draws=effective["n_window_draws"],
                random_state=seed,
            )
            windows[side, scope] = table
            save_table(table, output_dir, f"tables/window_diversity_{side}_{_slug(scope)}.csv", outputs)
        except (ImportError, KeyError, ValueError) as exc:
            unavailable(f"window_diversity:{side}:{scope}", str(exc), missing)
    # State is a native SVC question for each parent population.  It is not a
    # Raw/SVC difference and does not need common observation IDs.
    for scope in effective["scopes"]:
        if scope == "All":
            continue
        state = windows.get(("svc", scope))
        if state is None:
            unavailable(f"state:{scope}", "SVC parent partition local diversity unavailable", missing)
            continue
        _write_region(state, output_dir, f"state_{_slug(scope)}", "neff", geometry, effective, outputs, missing)
    if has_raw_level2:
        _write_raw_level2_window_baseline(
            sample, output_dir, effective, geometry, side_data, windows, outputs, missing,
        )
    for scope in effective["scopes"]:
        raw, svc = windows.get(("raw", scope)), windows.get(("svc", scope))
        if raw is None or svc is None:
            continue
        raw, svc = raw.set_index("window_id"), svc.set_index("window_id")
        common = raw.index.intersection(svc.index)
        valid = raw.loc[common, "valid_window"].astype(bool) & svc.loc[common, "valid_window"].astype(bool)
        if not valid.any():
            unavailable(f"gain:{scope}", "no common valid spatial windows", missing)
            continue
        gain = pd.DataFrame({
            "raw_neff": raw.loc[common, "neff"], "svc_neff": svc.loc[common, "neff"],
            "raw_n_units": raw.loc[common, "n_units"], "svc_n_units": svc.loc[common, "n_units"],
            "window_x": raw.loc[common, "window_x"], "window_y": raw.loc[common, "window_y"],
        }).loc[valid]
        gain["gain_neff"] = gain["svc_neff"] - gain["raw_neff"]
        # This explicit support count is only for the common-window extent;
        # diversity on each side remains native and independently sampled.
        gain["n_units"] = gain[["raw_n_units", "svc_n_units"]].min(axis=1).astype(int)
        gain = gain.reset_index()
        gain["valid_window"] = True
        name = f"gain_{_slug(scope)}"
        save_table(gain, output_dir, f"tables/{name}_common_valid_windows.csv", outputs)
        _write_region(gain, output_dir, name, "gain_neff", geometry, effective, outputs, missing, positive_only=True)


def _write_raw_level2_window_baseline(sample, output_dir, effective, geometry, side_data, svc_windows,
                                      outputs, missing) -> None:
    """Save the supplied Raw Level2 local baseline beside native SVC parent state.

    Raw Level2 and a de-novo SVC Leiden partition are different label systems.
    Their common-window delta is therefore a descriptive reference table, not a
    Gain field and never feeds State/Gain region selection.
    """
    raw_adata, raw_broad = side_data["raw"]
    raw_level2 = sample.labels("raw", sample.subtype_key).reindex(raw_adata.obs_names)
    for scope in effective["scopes"]:
        if scope == "All":
            selected = raw_adata.obs_names
        else:
            selected = raw_broad.index[raw_broad == scope]
        if len(selected) == 0:
            unavailable(f"raw_level2_baseline:{scope}", f"no Raw {scope!r} units", missing)
            continue
        component = f"raw_level2_baseline:{scope}"
        try:
            baseline = native_window_diversity(
                sample, "raw", raw_level2.loc[selected], adata=raw_adata[selected],
                window_side_length=geometry["window_side_length"], origin=geometry["origin"],
                min_units=effective["min_window_units"], n_draws=effective["n_window_draws"],
                random_state=effective["random_state"] + 17,
            )
        except (ImportError, KeyError, ValueError) as exc:
            unavailable(component, str(exc), missing)
            continue
        save_table(baseline, output_dir, f"tables/window_diversity_raw_level2_baseline_{_slug(scope)}.csv", outputs)
        svc = svc_windows.get(("svc", scope))
        if svc is None:
            unavailable(f"raw_level2_vs_svc:{scope}", "SVC parent local diversity unavailable", missing)
            continue
        raw = baseline.set_index("window_id")
        svc = svc.set_index("window_id")
        common = raw.index.intersection(svc.index)
        valid = raw.loc[common, "valid_window"].astype(bool) & svc.loc[common, "valid_window"].astype(bool)
        if not valid.any():
            unavailable(f"raw_level2_vs_svc:{scope}", "no common valid spatial windows", missing)
            continue
        comparison = pd.DataFrame({
            "raw_level2_baseline_neff": raw.loc[common, "neff"],
            "svc_parent_leiden_neff": svc.loc[common, "neff"],
            "raw_n_units": raw.loc[common, "n_units"], "svc_n_units": svc.loc[common, "n_units"],
            "window_x": raw.loc[common, "window_x"], "window_y": raw.loc[common, "window_y"],
        }).loc[valid]
        comparison["descriptive_delta_neff"] = (
            comparison["svc_parent_leiden_neff"] - comparison["raw_level2_baseline_neff"]
        )
        comparison = comparison.reset_index()
        comparison["valid_window"] = True
        save_table(comparison, output_dir,
                   f"tables/raw_level2_baseline_vs_svc_{_slug(scope)}_common_valid_windows.csv", outputs)


def _write_region(table, output_dir, name, value_column, geometry, effective, outputs, missing, positive_only=False) -> None:
    from revise_analysis.methods.regions import flag_region_windows, select_region_threshold
    values = pd.to_numeric(table[value_column], errors="coerce")
    threshold_info, audit = select_region_threshold(
        values[values > 0] if positive_only else values,
        n_bootstrap=effective["region_n_bootstrap"], random_state=effective["random_state"],
    )
    if threshold_info.get("threshold") is None:
        save_json(threshold_info, output_dir, f"tables/{name}_threshold.json", outputs)
        if not audit.empty:
            save_table(audit, output_dir, f"tables/{name}_threshold_bootstrap.csv", outputs)
        unavailable(name, f"stable region threshold unavailable: {threshold_info['status']}", missing)
        return
    summary, flagged = flag_region_windows(
        table, threshold=float(threshold_info["threshold"]), value_column=value_column,
        window_side_length=effective["window_side_microns"],
    )
    summary["window_side_microns"] = effective["window_side_microns"]
    summary["region_area_square_microns"] = summary["region_area"]
    save_table(flagged, output_dir, f"tables/{name}_region_windows.csv", outputs)
    save_table(summary, output_dir, f"tables/{name}_region_extent.csv", outputs)
    if not audit.empty:
        save_table(audit, output_dir, f"tables/{name}_threshold_bootstrap.csv", outputs)
    save_json({**threshold_info, "window_side_microns": effective["window_side_microns"]}, output_dir, f"tables/{name}_threshold.json", outputs)


def _write_moran_and_pathways(sample, output_dir, effective, side_data, outputs, missing) -> None:
    from .spatial_autocorrelation import compute_moran, record_uncomputed
    for side in ("raw", "svc"):
        try:
            adata, _ = side_data[side]
            table = compute_moran(adata, spatial_key=sample.spatial_key, n_neighbors=effective["moran_n_neighbors"])
            save_table(table, output_dir, f"tables/moran_{side}.csv", outputs)
            record_uncomputed(table, side, missing)
        except (ImportError, KeyError, ValueError) as exc:
            unavailable(f"moran:{side}", str(exc), missing)
    from .pathway_activity import compute_pathway_scores, resolve_gene_sets
    gene_sets, resource_error = resolve_gene_sets(effective)
    if resource_error:
        unavailable("pathway_resource", resource_error, missing)
        return
    for side in ("raw", "svc"):
        for name, genes in gene_sets.items():
            try:
                adata, _ = side_data[side]
                result = compute_pathway_scores(adata, genes=list(genes), score_name=name,
                                                 auc_threshold=effective["pathway_auc_threshold"], seed=effective["random_state"])
                save_table(result["scores"], output_dir, f"tables/pathway_{side}_{_slug(name)}.csv", outputs)
                save_json({"coverage": result["coverage"], "metadata": result["metadata"]}, output_dir,
                          f"tables/pathway_{side}_{_slug(name)}_metadata.json", outputs)
            except ImportError as exc:
                unavailable(f"pathway:{side}:{name}", str(exc), missing)
            except ValueError as exc:
                unavailable(f"pathway:{side}:{name}", str(exc), missing)


def _write_figures(output_dir, outputs, missing) -> None:
    """Figures are derived solely from just-saved tables, like the report."""
    try:
        from revise_analysis.plotting import (
            plot_anatomy_context, plot_gain, plot_membership_change,
            plot_moran_distribution, plot_partition_sizes, plot_pathway_scores,
            plot_window_field,
        )
        summary_path = output_dir / "tables/partition_summary.csv"
        if summary_path.exists():
            destination = output_dir / "figures/partition_sizes.png"
            plot_partition_sizes(pd.read_csv(summary_path), destination)
            outputs["figures/partition_sizes.png"] = "figures/partition_sizes.png"
        for side in ("raw", "svc"):
            path = output_dir / f"tables/moran_{side}.csv"
            if path.exists():
                destination = output_dir / f"figures/moran_{side}.png"
                plot_moran_distribution(pd.read_csv(path), destination, title=f"{side.upper()} Moran I")
                outputs[str(destination.relative_to(output_dir))] = str(destination.relative_to(output_dir))
        for path in (output_dir / "tables").glob("gain_*_common_valid_windows.csv"):
            destination = output_dir / "figures" / f"{path.stem}.png"
            plot_gain(pd.read_csv(path), destination)
            outputs[str(destination.relative_to(output_dir))] = str(destination.relative_to(output_dir))
        anatomy = output_dir / "tables/raw_anatomy_context.csv"
        if anatomy.exists():
            destination = output_dir / "figures/raw_anatomy_context.png"
            plot_anatomy_context(pd.read_csv(anatomy), destination)
            outputs[str(destination.relative_to(output_dir))] = str(destination.relative_to(output_dir))
        for path in (output_dir / "tables").glob("spatial_changed_map_*.csv"):
            destination = output_dir / "figures" / f"{path.stem}.png"
            plot_membership_change(pd.read_csv(path), destination, title=path.stem.replace("_", " "))
            outputs[str(destination.relative_to(output_dir))] = str(destination.relative_to(output_dir))
        for path in (output_dir / "tables").glob("*_region_windows.csv"):
            destination = output_dir / "figures" / f"{path.stem}.png"
            value_column = "gain_neff" if "gain_" in path.name else "neff"
            plot_window_field(pd.read_csv(path), destination, value_column=value_column,
                              title=path.stem.replace("_", " "))
            outputs[str(destination.relative_to(output_dir))] = str(destination.relative_to(output_dir))
        for path in (output_dir / "tables").glob("pathway_*.csv"):
            destination = output_dir / "figures" / f"{path.stem}.png"
            score = pd.read_csv(path, index_col=0).iloc[:, 0]
            plot_pathway_scores(score, destination, title=path.stem.replace("_", " "))
            outputs[str(destination.relative_to(output_dir))] = str(destination.relative_to(output_dir))
    except ImportError as exc:
        unavailable("figures", str(exc), missing)


def _anatomy_class(label: object, *, tumor_label: str = "Tumor",
                   normal_label: str = "Intestinal Epithelial") -> str:
    value = str(label)
    if value == tumor_label:
        return "Tumor"
    if value == normal_label:
        return "Normal"
    if value == "Interface":
        return "Interface"
    return "Other"


def _slug(value: str) -> str:
    return "".join(char if char.isalnum() else "_" for char in value).strip("_") or "scope"
