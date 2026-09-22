from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from revise_analysis.io import Sample
from revise_analysis.analyses import reconstruction_impact as impact
from revise_analysis.analyses._shared import unit_id_digest
from revise_analysis.analyses.impact_tables import anatomy_conditionals, common_valid_delta
from revise_analysis.methods import regions


def sample():
    raw = ad.AnnData(np.ones((12, 3)))
    raw.obs_names = [f'u{i}' for i in range(12)]
    raw.obs['Level1'] = ['Fibroblast'] * 6 + ['Tumor'] * 6
    raw.obsm['spatial'] = np.c_[np.arange(12) * 5, np.zeros(12)]
    svc = raw.copy()
    svc.obs['SVC_cluster'] = ['a', 'b'] * 6
    return Sample('test', raw, svc, {'spatial': {'microns_per_coordinate': 1},
                  'expression': {side: {'identity': 'unknown'} for side in ('raw', 'svc')}}, Path('sample.yaml'))


def workflow(tmp_path):
    return impact.ImpactWorkflow(sample(), tmp_path, {'scopes': ['All'], 'min_window_units': 2,
                                                     'n_window_draws': 3, 'region_n_bootstrap': 2})


def through(w, stop):
    for stage in impact.STAGE_ORDER:
        w.run_stage(stage)
        if stage == stop:
            break


def test_rerun_stable_unstable_stable_replaces_manifest_and_status(tmp_path, monkeypatch):
    threshold = [1.0]
    monkeypatch.setattr(regions, 'select_region_threshold', lambda *a, **k: (
        {'status': 'ok' if threshold[0] else 'no_stable_threshold', 'threshold': threshold[0]},
        pd.DataFrame({'threshold': [1.0]}) if threshold[0] else pd.DataFrame()))
    w = workflow(tmp_path)
    through(w, 'regions')
    extent = 'tables/state_All_region_extent.csv'
    bootstrap = 'tables/state_All_threshold_bootstrap.csv'
    assert extent in w.outputs
    # A rendered figure from the previous accepted region must also disappear.
    w.outputs['figures/old.png'] = 'figures/old.png'
    w._figure_sections['figures/old.png'] = 'regions'
    threshold[0] = None
    w.run_stage('regions')
    assert extent not in w.outputs and bootstrap not in w.outputs
    assert 'figures/old.png' not in w.outputs
    assert (tmp_path / extent).exists()  # disk retention is not current registration
    assert pd.read_csv(tmp_path / 'tables/state_All_region_windows.csv').in_region.isna().all()
    assert len([r for r in w.stage_records if r['stage'] == 'regions']) == 1
    assert 'tables/state_All_region_windows.csv' in w.stage_records[4]['artifacts']
    threshold[0] = 1.0
    w.run_stage('regions')
    assert extent in w.outputs
    assert not any(r['component'] == 'state_All' for r in w.missing)


def test_parameters_invalidate_region_associations_not_molecular_or_membership(tmp_path, monkeypatch):
    w = workflow(tmp_path)
    calls = {'molecular': 0, 'membership': 0}
    def molecular():
        calls['molecular'] += 1
        w.state['program_scores'] = [pd.DataFrame({'score': [.3], 'side': ['svc'], 'program': ['P']}, index=['u0'])]
    def membership():
        calls['membership'] += 1
        w.state['membership'] = {}
    monkeypatch.setattr(w, 'stage_molecular', molecular)
    monkeypatch.setattr(w, 'stage_membership', membership)
    through(w, 'integration')
    scores = w.state['program_scores']
    p = dict(w.parameters, region_n_bootstrap=3)
    invalidated = w.apply_parameters(p)
    assert set(invalidated) == {'regions', 'integration', 'figures'}
    assert w.state['program_scores'] is scores
    with pytest.raises(RuntimeError, match='requires current stages: regions'):
        w.run_stage('integration')
    w.run_stage('regions')
    w.run_stage('integration')
    assert calls == {'molecular': 1, 'membership': 1}
    assert not w.stage_errors
    w.parameters['region_n_bootstrap'] = 7
    with pytest.raises(RuntimeError, match='apply_parameters'):
        w.run_stage('regions')


def test_scope_before_sampling_and_unavailable_recovers(tmp_path, monkeypatch):
    s = sample()
    w = impact.ImpactWorkflow(s, tmp_path, {'scopes': ['Fibroblast'], 'sample_n_units': 6})
    w.run_stage('baseline')
    assert any(r['component'] == 'partition:raw' for r in w.missing)
    cohorts = []
    def partition(cohort, **kwargs):
        cohorts.append(list(cohort.obs_names))
        return {'labels': pd.Series('a', index=cohort.obs_names), 'summary': {'n_units': cohort.n_obs}}
    monkeypatch.setattr(impact, 'compute_partitions', partition)
    s.config['expression']['raw']['identity'] = 'measured'
    with pytest.raises(RuntimeError, match='rerun input'):
        w.run_stage('baseline')
    w.run_stage('input')
    w.run_stage('baseline')
    assert cohorts == [[f'u{i}' for i in range(6)]]
    assert not any(r['component'] == 'partition:raw' for r in w.missing)
    extras = s.raw[6:].copy()
    extras.obs_names = [f'extra{i}' for i in range(len(extras))]
    s.raw = ad.concat([s.raw, extras])
    w.run_stage('baseline')
    assert cohorts[-1] == cohorts[0]


def test_cohort_manifest_replaces_invalidated_sampling_records(tmp_path, monkeypatch):
    s = sample()
    s.config['expression'] = {side: {'identity': 'fixture'} for side in ('raw', 'svc')}
    w = impact.ImpactWorkflow(s, tmp_path, {'scopes': ['All'], 'sample_n_units': 5})
    monkeypatch.setattr(impact, 'compute_partitions', lambda data, **kw: {
        'labels': pd.Series('a', index=data.obs_names), 'summary': {'n_units': data.n_obs}})
    w.run_stage('input')
    w.run_stage('baseline')
    manifest = json.loads((tmp_path / 'tables/cohort_manifest.json').read_text())
    records = manifest['records']
    baseline = next(row for row in records if row['stage'] == 'baseline')
    assert baseline['selected_n_units'] == 5
    assert baseline['selection_random_state'] == 42
    assert baseline['method_random_state'] == 42
    assert baseline['selected_unit_ids_sha256'] == unit_id_digest(
        w.state['partition_cohorts']['raw', 'All'].obs_names)
    w.apply_parameters({'sample_n_units': 4})
    manifest = json.loads((tmp_path / 'tables/cohort_manifest.json').read_text())
    assert not any(row['stage'] == 'baseline' for row in manifest['records'])
    w.run_stage('baseline')
    manifest = json.loads((tmp_path / 'tables/cohort_manifest.json').read_text())
    baseline = next(row for row in manifest['records'] if row['stage'] == 'baseline')
    assert baseline['selected_n_units'] == 4


def test_gain_extent_keeps_two_real_unit_denominators():
    raw = pd.DataFrame({'window_id': ['a', 'b'], 'window_x': [0, 40], 'window_y': [0, 0],
                        'valid_window': [True, True], 'n_units': [4, 100], 'neff': [1., 2.]})
    svc = raw.assign(n_units=[100, 4], neff=[3., 2.])
    gain = common_valid_delta(raw, svc)
    assert 'n_units' not in gain
    extent, _ = regions.flag_region_windows(gain, threshold=1, value_column='delta_neff', window_side_length=40)
    assert extent.raw_region_unit_fraction.iloc[0] == pytest.approx(4 / 104)
    assert extent.svc_region_unit_fraction.iloc[0] == pytest.approx(100 / 104)
    assert 'region_unit_fraction' not in extent


def test_anatomy_conditional_denominators_keep_unknown():
    units = pd.DataFrame({'anatomy_region': ['Interface'] * 4 + ['Unknown'],
                          'state_region': [True, False, None, None, None], 'gain_region': [None] * 5})
    out = anatomy_conditionals(units)
    row = out[(out.region_kind == 'state') & (out.anatomy_region == 'Interface')].iloc[0]
    assert (row.n_inside, row.n_outside, row.n_unknown, row.n_known) == (1, 1, 2, 2)
    assert row.fraction_inside_known == .5 and row.fraction_unknown == .5
    assert row.fraction_inside_parent == .2
    assert out.loc[out.region_kind == 'gain', ['fraction_inside_known', 'fraction_inside_parent']].isna().all().all()


def test_region_rerun_reuses_actual_molecular_stage_scores(tmp_path, monkeypatch):
    from revise_analysis.analyses import pathway_activity, spatial_autocorrelation
    s = sample()
    s.config['expression'] = {side: {'identity': 'fixture'} for side in ('raw', 'svc')}
    calls = {'auc': 0, 'moran': 0}
    def auc(data, **kwargs):
        calls['auc'] += 1
        return {'scores': pd.Series(.4, index=data.obs_names), 'metadata': {'rank_cutoff': 1},
                'coverage': {'fraction_present': 1.}}
    def moran(data, **kwargs):
        calls['moran'] += 1
        return pd.DataFrame({'gene_id': data.var_names, 'moran': .1, 'status': 'computed'})
    monkeypatch.setattr(pathway_activity, 'compute_pathway_scores', auc)
    monkeypatch.setattr(spatial_autocorrelation, 'compute_moran', moran)
    monkeypatch.setattr(impact, 'compute_partitions', lambda data, **kw: {
        'labels': pd.Series(['0', '1'] * 6, index=data.obs_names), 'summary': {'n_units': len(data)}})
    w = impact.ImpactWorkflow(s, tmp_path, {'scopes': ['All'], 'min_window_units': 2,
        'n_window_draws': 3, 'region_n_bootstrap': 2, 'gene_sets': {'P': ['0']},
        'membership_same_units': True})
    through(w, 'integration')
    scores = (tmp_path / 'tables/pathway_svc_P.csv').read_bytes()
    membership = (tmp_path / 'tables/membership_All.csv').read_bytes()
    w.apply_parameters(dict(w.parameters, region_n_bootstrap=3))
    w.run_stage('regions')
    w.run_stage('integration')
    assert calls == {'auc': 2, 'moran': 2}
    assert (tmp_path / 'tables/pathway_svc_P.csv').read_bytes() == scores
    assert (tmp_path / 'tables/membership_All.csv').read_bytes() == membership
    assert 'tables/integrated_changed_unit_context_All.csv' in w.outputs


def test_raw_baseline_failure_does_not_block_native_label_state(tmp_path, monkeypatch):
    w = workflow(tmp_path)
    w.continue_on_error = True
    monkeypatch.setattr(w, 'stage_baseline', lambda: (_ for _ in ()).throw(RuntimeError('baseline failed')))
    through(w, 'integration')
    assert 'tables/state_All_region_windows.csv' in w.outputs
    assert 'tables/integrated_anatomy_conditionals_All.csv' in w.outputs
    assert [error['stage'] for error in w.stage_errors] == ['baseline']


def test_small_parameter_patch_preserves_other_explicit_parameters(tmp_path):
    w = workflow(tmp_path)
    through(w, 'diversity')
    changed = w.apply_parameters({'region_n_bootstrap': 3})
    assert w.parameters['scopes'] == ['All']
    assert w.parameters['min_window_units'] == 2
    assert set(changed) == {'regions', 'integration', 'figures'}
    w.parameters['region_n_bootstrap'] = 99
    for consumer in (w.result, w.render_figures):
        with pytest.raises(RuntimeError, match='apply_parameters'):
            consumer()


def test_first_region_run_preserves_pure_upstream_figures(tmp_path):
    w = workflow(tmp_path)
    through(w, 'diversity')
    retained = {'figures/window_support_grid_All.png': 'support',
                'figures/window_diversity_svc_All.png': 'diversity',
                'figures/svc_anatomy_context.png': 'anatomy'}
    for path, section in retained.items():
        w.outputs[path] = path
        w._figure_sections[path] = section
    w.run_stage('regions')
    assert all(path in w.outputs for path in retained)
