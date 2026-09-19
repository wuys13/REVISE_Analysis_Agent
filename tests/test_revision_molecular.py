from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from revise_analysis.io import Sample
from revise_analysis.analyses.impact_molecular import run_molecular, aggregate_programs


def make_workflow():
    raw = AnnData(np.ones((8, 4)), obs=pd.DataFrame({'Level1': ['Fibroblast'] * 4 + ['T'] * 4}, index=[f'u{i}' for i in range(8)]), var=pd.DataFrame(index=list('abcd')))
    raw.obsm['spatial'] = np.c_[np.arange(8), np.zeros(8)]
    sample = Sample('test', raw, raw.copy(), {'expression': {s: {'identity': 'test'} for s in ('raw', 'svc')}}, Path('sample.yaml'))
    params = {'scopes': ['All', 'Fibroblast'], 'sample_n_units': 3, 'random_state': 42,
              'moran_n_neighbors': 2, 'gene_sets': {'P': ['a']}, 'pathway_auc_threshold': .05,
              'parent_window_side_microns': {'All': 2, 'Fibroblast': 2}}
    saved = {}
    workflow = SimpleNamespace(sample=sample, parameters=params, missing=[],
                               state={'origin': (0, 0), 'microns_per_coordinate': 1})
    workflow._save_table = lambda frame, name: saved.update({name: frame.copy()})
    workflow._save_json = lambda obj, name: saved.update({name: obj.copy()})
    return workflow, saved


def test_fixed_program_cohort_scored_once_and_aggregated_with_full_denominators(monkeypatch):
    from revise_analysis.analyses import pathway_activity, spatial_autocorrelation
    calls = []
    def scores(adata, **kwargs):
        calls.append(adata.obs_names.tolist())
        return {'scores': pd.Series(1., index=adata.obs_names), 'metadata': {'rank_cutoff': 2}, 'coverage': {'fraction_present': 1.}}
    monkeypatch.setattr(pathway_activity, 'compute_pathway_scores', scores)
    monkeypatch.setattr(spatial_autocorrelation, 'compute_moran', lambda adata, **kw: pd.DataFrame({'gene_id': adata.var_names, 'moran': .1, 'status': 'computed'}))
    workflow, saved = make_workflow()
    run_molecular(workflow)
    assert len(calls) == 2  # one per side, not per scope
    assert all(len(ids) == 3 for ids in calls)
    aggregate_programs(workflow)
    assert len(calls) == 2
    table = saved['tables/integrated_program_svc_All_P_anatomy.csv']
    assert table.n_units.sum() == 8
    assert table.n_scored.sum() == 3
    units = saved['tables/integrated_program_svc_All_P_units.csv']
    assert units.score.isna().sum() == 5
    assert set(saved['tables/moran_shared_genes_All.csv'].raw_n_units) == {3}
    # Changing region definitions consumes the same score table, no new scoring.
    workflow.state['region_tables'] = {('state', 'All'): pd.DataFrame({'window_id': ['0_0', '1_0'], 'in_region': [True, False]})}
    aggregate_programs(workflow)
    assert len(calls) == 2


def test_molecular_provider_execution_errors_propagate_with_prior_native_tables(monkeypatch):
    from revise_analysis.analyses import pathway_activity, spatial_autocorrelation
    monkeypatch.setattr(spatial_autocorrelation, 'compute_moran', lambda adata, **kw: pd.DataFrame({'gene_id': adata.var_names, 'moran': .1, 'status': 'computed'}))
    def broken(*args, **kwargs):
        raise ValueError('unexpected provider error')
    monkeypatch.setattr(pathway_activity, 'compute_pathway_scores', broken)
    workflow, saved = make_workflow()
    with pytest.raises(ValueError, match='unexpected provider error'):
        run_molecular(workflow)
    assert 'tables/moran_raw_All.csv' in saved
    assert not any('unexpected provider' in item['reason'] for item in workflow.missing)
