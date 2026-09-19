"""Native molecular results and aggregation of fixed, saved unit scores."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from ._shared import coordinates, deterministic_subset, unavailable
from .impact_tables import aggregate_scores
from revise_analysis.methods.errors import PrerequisiteUnavailable


def _slug(value):
    return ''.join(c if c.isalnum() else '_' for c in str(value)).strip('_') or 'scope'


def _cohort(workflow, side, scope):
    data = getattr(workflow.sample, side)
    if scope == 'All':
        return data
    if workflow.sample.broad_key not in data.obs:
        raise PrerequisiteUnavailable(f'{side}: broad-label column is unavailable for {scope}')
    labels = workflow.sample.labels(side)
    return data[labels.index[labels == scope]]


def _input_reason(adata, spatial_key):
    if adata.n_obs == 0 or adata.n_vars == 0 or adata.X is None:
        return 'nonempty expression cohort and gene axis required'
    values = adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X)
    if not np.isfinite(values).all() or (values < 0).any():
        return 'expression must be finite and nonnegative'
    try:
        points = coordinates(adata, spatial_key)
    except (KeyError, ValueError) as exc:
        return str(exc)
    if not np.isfinite(points.to_numpy()).all():
        return 'coordinates must be finite'
    return None


def run_molecular(workflow):
    """Score a fixed union cohort once per side/program; never score a Region."""
    from .pathway_activity import compute_pathway_scores, resolve_gene_sets
    from .spatial_autocorrelation import compute_moran, record_uncomputed

    sample, params = workflow.sample, workflow.parameters
    moran, score_tables = {}, []
    workflow.state['moran'], workflow.state['program_scores'] = moran, score_tables
    for scope in params['scopes']:
        for side in ('raw', 'svc'):
            component = f'moran:{side}:{scope}'
            if reason := sample.expression_unavailable(side):
                unavailable(component, reason, workflow.missing)
                continue
            try:
                cohort = deterministic_subset(_cohort(workflow, side, scope), params['sample_n_units'],
                                              params['random_state'] + (side == 'svc'))
            except PrerequisiteUnavailable as exc:
                unavailable(component, str(exc), workflow.missing)
                continue
            if reason := _input_reason(cohort, sample.spatial_key):
                unavailable(component, reason, workflow.missing)
                continue
            try:
                table = compute_moran(cohort, spatial_key=sample.spatial_key,
                                      n_neighbors=params['moran_n_neighbors'],
                                      transformation_state=sample.expression(side)['scale']).copy()
            except ImportError as exc:
                unavailable(component, str(exc), workflow.missing)
                continue
            table['side'], table['scope'] = side, scope
            table['n_units'], table['n_native_genes'] = cohort.n_obs, cohort.n_vars
            moran[side, scope] = table
            workflow._save_table(table, f'tables/moran_{side}_{_slug(scope)}.csv')
            record_uncomputed(table, f'{side}:{scope}', workflow.missing)
        raw, svc = moran.get(('raw', scope)), moran.get(('svc', scope))
        if raw is None or svc is None:
            unavailable(f'moran_shared_genes:{scope}', 'both native Moran results required', workflow.missing)
            continue
        def index_genes(frame):
            key = next((key for key in ('gene_id', 'gene') if key in frame), None)
            return frame.set_index(key) if key else frame
        left, right = index_genes(raw), index_genes(svc)
        common = left.index.intersection(right.index)
        table = pd.DataFrame({'gene_id': common,
                              'raw_moran': left.loc[common, 'moran'].to_numpy(),
                              'svc_moran': right.loc[common, 'moran'].to_numpy()})
        table['delta_moran'] = table.svc_moran - table.raw_moran
        table['comparison_available'] = np.isfinite(table.raw_moran) & np.isfinite(table.svc_moran)
        table['n_computable_shared_genes'] = int(table.comparison_available.sum())
        for side, source in (('raw', left), ('svc', right)):
            table[f'{side}_n_units'] = int(source.n_units.iloc[0])
            table[f'{side}_n_native_genes'] = len(source)
            if 'status' in source:
                table[f'{side}_status'] = source.loc[common, 'status'].to_numpy()
        table['n_common_genes'], table['scope'] = len(common), scope
        workflow._save_table(table, f'tables/moran_shared_genes_{_slug(scope)}.csv')
        workflow._save_json({'scope': scope, 'n_common_genes': len(common),
                             'n_computable_shared_genes': int(table.comparison_available.sum()),
                             'raw_n_native_genes': len(left), 'svc_n_native_genes': len(right),
                             'raw_n_units': int(left.n_units.iloc[0]),
                             'svc_n_units': int(right.n_units.iloc[0])},
                            f'tables/moran_shared_genes_{_slug(scope)}_metadata.json')
        if not table.comparison_available.any():
            unavailable(f'moran_shared_genes:{scope}', 'no jointly computable shared genes', workflow.missing)

    gene_sets, error = resolve_gene_sets(params)
    if error:
        unavailable('pathway_resource', error, workflow.missing)
        return
    for side in ('raw', 'svc'):
        if reason := sample.expression_unavailable(side):
            for name in gene_sets:
                unavailable(f'pathway:{side}:{name}', reason, workflow.missing)
            continue
        full = getattr(sample, side)
        if sample.broad_key in full.obs:
            broad = sample.labels(side)
        elif 'All' in params['scopes']:
            broad = pd.Series('Unknown', index=full.obs_names)
        else:
            unavailable(f'pathway:{side}', 'broad-label column unavailable for scoring scopes', workflow.missing)
            continue
        ids = full.obs_names if 'All' in params['scopes'] else broad.index[broad.isin(params['scopes'])]
        cohort = deterministic_subset(full[ids], params['sample_n_units'], params['random_state'] + (side == 'svc'))
        if reason := _input_reason(cohort, sample.spatial_key):
            unavailable(f'pathway:{side}', reason, workflow.missing)
            continue
        for name, genes in gene_sets.items():
            try:
                result = compute_pathway_scores(cohort, genes=list(genes), score_name=name,
                                                auc_threshold=params['pathway_auc_threshold'], seed=params['random_state'])
            except (ImportError, PrerequisiteUnavailable) as exc:
                unavailable(f'pathway:{side}:{name}', str(exc), workflow.missing)
                continue
            # Unexpected provider ValueError/RuntimeError deliberately propagate.
            frame = coordinates(cohort, sample.spatial_key).join(result['scores'].rename('score'))
            frame['side'], frame['program'] = side, name
            frame['scope'] = '|'.join(params['scopes'])
            frame['broad_label'] = broad.reindex(frame.index)
            frame['rank_cutoff'] = result['metadata'].get('rank_cutoff')
            frame['coverage_fraction'] = result['coverage'].get('fraction_present')
            score_tables.append(frame)
            relative = f'tables/pathway_{side}_{_slug(name)}'
            workflow._save_table(frame, f'{relative}.csv')
            workflow._save_json({'coverage': result['coverage'], 'metadata': result['metadata'],
                                 'scoring_scopes': params['scopes'], 'n_eligible_units': len(ids),
                                 'n_scored_units': cohort.n_obs, 'gene_axis': cohort.var_names.tolist(),
                                 'expression': sample.expression(side)}, f'{relative}_metadata.json')


def aggregate_programs(workflow):
    """Join fixed scores onto full eligible populations, retaining missing values."""
    from revise_analysis.methods.regions import assign_square_windows
    for scores in workflow.state.get('program_scores', []):
        side, program = str(scores.side.iloc[0]), str(scores.program.iloc[0])
        for scope in workflow.parameters['scopes']:
            try:
                cohort = _cohort(workflow, side, scope)
            except PrerequisiteUnavailable as exc:
                unavailable(f'program_aggregation:{side}:{scope}', str(exc), workflow.missing)
                continue
            if not cohort.n_obs:
                continue
            joined = coordinates(cohort, workflow.sample.spatial_key).join(scores[['score']])
            joined['side'], joined['scope'], joined['program'] = side, scope, program
            points = workflow.state.get('point_anatomy', {}).get(side)
            joined['anatomy_region'] = (points.anatomy_region.reindex(joined.index)
                                        if points is not None else 'Unknown')
            joined['score_status'] = np.where(joined.score.notna(), 'scored', 'not_in_scoring_cohort')
            if 'origin' in workflow.state:
                assignments = assign_square_windows(joined[['x', 'y']],
                    window_side_length=workflow.parameters['parent_window_side_microns'][scope] / workflow.state['microns_per_coordinate'],
                    origin=workflow.state['origin'])
                joined['window_id'] = assignments.window_id
                for kind in ('state', 'gain'):
                    region = (workflow.state.get('region_tables', {}).get(('state', scope)) if kind == 'state'
                              else workflow.state.get('gains', {}).get(scope))
                    joined[f'{kind}_region'] = (joined.window_id.map(region.set_index('window_id').in_region)
                                                if region is not None else pd.NA)
            prefix = f'tables/integrated_program_{side}_{_slug(scope)}_{_slug(program)}'
            workflow._save_table(joined, f'{prefix}_units.csv')
            groupings = [('anatomy', ['side', 'scope', 'anatomy_region'])]
            if 'window_id' in joined:
                groupings += [('window', ['side', 'scope', 'window_id']),
                              ('region', ['side', 'scope', 'state_region', 'gain_region', 'anatomy_region'])]
            for suffix, keys in groupings:
                summary = aggregate_scores(joined, group_columns=keys)
                summary['evidence_table'] = f'{prefix}_units.csv'
                workflow._save_table(summary, f'{prefix}_{suffix}.csv')
