import importlib.util
from pathlib import Path

import pytest


def launcher():
    path = Path(__file__).parents[1] / 'scripts' / 'run_p2_impact.py'
    spec = importlib.util.spec_from_file_location('run_p2_impact', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_formal_config_is_impact_only_and_matches_frozen_parameters():
    module = launcher()
    repo = Path(__file__).parents[1]
    value = module.validate_project(repo / module.EXPECTED_CONFIG, repo, require_inputs=False,
                                    iteration='iter-001', producer_iteration='iter-001')
    assert list(value['config']['analyses']) == ['reconstruction_impact']
    assert value['parameters']['scopes'] == ['All', 'Fibroblast', 'Mono_Macro', 'T']
    for changed in (
        dict(value['parameters'], sample_n_units=29_999),
        dict(value['parameters'], random_state=7),
        dict(value['parameters'], raw_k_control=True),
        dict(value['parameters'], unexpected_parameter=True),
    ):
        with pytest.raises(ValueError, match='frozen P2 full baseline'):
            module.require_frozen_parameters(changed, repo)


def test_raw_k_control_diagnostic_only_allows_its_single_parameter_difference(tmp_path):
    module = launcher()
    repo = Path(__file__).parents[1]
    config = repo / module.DIAGNOSTICS['raw-k-control']['config']
    assert module.selected_config_path(repo, None, None) == repo / module.EXPECTED_CONFIG
    assert module.selected_config_path(repo, None, 'raw-k-control') == config
    value = module.validate_project(config, repo, require_inputs=False,
                                    iteration='iter-002', producer_iteration='iter-001',
                                    diagnostic='raw-k-control')
    assert value['output_dir'] == repo / 'output/p2-xenium-impact-20260922/iter-002'
    assert value['parameters']['raw_k_control'] is True
    assert value['diagnostic']['allowed_parameter_difference'] == {
        'raw_k_control': {'frozen_baseline': False, 'diagnostic': True},
    }
    with pytest.raises(ValueError, match='only'):
        module.require_frozen_parameters(dict(value['parameters'], raw_k_control=1), repo,
                                         diagnostic='raw-k-control')
    changed_seed = dict(value['parameters'], random_state=7)
    with pytest.raises(ValueError, match='only'):
        module.require_frozen_parameters(changed_seed, repo, diagnostic='raw-k-control')
    with pytest.raises(ValueError, match='analysis iteration'):
        module.validate_project(config, repo, require_inputs=False,
                                iteration='iter-001', producer_iteration='iter-001',
                                diagnostic='raw-k-control')
    with pytest.raises(ValueError, match='producer iteration'):
        module.validate_project(config, repo, require_inputs=False,
                                iteration='iter-002', producer_iteration='iter-002',
                                diagnostic='raw-k-control')
    cases = (
        ('scopes', '    - T\n', '    - Other\n', 'Impact scopes'),
        ('output', 'output_dir: ../output/p2-xenium-impact-20260922/iter-002',
         'output_dir: ../output/p2-xenium-impact-20260922/iter-003', 'output_dir must resolve'),
        ('input', 'results/p2-xenium-impact-20260922/iter-001',
         'results/p2-xenium-impact-20260922/iter-002', 'sample must resolve'),
    )
    source = config.read_text(encoding='utf-8')
    for name, old, new, error in cases:
        changed = repo / 'configs' / f'_test_raw_k_control_{tmp_path.name}_{name}.yaml'
        changed.write_text(source.replace(old, new, 1), encoding='utf-8')
        try:
            with pytest.raises(ValueError, match=error):
                module.validate_project(changed, repo, require_inputs=False,
                                        iteration='iter-002', producer_iteration='iter-001',
                                        diagnostic='raw-k-control')
        finally:
            changed.unlink()


def test_raw_k_control_diagnostic_receipt_is_explicit_and_non_resumable():
    module = launcher()
    repo = Path(__file__).parents[1]
    value = module.validate_project(repo / module.DIAGNOSTICS['raw-k-control']['config'], repo,
                                    require_inputs=False, iteration='iter-002',
                                    producer_iteration='iter-001', diagnostic='raw-k-control')
    args = type('Args', (), {
        'iteration': 'iter-002', 'producer_iteration': 'iter-001',
        'diagnostic': 'raw-k-control', 'execution_location': 'test', 'qz_notebook': None,
    })()
    receipt = module.base_receipt(args, value, {'commit': 'test'})
    assert receipt['diagnostic'] == {
        'name': 'raw-k-control',
        'frozen_baseline_config': str((repo / module.BASELINE_CONFIG).resolve()),
        'allowed_parameter_difference': {
            'raw_k_control': {'frozen_baseline': False, 'diagnostic': True},
        },
        'execution_scope': 'full Impact run; resume is unsupported',
    }
    assert 'receipt or Impact result' in receipt['resubmission_guard']
    with pytest.raises(SystemExit):
        module.main(['preflight', '--resume'])


def test_formal_run_requires_clean_commit_and_complete_artifact_registration():
    module = launcher()
    with pytest.raises(RuntimeError, match='clean committed worktree'):
        module.require_clean_worktree({'dirty': True})
    assert module.execution_status('partial', []) == 'execution_complete'
    assert module.execution_status('succeeded', ['missing.csv']) == 'artifact_validation_failed'


def test_analysis_and_producer_iterations_are_independent_safe_segments():
    module = launcher()
    repo = Path('/workspace/REVISE_Analysis_Agent')
    output, sample = module.expected_run_paths(repo, 'iter-002', 'iter-001')
    assert output == Path('/workspace/REVISE_Analysis_Agent/output/p2-xenium-impact-20260922/iter-002')
    assert sample == Path('/workspace/REVISE/results/p2-xenium-impact-20260922/iter-001/P2CRC_Xenium/random/sample.yaml')
    with pytest.raises(ValueError, match='iter-NNN'):
        module.expected_run_paths(repo, '../iter-002', 'iter-001')
