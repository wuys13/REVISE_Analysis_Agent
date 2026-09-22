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
    changed = dict(value['parameters'], sample_n_units=29_999)
    with pytest.raises(ValueError, match='frozen P2 full baseline'):
        module.require_frozen_parameters(changed, repo)


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
