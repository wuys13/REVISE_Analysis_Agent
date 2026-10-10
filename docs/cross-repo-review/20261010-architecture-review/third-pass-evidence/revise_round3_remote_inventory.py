import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path('/inspire/hdd/project/sais-phenotype-fm/public/wuyushuai/zz/0_REVISE_developing/REVISE_code')
SAMPLE = ROOT / 'REVISE/output/CRC/P2_CRC_Xenium'
PARENT_SAMPLE = ROOT / 'output/CRC/P2_CRC_Xenium'
HASH_LIMIT = 64 * 1024 * 1024


def run(cmd, cwd=None):
    try:
        p = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           timeout=45)
        return {'returncode': p.returncode, 'stdout': p.stdout.strip(),
                'stderr': p.stderr.strip()}
    except Exception as e:
        return {'returncode': None, 'stdout': '', 'stderr': repr(e)}


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def file_record(path, root=ROOT, include_hash=True):
    try:
        st = path.stat()
    except OSError as e:
        return {'path': str(path), 'relative_to_root': str(path.relative_to(root)) if path.is_relative_to(root) else None,
                'exists': False, 'error': repr(e)}
    rec = {
        'path': str(path),
        'relative_to_root': str(path.relative_to(root)) if path.is_relative_to(root) else None,
        'exists': True,
        'is_file': path.is_file(),
        'size': st.st_size,
        'mtime_ns': st.st_mtime_ns,
    }
    if include_hash and path.is_file():
        if path.suffix.lower() == '.h5ad':
            rec['sha256'] = None
            rec['hash_status'] = 'skipped_h5ad'
        elif st.st_size > HASH_LIMIT:
            rec['sha256'] = None
            rec['hash_status'] = f'skipped_size_gt_{HASH_LIMIT}'
        else:
            try:
                rec['sha256'] = sha256(path)
                rec['hash_status'] = 'sha256'
            except OSError as e:
                rec['sha256'] = None
                rec['hash_status'] = f'error:{e!r}'
    return rec


def git_identity(path):
    if not path.exists():
        return {'path': str(path), 'exists': False}
    top = run(['git', 'rev-parse', '--show-toplevel'], path)
    if top['returncode'] != 0:
        return {'path': str(path), 'exists': True, 'is_git': False,
                'git_probe': top}
    head = run(['git', 'rev-parse', 'HEAD'], path)
    branch = run(['git', 'branch', '--show-current'], path)
    status = run(['git', 'status', '--short', '--branch'], path)
    upstream = run(['git', 'rev-parse', '--abbrev-ref', '--symbolic-full-name', '@{u}'], path)
    tracking_head = run(['git', 'rev-parse', '@{u}'], path)
    log1 = run(['git', 'log', '-1', '--format=%H%x09%ad%x09%s', '--date=iso-strict'], path)
    refs = {}
    b = branch['stdout']
    if b:
        live = run(['git', 'ls-remote', 'origin', f'refs/heads/{b}'], path)
        refs[f'origin/{b}'] = live['stdout']
        refs[f'origin/{b}_stderr'] = live['stderr'] if live['returncode'] != 0 else ''
    return {
        'path': str(path),
        'exists': True,
        'is_git': True,
        'show_toplevel': top['stdout'],
        'head': head['stdout'],
        'branch': b,
        'status': status['stdout'],
        'upstream': upstream['stdout'],
        'upstream_returncode': upstream['returncode'],
        'tracking_head': tracking_head['stdout'],
        'tracking_returncode': tracking_head['returncode'],
        'last_commit': log1['stdout'],
        'live_origin_refs': refs,
    }


def output_inventory(path):
    if not path.exists():
        return {'path': str(path), 'exists': False, 'file_count': 0, 'files': []}
    files = []
    for p in sorted(path.rglob('*')):
        if p.is_file():
            rel = p.relative_to(path)
            rels = str(rel)
            is_curated = rels.startswith('analysis/') or rels in {
                'README.md', 'summary.json', 'agent_input_payload.json',
                'report.html', 'local-curated-manifest.json'
            }
            rec = file_record(p, ROOT, include_hash=True)
            rec['relative_to_sample'] = rels
            rec['classification'] = 'curated_candidate' if is_curated else 'other'
            files.append(rec)
    return {
        'path': str(path),
        'exists': True,
        'file_count': len(files),
        'curated_candidate_count': sum(x['classification'] == 'curated_candidate' for x in files),
        'old_slot_dirs': {
            'analysis/00_benchmarks': (path / 'analysis/00_benchmarks').exists(),
            'analysis/01_discoveries': (path / 'analysis/01_discoveries').exists(),
            'analysis/00_benchmark': (path / 'analysis/00_benchmark').exists(),
        },
        'files': files,
    }


def path_exists_records(paths):
    out = []
    for p in paths:
        p = Path(p)
        rec = file_record(p, ROOT, include_hash=True)
        if p.is_dir():
            rec['is_dir'] = True
        out.append(rec)
    return out


def scan_named(root, names):
    found = []
    if not root.exists():
        return found
    for p in sorted(root.rglob('*')):
        if p.is_file() and any(tok in p.name.lower() for tok in names):
            found.append(file_record(p, ROOT, include_hash=True))
    return found

# Explicit deployment contract paths.
script_doc_paths = [
    ROOT / 'REVISE_Analysis_Agent/scripts/run_unbiased_all_subcluster_analysis.py',
    ROOT / 'REVISE_Analysis_Agent/scripts/organize_three_tier_analysis.py',
    ROOT / 'REVISE_Analysis_Agent/scripts/prepare_agent_payload.py',
    ROOT / 'REVISE_Analysis_Agent/scripts/build_sample_analysis_delivery.py',
    ROOT / 'REVISE_Analysis_Agent/docs/ANALYSIS_OUTPUT_PARADIGM.md',
    ROOT / 'REVISE_Analysis_Agent/docs/AGENT_GENERATION_SPEC.md',
    ROOT / 'REVISE_Analysis_Agent/docs/EXECUTION_FEEDBACK.md',
    ROOT / 'REVISE_Analysis_Agent/docs/cross-repo-review/20261010-architecture-review/local-curated-manifest.json',
    ROOT / 'REVISE_Analysis_Agent/docs/cross-repo-review/20261010-architecture-review/remote-observation-index.txt',
    ROOT / 'REVISE_Analysis_Agent/docs/cross-repo-review/20261010-architecture-review/remote-single-file-supplement.json',
]

source_paths = [
    ROOT / 'REVISE/reproduce/case/CCI/cellphonedb.zip',
    ROOT / 'REVISE/reproduce/case/pathway/HYPOXIA_related_pathways.gmt',
    ROOT / 'REVISE/reproduce/case/pathway/KEGG_2021_Human.gmt',
    ROOT / 'REVISE/reproduce/case/pathway/h.all.v2025.1.Hs.symbols.gmt',
    ROOT / 'REVISE/reproduce/case/pathway/GO_Biological_Process_2025.gmt',
    ROOT / 'REVISE/revise/reference_preparation/SOURCE_PROVENANCE.md',
    ROOT / 'REVISE/reproduce/case/P2CRC_Xenium.h5ad',
    ROOT / 'REVISE/reproduce/case/adata_sc_all_reanno.h5ad',
]

repos = {
    'workspace_root': ROOT,
    'REVISE': ROOT / 'REVISE',
    'REVISE_Analysis_Agent': ROOT / 'REVISE_Analysis_Agent',
}

ps = run(['ps', '-eo', 'pid=,etime=,args='])
keywords = ('reconstruct.py', 'run_unbiased_all_subcluster_analysis.py',
            'organize_three_tier_analysis.py', 'prepare_agent_payload.py',
            'build_sample_analysis_delivery.py', 'P2_CRC_Xenium')
processes = []
for line in ps['stdout'].splitlines():
    if any(k in line for k in keywords):
        processes.append(line.strip())

result = {
    'schema': 'revise.remote_inventory.round3.v1',
    'observed_at_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
    'notebook_claim_checked': {'name': 'cpu-2', 'id': '58bfcc6c-e848-44a7-a227-2ee82d25da01'},
    'remote_root': str(ROOT),
    'git': {k: git_identity(v) for k, v in repos.items()},
    'script_docs': path_exists_records(script_doc_paths),
    'source_inputs_and_provenance': path_exists_records(source_paths),
    'sample_output': output_inventory(SAMPLE),
    'parent_sibling_output': output_inventory(PARENT_SAMPLE),
    'manifest_source_receipt_scan': scan_named(ROOT, ('manifest', 'source', 'provenance', 'receipt', 'index')),
    'old_slot_scan': {
        '00_benchmarks': [str(p) for p in ROOT.rglob('00_benchmarks')],
        '01_discoveries': [str(p) for p in ROOT.rglob('01_discoveries')],
    },
    'matching_processes': processes,
    'ps_probe': {'returncode': ps['returncode'], 'stderr': ps['stderr']},
}
print(json.dumps(result, ensure_ascii=False, sort_keys=True))
