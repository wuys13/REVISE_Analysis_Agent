"""Scientific availability and navigation regressions in the saved-result reader."""
from html.parser import HTMLParser
import json
from pathlib import Path

from revise_analysis.reporting import render_report
from revise_analysis.reporting.report import _topic_unavailable


class _Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = []
        self.anchors = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if 'id' in attributes:
            self.ids.append(attributes['id'])
        href = attributes.get('href', '')
        if href.startswith('#'):
            self.anchors.append(href[1:])


def test_availability_preserves_scope_slug_and_global_side_dependencies():
    result = {'parameters': {'scopes': ['All', 'Fibroblast', 'Mono/Macro', 'T']},
              'unavailable': [
                  {'component': 'state_Mono_Macro', 'reason': 'no_stable_threshold'},
                  {'component': 'partition:raw', 'reason': 'unknown expression'},
                  {'component': 'moran:svc:T', 'reason': 'unknown expression'},
              ]}
    assert _topic_unavailable('state', 'Mono/Macro', result)
    for scope in ('All', 'Fibroblast', 'T'):
        assert not _topic_unavailable('state', scope, result)
    assert _topic_unavailable('rawbaseline', 'Fibroblast', result)
    assert not _topic_unavailable('genes', 'Fibroblast', result)
    assert _topic_unavailable('genes', 'T', result)


def test_unknown_region_is_not_zero_and_every_summary_anchor_resolves(tmp_path: Path):
    tables = tmp_path / 'tables'
    figures = tmp_path / 'figures'
    tables.mkdir()
    figures.mkdir()
    (tables / 'input_overview.csv').write_text('side,n_units,n_genes\nraw,10,5\nsvc,6,5\n')
    (tables / 'reconstruction_label_summary.csv').write_text('reconstruction_label,n_units\na,3\nb,3\n')
    (tables / 'integrated_region_anatomy_summary_All.csv').write_text(
        'region_kind,in_region,anatomy_region,n_units,denominator_units,region_definition_status\n'
        'state,,Tumor,6,6,unavailable\n')
    (figures / 'reconstruction_labels.png').write_bytes(b'presentation-only-fixture')
    outputs = {p.name: p.relative_to(tmp_path).as_posix()
               for directory in (tables, figures) for p in directory.iterdir()}
    record = {'analysis': 'reconstruction_impact', 'sample_id': 'fixture', 'status': 'partial',
              'parameters': {'scopes': ['All']}, 'outputs': outputs}
    (tmp_path / 'result.json').write_text(json.dumps(record))
    html = render_report(tmp_path).read_text()
    assert 'State 区域状态尚不可判定' in html
    assert 'State 区域含 0/' not in html
    assert 'scope-population-All' not in html
    links = _Links()
    links.feed(html)
    assert len(links.ids) == len(set(links.ids))
    assert set(links.anchors) <= set(links.ids)
