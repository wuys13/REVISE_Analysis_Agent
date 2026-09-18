import json
from pathlib import Path

from revise_analysis.reporting import render_report


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
