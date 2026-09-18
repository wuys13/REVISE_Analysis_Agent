"""Command-line access to the same Python execution functions."""
import argparse
import json
from pathlib import Path


def main(argv=None):
    parser = argparse.ArgumentParser(prog="revise-analysis")
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("run", help="Run one analysis")
    one.add_argument("--sample", required=True)
    one.add_argument("--analysis", required=True)
    one.add_argument("--output", required=True)
    one.add_argument("--parameters", help="JSON object of analysis parameters")
    batch = commands.add_parser("batch", help="Run a project YAML")
    batch.add_argument("--config", required=True)
    report = commands.add_parser("report", help="Re-render saved results without scientific computation")
    report.add_argument("result_dir")
    args = parser.parse_args(argv)
    try:
        if args.command == "run":
            from .runner import run_analysis
            result = run_analysis(args.sample, args.analysis, args.output, json.loads(args.parameters) if args.parameters else {})
        elif args.command == "batch":
            from .batch import run_batch
            result = run_batch(args.config)
        else:
            from .reporting import render_report
            render_report(Path(args.result_dir))
            return 0
    except Exception as exc:
        parser.exit(2, f"{type(exc).__name__}: {exc}\n")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "succeeded" else 1


if __name__ == "__main__":
    raise SystemExit(main())
