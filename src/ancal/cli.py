"""`ancal` command-line entry point."""

from __future__ import annotations

import argparse
import importlib
import shutil
import subprocess
import sys

from . import budget, status
from .settings import load_settings

REQUIRED_TOOLS = ("bcftools", "tabix", "bgzip", "samtools", "snakemake")
REQUIRED_MODULES = ("polars", "pyarrow", "pandas", "pydantic", "yaml", "requests", "fsspec")


def _tool_version(tool: str) -> str:
    out = subprocess.run([tool, "--version"], capture_output=True, text=True)
    return (out.stdout or out.stderr).splitlines()[0].strip() if (out.stdout or out.stderr) else "?"


def cmd_budget_check(args) -> int:
    s = load_settings(scope_override=args.scope)
    u = budget.usage(s)
    print(f"scope={s.scope}  budget={s.local_budget_gb} GB  used={u.used_mb:.1f} MB")
    for d, b in sorted(u.by_dir.items()):
        print(f"  {d:<32} {b / 1024**2:10.2f} MB")
    if not u.ok:
        print("OVER BUDGET - see SCALING_UP.md", file=sys.stderr)
        return 1
    return 0


def cmd_check_env(args) -> int:
    s = load_settings()
    checks = []
    for tool in REQUIRED_TOOLS:
        path = shutil.which(tool)
        checks.append(status.Check(f"tool:{tool}", path is not None, _tool_version(tool) if path else "not on PATH"))
    for mod in REQUIRED_MODULES:
        try:
            m = importlib.import_module(mod)
            checks.append(status.Check(f"python:{mod}", True, getattr(m, "__version__", "")))
        except ImportError as e:
            checks.append(status.Check(f"python:{mod}", False, str(e)))
    u = budget.usage(s)
    checks.append(status.Check("disk_budget", u.ok, f"{u.used_mb:.1f} MB of {s.local_budget_gb} GB"))
    if args.pytest_ok is not None:
        checks.append(status.Check("pytest", args.pytest_ok == "1", "unit tests"))

    st = status.write(s, "stage0", checks, notes="Environment and repository scaffold (item 1).")
    for c in checks:
        print(f"  [{'ok' if c.ok else 'FAIL'}] {c.name:<20} {c.detail}")
    print(f"stage0: {st.status}")
    return 0 if st.status == "pass" else 1


def cmd_status(args) -> int:
    s = load_settings()
    rows = status.read_all(s)
    if not rows:
        print("No stage has run yet.")
        return 0
    for r in rows:
        extra = f"  ({r['blocked_reason']})" if r.get("blocked_reason") else ""
        print(f"{r['stage']:<8} {r['status']:<8} scope={r['scope']:<10} disk={r['disk_used_mb']} MB{extra}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="ancal")
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("budget-check", help="Report local disk usage against the budget")
    b.add_argument("--scope", default=None)
    b.set_defaults(func=cmd_budget_check)

    e = sub.add_parser("check-env", help="Stage 0 gate: tools, imports, budget")
    e.add_argument("--pytest-ok", choices=("0", "1"), default=None)
    e.set_defaults(func=cmd_check_env)

    st = sub.add_parser("status", help="Show the status of every stage run so far")
    st.set_defaults(func=cmd_status)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
