#!/usr/bin/env python3
"""Regression validation for updated_docking parser against legacy exp06/07 outputs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

PUBCHEM_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PUBCHEM_ROOT / "scripts"
LEGACY_LIMIT = 5


def main() -> int:
    checks = [
        (
            "CTNNB1 exp06 (limit=5)",
            [
                sys.executable,
                str(SCRIPTS / "09_docking_analysis.py"),
                "--exp",
                "06",
                "--limit",
                str(LEGACY_LIMIT),
                "--legacy-logs",
                str(PUBCHEM_ROOT / "docking" / "ctnnb1_docking" / "logs"),
                "--legacy-gene",
                "CTNNB1",
            ],
        ),
        (
            "MYC exp07 (limit=5, legacy exp06 filename)",
            [
                sys.executable,
                str(SCRIPTS / "09_docking_analysis.py"),
                "--exp",
                "07",
                "--limit",
                str(LEGACY_LIMIT),
                "--legacy-logs",
                str(PUBCHEM_ROOT / "docking" / "myc_docking" / "logs"),
                "--legacy-gene",
                "MYC",
            ],
        ),
    ]

    print("=" * 70)
    print("  updated_docking seed validation (legacy top-5 regression)")
    print("=" * 70)
    failed = 0
    for label, cmd in checks:
        print(f"\n>>> {label}")
        result = subprocess.run(cmd, cwd=SCRIPTS)
        if result.returncode != 0:
            failed += 1
            print(f"  FAILED: {label}")
        else:
            print(f"  OK: {label}")

    print("\n" + "=" * 70)
    if failed:
        print(f"  {failed} check(s) failed")
        return 1
    print("  All seed validation checks passed")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
