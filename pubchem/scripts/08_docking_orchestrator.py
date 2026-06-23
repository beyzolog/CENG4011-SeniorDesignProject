#!/usr/bin/env python3
"""
Aşama 8 — Molecular Docking Orchestrator

Upstream : 07_toxicity_and_grouping.py → final_1000_cids.tmp / *_final_1000.csv
Engine   : AutoDock Vina via docking/core/

Usage:
  TARGET_GENE=CTNNB1 EXPERIMENT_SUFFIX=06 python3 08_docking_orchestrator.py
  TARGET_GENE=MYC EXPERIMENT_SUFFIX=07 python3 08_docking_orchestrator.py --limit 100 --cpu 16
  TARGET_GENE=both EXPERIMENT_SUFFIX=06 python3 08_docking_orchestrator.py --parallel
  DOCKING_LIMIT=50 TARGET_GENE=CTNNB1 python3 08_docking_orchestrator.py --dry-run
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
PUBCHEM_ROOT = SCRIPT_DIR.parent
if str(PUBCHEM_ROOT) not in sys.path:
    sys.path.insert(0, str(PUBCHEM_ROOT))

from docking.core.champion_loader import (
    DEFAULT_DOCKING_LIMIT,
    load_docking_cids,
    resolve_docking_source,
)
from docking.core.paths import RunPaths, ensure_run_dirs

DEFAULT_GENES = ("MYC", "CTNNB1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stage 8: Vina docking orchestrator")
    parser.add_argument(
        "--gene",
        default=os.environ.get("TARGET_GENE", "CTNNB1"),
        choices=[*DEFAULT_GENES, "both"],
    )
    parser.add_argument(
        "--exp",
        default=os.environ.get("EXPERIMENT_SUFFIX", "06"),
        help="Experiment suffix matching screening/predictions_exp_{NN}/",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=int(os.environ.get("DOCKING_LIMIT", str(DEFAULT_DOCKING_LIMIT))),
        help=f"Number of top-scoring ligands to dock (default: {DEFAULT_DOCKING_LIMIT})",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        default=int(os.environ.get("DOCKING_CPU", "0")),
    )
    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Run MYC+CTNNB1 as separate processes",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve CIDs and paths without running Vina",
    )
    parser.add_argument(
        "--analyze",
        action="store_true",
        help="Run 09_docking_analysis.py after docking completes",
    )
    return parser.parse_args()


def orchestrate(gene: str, exp: str, limit: int, cpu: int, dry_run: bool) -> Path:
    gene = gene.upper()
    cids = load_docking_cids(gene, exp, limit=limit)
    source = resolve_docking_source(gene, exp)
    paths = RunPaths.for_gene(gene, exp)
    ensure_run_dirs(paths)

    if dry_run:
        preview = cids[:5]
        suffix = "..." if len(cids) > 5 else ""
        print(f"[DRY-RUN] {gene} exp{exp} → N={len(cids)} CIDs from {source.name}")
        print(f"          first CIDs: {preview}{suffix}")
        print(f"          run_dir={paths.run_dir}")
        print(f"          summary={paths.summary_json}")
        return paths.summary_json

    from docking.core.pipeline import run_gene_docking

    run_gene_docking(gene=gene, cids=cids, paths=paths, cpu=cpu)
    return paths.summary_json


def main() -> None:
    args = parse_args()
    genes = DEFAULT_GENES if args.gene == "both" else (args.gene.upper(),)

    if args.parallel and len(genes) > 1 and not args.dry_run:
        from multiprocessing import Process

        procs = [
            Process(
                target=orchestrate,
                args=(g, args.exp, args.limit, args.cpu, args.dry_run),
            )
            for g in genes
        ]
        for proc in procs:
            proc.start()
        for proc in procs:
            proc.join()
    else:
        for gene in genes:
            orchestrate(gene, args.exp, args.limit, args.cpu, args.dry_run)

    if args.analyze and not args.dry_run:
        import subprocess

        analysis_script = SCRIPT_DIR / "09_docking_analysis.py"
        subprocess.run(
            [
                sys.executable,
                str(analysis_script),
                "--exp",
                args.exp,
                "--limit",
                str(args.limit),
            ],
            check=False,
        )


if __name__ == "__main__":
    main()
