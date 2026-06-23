#!/usr/bin/env python3
"""
Aşama 9 — Docking Analysis & Reporting

Aggregates docking run summaries into cross-gene leaderboard,
per-gene markdown reports, and a reproducibility manifest.

Usage:
  python3 09_docking_analysis.py --exp 06 --limit 100
  python3 09_docking_analysis.py --exp 07 --genes MYC --limit 100
  python3 09_docking_analysis.py --exp 06 --legacy-logs ../docking/ctnnb1_docking/logs --limit 5
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PUBCHEM_ROOT = SCRIPT_DIR.parent
if str(PUBCHEM_ROOT) not in sys.path:
    sys.path.insert(0, str(PUBCHEM_ROOT))

from docking.core.champion_loader import (
    DEFAULT_DOCKING_LIMIT,
    load_docking_cids,
    load_prediction_scores,
)
from docking.core.gene_config import load_gene_config
from docking.core.paths import PACKAGE_ROOT, RunPaths
from docking.core.result_parser import parse_all_results, save_summary_json

DEFAULT_GENES = ("MYC", "CTNNB1")
REPORTS_DIR = PACKAGE_ROOT / "reports"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stage 9: docking analysis")
    parser.add_argument("--exp", default="06", help="Experiment suffix")
    parser.add_argument(
        "--limit",
        type=int,
        default=DEFAULT_DOCKING_LIMIT,
        help=f"Expected docking candidate count (default: {DEFAULT_DOCKING_LIMIT})",
    )
    parser.add_argument(
        "--genes",
        nargs="+",
        default=list(DEFAULT_GENES),
        help="Genes to include",
    )
    parser.add_argument(
        "--legacy-logs",
        type=Path,
        default=None,
        help="Optional legacy logs dir for regression parsing",
    )
    parser.add_argument(
        "--legacy-gene",
        default=None,
        help="Gene name when using --legacy-logs",
    )
    return parser.parse_args()


def load_summary_records(gene: str, exp: str, limit: int) -> list[dict]:
    paths = RunPaths.for_gene(gene, exp)
    if paths.summary_json.exists():
        with open(paths.summary_json, encoding="utf-8") as fh:
            return json.load(fh)
    if not any(paths.logs_dir.glob("*_cmd.log")):
        return []
    cids = load_docking_cids(gene, exp, limit=limit)
    return parse_all_results(paths, gene.upper(), cids)


def build_leaderboard_rows(gene: str, exp: str, records: list[dict]) -> list[dict]:
    gene = gene.upper()
    scores = load_prediction_scores(gene, exp, [r["cid"] for r in records])
    rows = []
    for rank, record in enumerate(records, start=1):
        rows.append(
            {
                "gene": gene,
                "exp": exp,
                "rank": rank,
                "cid": record["cid"],
                "prediction_score": scores.get(record["cid"]),
                "best_affinity_kcal_mol": record.get("best_affinity"),
                "best_mode_rmsd_ub": record.get("best_mode_rmsd"),
                "num_poses": record.get("num_poses", 0),
            }
        )
    return rows


def write_gene_report(gene: str, exp: str, records: list[dict], out_path: Path) -> None:
    gene = gene.upper()
    lines = [
        f"# Docking Report — {gene} exp_{exp}",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        f"Candidates: {len(records)}",
        "",
        "| Rank | CID | Affinity (kcal/mol) | RMSD ub | Poses |",
        "|------|-----|---------------------|---------|-------|",
    ]
    for rank, record in enumerate(records, start=1):
        aff = record.get("best_affinity")
        aff_str = f"{aff:.3f}" if aff is not None else "N/A"
        rmsd = record.get("best_mode_rmsd")
        rmsd_str = f"{rmsd:.3f}" if rmsd is not None else "N/A"
        lines.append(
            f"| {rank} | {record['cid']} | {aff_str} | {rmsd_str} | {record.get('num_poses', 0)} |"
        )
    if records and records[0].get("best_affinity") is not None:
        winner = records[0]
        lines.extend(
            [
                "",
                f"**Best candidate:** CID {winner['cid']} "
                f"({winner['best_affinity']:.3f} kcal/mol)",
            ]
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  [OK] Report: {out_path}")


def build_manifest(gene_entries: list[dict], exp: str, limit: int) -> dict:
    return {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "experiment_suffix": exp,
        "docking_limit": limit,
        "genes": gene_entries,
    }


def analyze_genes(genes: list[str], exp: str, limit: int) -> None:
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    manifest_entries: list[dict] = []

    for gene in genes:
        gene = gene.upper()
        records = load_summary_records(gene, exp, limit)
        if not records:
            print(f"  [UYARI] No results for {gene} exp_{exp}")
            continue

        rows = build_leaderboard_rows(gene, exp, records)
        all_rows.extend(rows)

        report_path = REPORTS_DIR / f"docking_report_{gene.lower()}_exp{exp}.md"
        write_gene_report(gene, exp, records, report_path)

        config = load_gene_config(gene)
        cids = load_docking_cids(gene, exp, limit=limit)
        manifest_entries.append(
            {
                "gene": gene,
                "exp_suffix": exp,
                "docking_limit": limit,
                "cids": cids,
                "docked_count": len(records),
                "pdb_id": config["pdb_id"],
                "summary_json": str(RunPaths.for_gene(gene, exp).summary_json),
                "best_cid": records[0]["cid"] if records else None,
                "best_affinity": records[0].get("best_affinity") if records else None,
            }
        )

    if all_rows:
        df = pd.DataFrame(all_rows)
        df = df.sort_values("best_affinity_kcal_mol", na_position="last")
        leaderboard_path = REPORTS_DIR / f"docking_leaderboard_exp{exp}.csv"
        df.to_csv(leaderboard_path, index=False)
        print(f"  [OK] Leaderboard: {leaderboard_path} ({len(df)} rows)")

    manifest = build_manifest(manifest_entries, exp, limit)
    manifest_path = REPORTS_DIR / f"docking_manifest_exp{exp}.json"
    with open(manifest_path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"  [OK] Manifest: {manifest_path}")


def parse_legacy_logs(legacy_logs: Path, gene: str, cids: list[int]) -> list[dict]:
    """Parse legacy log files using result_parser regex."""
    from docking.core.result_parser import parse_vina_log

    gene = gene.upper()
    summary: list[dict] = []
    for cid in cids:
        label = f"{gene}_CID_{cid}"
        log_path = legacy_logs / f"{label}_cmd.log"
        poses = parse_vina_log(log_path)
        if not poses:
            summary.append(
                {
                    "protein": gene,
                    "cid": cid,
                    "best_affinity": None,
                    "best_mode_rmsd": None,
                    "num_poses": 0,
                    "all_poses": [],
                }
            )
            continue
        best = poses[0]
        summary.append(
            {
                "protein": gene,
                "cid": cid,
                "best_affinity": best["affinity"],
                "best_mode_rmsd": best["rmsd_ub"],
                "num_poses": len(poses),
                "all_poses": poses,
            }
        )
    summary.sort(
        key=lambda row: row["best_affinity"] if row["best_affinity"] is not None else 0
    )
    return summary


def compare_with_legacy_json(parsed: list[dict], legacy_json: Path) -> bool:
    with open(legacy_json, encoding="utf-8") as fh:
        legacy = json.load(fh)
    legacy_map = {row["cid"]: row["best_affinity"] for row in legacy}
    ok = True
    for row in parsed:
        cid = row["cid"]
        expected = legacy_map.get(cid)
        actual = row["best_affinity"]
        if expected is None and actual is None:
            continue
        if expected is None or actual is None or abs(expected - actual) > 1e-6:
            print(f"  [FAIL] CID {cid}: expected={expected}, parsed={actual}")
            ok = False
        else:
            print(f"  [PASS] CID {cid}: {actual} kcal/mol")
    return ok


def main() -> None:
    args = parse_args()

    if args.legacy_logs:
        gene = (args.legacy_gene or args.genes[0]).upper()
        cids = load_docking_cids(gene, args.exp, limit=args.limit)
        parsed = parse_legacy_logs(args.legacy_logs, gene, cids)
        legacy_dir = PUBCHEM_ROOT / "docking" / f"{gene.lower()}_docking" / "results"
        legacy_json = legacy_dir / f"docking_summary_{gene.lower()}_exp{args.exp}.json"
        if not legacy_json.exists():
            fallback = legacy_dir / f"docking_summary_{gene.lower()}_exp06.json"
            if fallback.exists():
                legacy_json = fallback
        if legacy_json.exists():
            print(f"Regression check: {gene} exp_{args.exp} limit={args.limit} vs {legacy_json.name}")
            ok = compare_with_legacy_json(parsed, legacy_json)
            sys.exit(0 if ok else 1)
        out = REPORTS_DIR / f"legacy_parsed_{gene.lower()}_exp{args.exp}.json"
        save_summary_json(parsed, out)
        print(f"  [OK] Legacy parse output: {out}")
        return

    print(f"Analyzing genes={args.genes} exp_{args.exp} limit={args.limit}")
    analyze_genes(args.genes, args.exp, args.limit)


if __name__ == "__main__":
    main()
