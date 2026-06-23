"""Parse Vina logs and produce docking summary JSON."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

from docking.core.paths import RunPaths


def parse_vina_log(log_path: Path) -> list[dict]:
    if not log_path.exists():
        return []
    poses: list[dict] = []
    pattern = re.compile(
        r"^\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)", re.MULTILINE
    )
    for match in pattern.finditer(log_path.read_text(errors="ignore")):
        poses.append(
            {
                "mode": int(match.group(1)),
                "affinity": float(match.group(2)),
                "rmsd_lb": float(match.group(3)),
                "rmsd_ub": float(match.group(4)),
            }
        )
    return poses


def parse_all_results(
    paths: RunPaths, protein_name: str, cids: list[int]
) -> list[dict]:
    summary: list[dict] = []
    for cid in cids:
        label = f"{protein_name}_CID_{cid}"
        poses = parse_vina_log(paths.logs_dir / f"{label}_cmd.log")
        if not poses:
            summary.append(
                {
                    "protein": protein_name,
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
                "protein": protein_name,
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


def print_summary_table(summary: list[dict], protein_name: str, exp_suffix: str):
    sep = "─" * 70
    display_limit = 20
    print(f"\n{'═' * 70}")
    print(f"  {protein_name} (exp_{exp_suffix}) — Docking Sonuçları ({len(summary)} aday)")
    print(f"{'═' * 70}")
    print(f"  {'CID':<15} {'Affinity (kcal/mol)':<22} {'RMSD ub':<12} {'Poz'}")
    print(sep)
    for row in summary[:display_limit]:
        aff = (
            f"{row['best_affinity']:.2f}"
            if row["best_affinity"] is not None
            else "N/A"
        )
        rmsd = (
            f"{row['best_mode_rmsd']:.3f}"
            if row["best_mode_rmsd"] is not None
            else "N/A"
        )
        print(f"  {row['cid']:<15} {aff:<22} {rmsd:<12} {row['num_poses']}")
    if len(summary) > display_limit:
        print(f"  ... +{len(summary) - display_limit} daha (tam liste: docking_summary.json)")
    print(sep)
    if summary and summary[0]["best_affinity"] is not None:
        winner = summary[0]
        print(
            f"\n  ★  En iyi aday: CID {winner['cid']}  "
            f"({winner['best_affinity']:.2f} kcal/mol)"
        )
    print()


def save_summary_json(summary: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  [OK] JSON kaydedildi: {out_path}")


@dataclass
class DockingSummary:
    gene: str
    exp_suffix: str
    cids: list[int]
    records: list[dict] = field(default_factory=list)

    def write_json(self, out_path: Path) -> None:
        save_summary_json(self.records, out_path)

    def to_manifest_entry(self) -> dict:
        return {
            "gene": self.gene,
            "exp_suffix": self.exp_suffix,
            "cids": self.cids,
            "best_cid": self.records[0]["cid"] if self.records else None,
            "best_affinity": self.records[0]["best_affinity"] if self.records else None,
        }
