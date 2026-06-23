"""Load docking candidate CIDs from upstream screening output (final 1000 pool)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from docking.core.paths import SCREENING_ROOT

DEFAULT_DOCKING_LIMIT = 100


def normalize_exp_suffix(exp_suffix: str) -> str:
    """Normalize '06', 'exp06', or '_06' to '06'."""
    suffix = str(exp_suffix).lstrip("_")
    if suffix.lower().startswith("exp"):
        suffix = suffix[3:]
    return suffix.lstrip("_")


def experiment_dir(gene: str, exp_suffix: str) -> Path:
    gene = gene.upper()
    suffix = normalize_exp_suffix(exp_suffix)
    return SCREENING_ROOT / gene / f"predictions_exp_{suffix}"


def final_1000_cids_path(gene: str, exp_suffix: str) -> Path:
    return experiment_dir(gene, exp_suffix) / "final_1000_cids.tmp"


def final_1000_csv_path(gene: str, exp_suffix: str) -> Path:
    gene = gene.upper()
    suffix = normalize_exp_suffix(exp_suffix)
    return experiment_dir(gene, exp_suffix) / f"{gene.lower()}_exp{suffix}_final_1000.csv"


def clusters_csv_path(gene: str, exp_suffix: str) -> Path:
    gene = gene.upper()
    suffix = normalize_exp_suffix(exp_suffix)
    return experiment_dir(gene, exp_suffix) / f"{gene.lower()}_exp{suffix}_clusters.csv"


def _load_cids_from_tmp(tmp_path: Path, limit: int) -> list[int]:
    lines = [
        line.strip()
        for line in tmp_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(lines) < limit:
        raise ValueError(
            f"Expected at least {limit} CIDs in {tmp_path}, found {len(lines)}"
        )
    return [int(line) for line in lines[:limit]]


def _load_cids_from_csv(csv_path: Path, limit: int) -> list[int]:
    df = pd.read_csv(csv_path)
    if "prediction_score" in df.columns:
        df = df.sort_values("prediction_score", ascending=False)
    df = df.head(limit)
    if len(df) < limit:
        raise ValueError(
            f"Expected at least {limit} compounds in {csv_path}, found {len(df)}"
        )
    return df["CID"].astype(int).tolist()


def load_docking_cids(
    gene: str,
    exp_suffix: str,
    limit: int = DEFAULT_DOCKING_LIMIT,
) -> list[int]:
    """Return top-N CIDs by prediction_score from the final-1000 candidate pool."""
    tmp_path = final_1000_cids_path(gene, exp_suffix)
    if tmp_path.exists():
        return _load_cids_from_tmp(tmp_path, limit)

    csv_path = final_1000_csv_path(gene, exp_suffix)
    if csv_path.exists():
        return _load_cids_from_csv(csv_path, limit)

    raise FileNotFoundError(
        f"No docking candidate list found for {gene} exp_{normalize_exp_suffix(exp_suffix)}.\n"
        f"  Expected: {tmp_path}\n"
        f"  Or:       {csv_path}\n"
        f"Run scripts/07_toxicity_and_grouping.py with TARGET_GENE={gene} first."
    )


def resolve_docking_source(gene: str, exp_suffix: str) -> Path:
    """Return the path actually used for CID loading (for dry-run logging)."""
    tmp_path = final_1000_cids_path(gene, exp_suffix)
    if tmp_path.exists():
        return tmp_path
    csv_path = final_1000_csv_path(gene, exp_suffix)
    if csv_path.exists():
        return csv_path
    return tmp_path


def load_prediction_scores(
    gene: str, exp_suffix: str, cids: list[int]
) -> dict[int, float | None]:
    """Map CID -> prediction_score from the final-1000 CSV."""
    csv_path = final_1000_csv_path(gene, exp_suffix)
    if not csv_path.exists():
        return {cid: None for cid in cids}
    df = pd.read_csv(csv_path)
    score_map = dict(zip(df["CID"].astype(int), df["prediction_score"]))
    return {cid: score_map.get(cid) for cid in cids}


def load_top5_cids(gene: str, exp_suffix: str, n: int = 5) -> list[int]:
    """Backward-compatible alias: top-N from final-1000 pool (default N=5)."""
    return load_docking_cids(gene, exp_suffix, limit=n)


def load_champion_scores(
    gene: str, exp_suffix: str, cids: list[int]
) -> dict[int, float | None]:
    """Backward-compatible alias for load_prediction_scores."""
    return load_prediction_scores(gene, exp_suffix, cids)
