"""Central path resolution for docking runs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent.parent
PUBCHEM_ROOT = PACKAGE_ROOT.parent
LEGACY_DOCKING_ROOT = PUBCHEM_ROOT / "docking"
SCREENING_ROOT = PUBCHEM_ROOT / "screening"


@dataclass(frozen=True)
class RunPaths:
    gene: str
    exp_suffix: str
    run_dir: Path
    configs_dir: Path
    logs_dir: Path
    results_dir: Path
    summary_json: Path
    proteins_dir: Path
    ligands_dir: Path
    bin_dir: Path

    @classmethod
    def for_gene(cls, gene: str, exp_suffix: str) -> RunPaths:
        from docking.core.champion_loader import normalize_exp_suffix

        gene = gene.upper()
        exp_suffix = normalize_exp_suffix(exp_suffix)
        run_dir = PACKAGE_ROOT / "runs" / gene / f"exp_{exp_suffix}"
        legacy_bin = LEGACY_DOCKING_ROOT / "bin"
        bin_dir = PACKAGE_ROOT / "bin"
        if not bin_dir.exists() and legacy_bin.exists():
            bin_dir = legacy_bin
        return cls(
            gene=gene,
            exp_suffix=exp_suffix,
            run_dir=run_dir,
            configs_dir=run_dir / "configs",
            logs_dir=run_dir / "logs",
            results_dir=run_dir / "results",
            summary_json=run_dir / "results" / "docking_summary.json",
            proteins_dir=PACKAGE_ROOT / "shared" / "proteins",
            ligands_dir=PACKAGE_ROOT / "shared" / "ligands",
            bin_dir=bin_dir,
        )

    def clusters_csv(self) -> Path:
        exp_dir = SCREENING_ROOT / self.gene / f"predictions_exp_{self.exp_suffix}"
        return exp_dir / f"{self.gene.lower()}_{self.exp_suffix}_clusters.csv"


def ensure_run_dirs(paths: RunPaths) -> None:
    for directory in (
        paths.run_dir,
        paths.configs_dir,
        paths.logs_dir,
        paths.results_dir,
        paths.proteins_dir,
        paths.ligands_dir,
        paths.bin_dir,
    ):
        directory.mkdir(parents=True, exist_ok=True)
