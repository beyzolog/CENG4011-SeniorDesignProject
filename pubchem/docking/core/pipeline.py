"""End-to-end docking pipeline for a single gene/experiment."""

from __future__ import annotations

from docking.core.gene_config import load_gene_config
from docking.core.ligand_prep import prepare_ligands
from docking.core.paths import RunPaths, ensure_run_dirs
from docking.core.protein_prep import prepare_protein
from docking.core.result_parser import (
    DockingSummary,
    parse_all_results,
    print_summary_table,
)
from docking.core.vina_runner import run_all_dockings


def run_gene_docking(
    gene: str,
    cids: list[int],
    paths: RunPaths,
    cpu: int = 0,
) -> DockingSummary:
    gene = gene.upper()
    config = load_gene_config(gene)

    print("=" * 70)
    print(f"  Docking Pipeline  ·  {gene} exp_{paths.exp_suffix}")
    print(f"  Protein   : {gene} ({config['pdb_id']})")
    print(f"  CIDs      : {cids}")
    print(f"  Box       : {config['size_x']}³ Å")
    print(f"  Exhaustiveness: {config['exhaustiveness']}")
    print(f"  CPU       : {'otomatik' if cpu == 0 else cpu}")
    print(f"  Çıktı     : {paths.run_dir}")
    print("=" * 70)

    ensure_run_dirs(paths)

    print("\n[ADIM 1] Protein hazırlanıyor...")
    receptor = prepare_protein(config["pdb_id"], paths.proteins_dir)

    print("\n[ADIM 2] Ligandlar hazırlanıyor...")
    pdbqt_map = prepare_ligands(cids, paths.ligands_dir)

    print("\n[ADIM 3–4] Grid config + Vina Docking...")
    run_all_dockings(paths, gene, receptor, config, pdbqt_map, cpu=cpu)

    print("\n[ADIM 5] Sonuçlar ayrıştırılıyor...")
    records = parse_all_results(paths, gene, cids)
    print_summary_table(records, gene, paths.exp_suffix)

    summary = DockingSummary(
        gene=gene,
        exp_suffix=paths.exp_suffix,
        cids=cids,
        records=records,
    )
    summary.write_json(paths.summary_json)

    print(f"\n{'=' * 70}")
    print(f"  [BİTTİ] {gene} exp_{paths.exp_suffix} docking tamamlandı.")
    print(f"  Sonuçlar : {paths.results_dir.resolve()}")
    print(f"  Loglar   : {paths.logs_dir.resolve()}")
    print(f"  JSON     : {paths.summary_json.resolve()}")
    print(f"{'=' * 70}")

    return summary
