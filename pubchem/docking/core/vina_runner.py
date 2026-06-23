"""AutoDock Vina config generation and subprocess execution."""

from __future__ import annotations

import os
import shutil
import subprocess
import time
from pathlib import Path

from docking.core.paths import RunPaths


def find_vina(bin_dir: Path) -> str | None:
    local_vina = bin_dir / "vina"
    if local_vina.exists() and os.access(local_vina, os.X_OK):
        return str(local_vina)
    return shutil.which("vina") or shutil.which("autodock_vina")


def write_vina_config(
    paths: RunPaths,
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    config: dict,
    label: str,
    cpu: int = 0,
) -> Path:
    cfg_path = paths.configs_dir / f"config_{label}.txt"
    out_pdbqt = paths.results_dir / f"{label}_out.pdbqt"

    lines = [
        f"receptor = {receptor_pdbqt.resolve()}",
        f"ligand = {ligand_pdbqt.resolve()}",
        "",
        f"center_x = {config['center_x']}",
        f"center_y = {config['center_y']}",
        f"center_z = {config['center_z']}",
        "",
        f"size_x = {config['size_x']}",
        f"size_y = {config['size_y']}",
        f"size_z = {config['size_z']}",
        "",
        f"exhaustiveness = {config['exhaustiveness']}",
        f"num_modes = {config['num_modes']}",
        f"energy_range = {config['energy_range']}",
        "",
        f"out = {out_pdbqt.resolve()}",
    ]
    if cpu > 0:
        lines.append(f"cpu = {cpu}")

    cfg_path.write_text("\n".join(lines))
    print(f"    [OK] Config yazıldı: {cfg_path.name}")
    return cfg_path


def run_vina_process(cmd: list[str], log_file: Path, label: str) -> bool:
    print(f"    $ {' '.join(str(c) for c in cmd)}")
    with open(log_file, "w") as fh:
        result = subprocess.run(cmd, stdout=fh, stderr=fh)
    if result.returncode != 0:
        print(f"    [HATA] {label}: returncode={result.returncode}")
        try:
            print(log_file.read_text(errors="ignore")[-600:])
        except OSError:
            pass
        return False
    return True


def run_all_dockings(
    paths: RunPaths,
    protein_name: str,
    receptor_pdbqt: Path,
    config: dict,
    pdbqt_map: dict[int, Path],
    cpu: int = 0,
) -> list[dict]:
    vina_bin = find_vina(paths.bin_dir)
    if vina_bin is None:
        print(
            "  [HATA] Vina binary bulunamadı.\n"
            "  Çözüm:\n"
            "    cd pubchem/docking/bin\n"
            "    wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/"
            "download/v1.2.5/vina_1.2.5_linux_x86_64\n"
            "    chmod +x vina_1.2.5_linux_x86_64 && mv vina_1.2.5_linux_x86_64 vina"
        )
        return []

    print(
        f"  [Vina] binary: {vina_bin}"
        + (f"  |  CPU sınırı: {cpu}" if cpu > 0 else "  |  CPU: otomatik")
    )

    results: list[dict] = []
    for cid, lig_pdbqt in pdbqt_map.items():
        label = f"{protein_name}_CID_{cid}"
        out_pdbqt = paths.results_dir / f"{label}_out.pdbqt"
        cmd_log = paths.logs_dir / f"{label}_cmd.log"

        if out_pdbqt.exists() and cmd_log.exists() and out_pdbqt.stat().st_size > 1024:
            print(f"  [ATLA] Docking mevcut: {label}")
            results.append({"protein": protein_name, "cid": cid, "success": True})
            continue

        cfg = write_vina_config(paths, receptor_pdbqt, lig_pdbqt, config, label, cpu)
        cmd = [vina_bin, "--config", str(cfg)]

        print(f"\n  [VINA] {label} başlıyor...")
        t0 = time.time()
        ok = run_vina_process(cmd, cmd_log, label)
        elapsed = time.time() - t0
        status = "TAMAMLANDI" if ok else "BAŞARISIZ"
        print(f"  [{status}] {label} — {elapsed:.1f}s")
        results.append({"protein": protein_name, "cid": cid, "success": ok})

    return results
