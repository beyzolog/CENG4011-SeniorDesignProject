"""Ligand preparation: PubChem SDF download and Meeko PDBQT conversion."""

from __future__ import annotations

import time
from pathlib import Path

import requests
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops


def fetch_sdf(cid: int, out_dir: Path) -> Path:
    sdf_path = out_dir / f"CID_{cid}.sdf"
    if sdf_path.exists():
        print(f"    [ATLA] SDF mevcut: {sdf_path.name}")
        return sdf_path

    response = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        f"?record_type=3d",
        timeout=30,
    )
    if response.status_code == 200 and len(response.content) > 100:
        sdf_path.write_bytes(response.content)
        print(f"    [OK] 3D SDF indirildi: {sdf_path.name}")
        return sdf_path

    print(f"    [UYARI] CID {cid}: kayıtlı 3D yok → 2D→3D (RDKit MMFF94)")
    response_2d = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
        f"?record_type=2d",
        timeout=30,
    )
    response_2d.raise_for_status()
    tmp = out_dir / f"CID_{cid}_2d.sdf"
    tmp.write_bytes(response_2d.content)

    mol = Chem.MolFromMolFile(str(tmp), removeHs=False)
    if mol is None:
        raise ValueError(f"CID {cid}: RDKit mol nesnesi oluşturulamadı.")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    writer = Chem.SDWriter(str(sdf_path))
    writer.write(mol)
    writer.close()
    tmp.unlink()
    print(f"    [OK] 3D konformer üretildi: {sdf_path.name}")
    return sdf_path


def sdf_to_pdbqt(sdf_path: Path, out_dir: Path) -> Path:
    pdbqt_path = out_dir / (sdf_path.stem + ".pdbqt")
    if pdbqt_path.exists():
        print(f"    [ATLA] PDBQT mevcut: {pdbqt_path.name}")
        return pdbqt_path

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next(iter(supplier))
    if mol is None:
        raise ValueError(f"{sdf_path.name} okunamadı.")

    if len(rdmolops.GetMolFrags(mol)) > 1:
        frag_count = len(rdmolops.GetMolFrags(mol))
        print(f"    [UYARI] Çoklu fragment ({frag_count}), en büyüğü seçiliyor...")
        frags = rdmolops.GetMolFrags(mol, asMols=True)
        mol = max(frags, key=lambda x: x.GetNumAtoms())

    if mol.GetNumAtoms() == mol.GetNumHeavyAtoms():
        mol = Chem.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)

    with open(pdbqt_path, "w") as fh:
        for setup in mol_setups:
            pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setup)
            if not is_ok:
                raise RuntimeError(f"Meeko PDBQT hatası ({sdf_path.name}): {err}")
            fh.write(pdbqt_str)

    print(f"    [OK] PDBQT oluşturuldu: {pdbqt_path.name}")
    return pdbqt_path


def prepare_ligands(cids: list[int], out_dir: Path) -> dict[int, Path]:
    pdbqt_map: dict[int, Path] = {}
    for cid in cids:
        print(f"\n  Ligand CID {cid} hazırlanıyor...")
        try:
            sdf = fetch_sdf(cid, out_dir)
            time.sleep(0.3)
            pdbqt = sdf_to_pdbqt(sdf, out_dir)
            pdbqt_map[cid] = pdbqt
        except Exception as exc:
            print(f"    [HATA] CID {cid}: {exc}")
    return pdbqt_map
