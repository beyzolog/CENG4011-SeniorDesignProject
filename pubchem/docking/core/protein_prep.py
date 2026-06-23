"""Protein preparation: RCSB download, cleaning, receptor PDBQT."""

from __future__ import annotations

from pathlib import Path

import requests
from Bio.PDB import PDBIO, PDBParser, Select

_AD4_ELEMENT_MAP: dict[str, str] = {
    "C": "C",
    "N": "N",
    "O": "OA",
    "S": "SA",
    "H": "H",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I": "I",
    "FE": "Fe",
    "ZN": "Zn",
    "MG": "Mg",
    "CA": "Ca",
    "MN": "Mn",
    "CU": "Cu",
}

_AROMATIC_RING_ATOMS: dict[str, set] = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
}


def _get_ad4_type(element: str, atom_name: str, res_name: str) -> str:
    el = element.strip().upper()
    nm = atom_name.strip().upper()
    rn = res_name.strip().upper()
    if el == "C" and nm in _AROMATIC_RING_ATOMS.get(rn, set()):
        return "A"
    return _AD4_ELEMENT_MAP.get(el, "C")


def _pdb_to_pdbqt_receptor(clean_pdb: Path, pdbqt_out: Path) -> None:
    parser_obj = PDBParser(QUIET=True)
    struct = parser_obj.get_structure("rec", clean_pdb)
    lines: list[str] = []
    serial = 0
    for model in struct:
        for chain in model:
            cid = chain.get_id()
            for res in chain:
                het, resseq, _ = res.get_id()
                if het.strip():
                    continue
                rname = res.get_resname().strip()
                for atom in res.get_atoms():
                    if atom.is_disordered() and atom.get_altloc() not in ("A", " "):
                        continue
                    serial += 1
                    aname = atom.get_name()
                    element = (atom.element or aname[0]).strip().upper()
                    x, y, z = atom.get_vector()
                    occ = atom.get_occupancy() or 1.0
                    bfac = atom.get_bfactor() or 0.0
                    ad4 = _get_ad4_type(element, aname, rname)
                    padded = f" {aname:<3s}" if len(aname) < 4 else f"{aname:<4s}"
                    line = (
                        f"ATOM  {serial:5d} {padded} {rname:<3s} "
                        f"{cid}{resseq:4d}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}"
                        f"{occ:6.2f}{bfac:6.2f}    "
                        f"{'0.000':>6s} {ad4:<2s}"
                    )
                    lines.append(line)
    pdbqt_out.write_text("\n".join(lines) + "\n")


class CleanProtein(Select):
    def accept_residue(self, residue):
        return not residue.get_id()[0].strip()

    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() in ("A", " ")


def download_pdb(pdb_id: str, out_dir: Path) -> Path:
    raw = out_dir / f"{pdb_id}_raw.pdb"
    if raw.exists():
        print(f"  [ATLA] {pdb_id} zaten mevcut: {raw.name}")
        return raw
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"  [İNDİR] {url}")
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    raw.write_bytes(response.content)
    return raw


def prepare_protein(pdb_id: str, out_dir: Path) -> Path:
    """Raw PDB -> clean PDB -> receptor PDBQT."""
    raw_pdb = download_pdb(pdb_id, out_dir)
    clean_pdb = out_dir / f"{pdb_id}_clean.pdb"
    pdbqt_out = out_dir / f"{pdb_id}_receptor.pdbqt"

    if not clean_pdb.exists():
        parser_obj = PDBParser(QUIET=True)
        struct = parser_obj.get_structure(pdb_id, raw_pdb)
        io = PDBIO()
        io.set_structure(struct)
        io.save(str(clean_pdb), CleanProtein())
        print(f"  [OK] Temizlendi → {clean_pdb.name}")
    else:
        print(f"  [ATLA] Temiz PDB mevcut: {clean_pdb.name}")

    if not pdbqt_out.exists():
        _pdb_to_pdbqt_receptor(clean_pdb, pdbqt_out)
        print(f"  [OK] PDBQT oluşturuldu → {pdbqt_out.name}")
    else:
        print(f"  [ATLA] PDBQT mevcut: {pdbqt_out.name}")

    return pdbqt_out
