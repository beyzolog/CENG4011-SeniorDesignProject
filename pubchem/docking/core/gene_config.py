"""Load gene docking configuration from YAML."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from docking.core.paths import PACKAGE_ROOT

_GENES_CONFIG = PACKAGE_ROOT / "config" / "genes.yaml"


def load_gene_config(gene: str) -> dict[str, Any]:
    gene = gene.upper()
    with open(_GENES_CONFIG, encoding="utf-8") as fh:
        all_genes = yaml.safe_load(fh)
    if gene not in all_genes:
        raise KeyError(f"Unknown gene '{gene}'. Available: {list(all_genes)}")
    raw = all_genes[gene]
    center = raw["center"]
    box = raw["box_size"]
    return {
        "pdb_id": raw["pdb_id"],
        "chain": raw.get("chain", "A"),
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": box,
        "size_y": box,
        "size_z": box,
        "exhaustiveness": raw["exhaustiveness"],
        "num_modes": raw["num_modes"],
        "energy_range": raw["energy_range"],
    }
