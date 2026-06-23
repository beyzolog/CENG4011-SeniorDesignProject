# updated_docking — Modular Molecular Docking Pipeline

Stage 8 of the PubChem screening pipeline. Runs AutoDock Vina on the top-N
scoring candidates from the step-07 final-1000 pool (default N=100).

## Layout

```
docking/
├── config/genes.yaml       # MYC + CTNNB1 Vina box parameters
├── core/                   # Shared pipeline modules
├── shared/proteins/        # Receptor PDB/PDBQT (shared across runs)
├── shared/ligands/         # Ligand SDF/PDBQT (shared across runs)
├── bin/vina                # AutoDock Vina binary (not in git)
└── runs/{GENE}/exp_{NN}/   # Per-run configs, logs, results (gitignored)
    ├── configs/
    ├── logs/
    ├── results/
    │   └── docking_summary.json
    └── ...
```

## Prerequisites

```bash
cd pubchem/updated_docking
pip install -r requirements.txt
```

Install AutoDock Vina 1.2.5:

```bash
mkdir -p bin && cd bin
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64 && mv vina_1.2.5_linux_x86_64 vina
```

Alternatively, reuse the legacy binary at `pubchem/docking/bin/vina` (auto-detected).

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `TARGET_GENE` | `CTNNB1` | `MYC`, `CTNNB1`, or `both` |
| `EXPERIMENT_SUFFIX` | `06` | Experiment folder suffix (e.g. `06`, `07`) |
| `DOCKING_LIMIT` | `100` | Top-N ligands from final-1000 pool to dock |
| `DOCKING_CPU` | `0` | Vina CPU limit (`0` = automatic) |

## Usage

From `pubchem/scripts/`:

```bash
# Dry-run: show top-100 CIDs without docking
TARGET_GENE=CTNNB1 EXPERIMENT_SUFFIX=06 python3 08_docking_orchestrator.py --dry-run

# Run CTNNB1 exp06 — top 100 ligands
TARGET_GENE=CTNNB1 EXPERIMENT_SUFFIX=06 python3 08_docking_orchestrator.py

# MYC exp07 (use exp07 for MYC final-1000 pool)
TARGET_GENE=MYC EXPERIMENT_SUFFIX=07 python3 08_docking_orchestrator.py --limit 100 --cpu 16

# Custom limit
DOCKING_LIMIT=50 TARGET_GENE=CTNNB1 python3 08_docking_orchestrator.py --dry-run

# Both genes in parallel (same --exp required; prefer separate runs per gene)
TARGET_GENE=both EXPERIMENT_SUFFIX=06 python3 08_docking_orchestrator.py --parallel
```

Validate parser regression against legacy top-5 outputs:

```bash
cd pubchem/updated_docking && python3 validate_seed.py
```

Cross-gene analysis:

```bash
python3 09_docking_analysis.py --exp 06 --limit 100
python3 09_docking_analysis.py --exp 07 --genes MYC --limit 100
```

## Upstream integration

CIDs are read from the step-07 final-1000 candidate pool (prediction_score order):

```
screening/{GENE}/predictions_exp_{NN}/final_1000_cids.tmp   # primary
screening/{GENE}/predictions_exp_{NN}/{gene}_exp{NN}_final_1000.csv  # fallback
```

First `DOCKING_LIMIT` entries are docked. Recommended experiment suffixes:

- **CTNNB1** → `EXPERIMENT_SUFFIX=06`
- **MYC** → `EXPERIMENT_SUFFIX=07`

## Legacy code

The previous `pubchem/docking/` tree is archived. See [docs/legacy_docking.md](../docs/legacy_docking.md).
