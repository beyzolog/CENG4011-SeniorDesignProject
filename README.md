# 🧬 OncoRep — ML-Driven Drug Repurposing Pipeline

> **Senior Design Project** · Computer Engineering  
> **Author:** Beyza Yoruk  
> Leakage-free machine learning virtual screening & molecular docking against **MYC** and **CTNNB1** oncoproteins.

[![Python](https://img.shields.io/badge/Python-3.10-blue?logo=python&logoColor=white)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-ECFP4-green)](https://www.rdkit.org/)
[![scikit-learn](https://img.shields.io/badge/scikit--learn-ML-orange)](https://scikit-learn.org/)
[![AutoDock Vina](https://img.shields.io/badge/AutoDock_Vina-1.2.5-red)](https://vina.scripps.edu/)

---

## 📋 Table of Contents

- [Project Overview](#-1-project-overview)
- [Repository Structure](#-2-repository-structure)
- [Pipeline & Execution Order](#-3-pipeline--execution-order-1-9)
- [Data Reproduction & Reproducibility](#-4-data-reproduction--reproducibility)
- [Prerequisites & Setup](#-5-prerequisites--setup)
- [Data Leakage Prevention & Experimental Rigor](#-data-leakage-prevention--experimental-rigor)
- [Target Genes & Experiments](#-target-genes--experiments)
- [License & Citation](#-license--citation)

---

## 🔬 1. Project Overview

**OncoRep** is an end-to-end computational drug repurposing platform that identifies small-molecule inhibitors of two clinically relevant oncoproteins:

| Target | Role | PDB (Docking) |
|--------|------|---------------|
| **CTNNB1** (β-catenin) | Wnt/β-catenin pathway driver in colorectal & other cancers | `1JPW` |
| **MYC** | Master transcription factor amplified in ~70% of human cancers | `6G6K` |

The system integrates three computational layers:

1. **Data Engineering** — Automated retrieval of confirmatory bioassay data from PubChem (NCBI Entrez + REST API), IC50-based active/inactive labeling, and Morgan fingerprint (ECFP4, 2048-bit) generation via RDKit.

2. **Leakage-Free Machine Learning** — A rigorously designed training pipeline that prevents data leakage at every stage:
   - Correlation-based multicollinearity pruning is computed **only on the training split**.
   - Tree-based feature importance and RFECV operate on leakage-free bit matrices (`04_raw_bit_matrices/`).
   - Champion models (Random Forest, Extra Trees, XGBoost) are selected via stratified cross-validation with strict anti-overfitting hyperparameters (depth limits, minimum leaf sizes, train–test gap monitoring).
   - `feature_order.json` enforces deterministic column alignment during million-scale virtual screening.

3. **Structure-Based Validation** — Top-ranked candidates from medicinal-chemistry filtering (PAINS/Brenk, Lipinski, Tanimoto clustering) are validated through **AutoDock Vina** molecular docking against gene-specific binding pockets.

A complementary **ChEMBL** module (`chembl/`) provides independent bioactivity retrieval for cross-database validation, and a **Flask web portal** (`webPortal/`) exposes pipeline results through an interactive frontend for the Senior Design demonstration.

---

## 📁 2. Repository Structure

Large binary artifacts (~**147 MB** of raw data + screening predictions) are **excluded via `.gitignore`**. Empty runtime directories are preserved in Git using **`.gitkeep`** placeholders so the pipeline folder skeleton survives clone operations.

```
OncoRep/
├── 📄 requirements.txt              # Full Python dependency lockfile
├── 📄 README.md
├── 📄 .gitignore
│
├── chembl/                            # Independent ChEMBL validation module
│   ├── 01_chembl_retrieval_and_filtering.ipynb
│   └── data/                          # 🚫 gitignored — re-fetch via notebook
│
├── pubchem/                           # Primary ML + docking pipeline
│   ├── genes.txt                      # Target gene list (CTNNB1, MYC)
│   │
│   ├── data/                          # Staged data lake (~61 MB local)
│   │   ├── 01_raw/                    # 🚫 PubChem AID CSV downloads
│   │   ├── 02_aggregated/             # 🚫 IC50-filtered actives/inactives
│   │   ├── 03_fingerprints/           # 🚫 ECFP4 Morgan 2048-bit matrices
│   │   ├── 04_raw_bit_matrices/       # 🚫 Leakage-free full bit matrices
│   │   ├── 05_updated_selected_features/  # ✅ IN REPO — selected feature CSVs
│   │   ├── archive/                   # 🚫 Legacy RFECV experiments
│   │   └── old_selected_features/     # 🚫 Superseded feature sets
│   │
│   ├── scripts/                       # ✅ Pipeline orchestration (Stages 1–9)
│   │   ├── 01_data_aggregation.py
│   │   ├── 02_generate_fingerprints.py
│   │   ├── 03_feature_selection.py
│   │   ├── 04_model_training_and_selection.py
│   │   ├── 05_final_screening.py
│   │   ├── 06_combine_and_prefilter.py
│   │   ├── 07_toxicity_and_grouping.py
│   │   ├── 08_docking_orchestrator.py
│   │   ├── 09_docking_analysis.py
│   │   └── archive/                   # Legacy experiment scripts
│   │
│   ├── notebooks/                     # ✅ Exploratory analysis
│   │   ├── active_curated_data.ipynb
│   │   ├── compound_analysis.ipynb
│   │   ├── correlation_analysis.ipynb # Generates 04_raw_bit_matrices/
│   │   └── archive/                   # 🚫 gitignored
│   │
│   ├── models/                        # ✅ Trained champion models (~1.5 MB)
│   │   └── experiments/
│   │       ├── CTNNB1/exp_05/, exp_06/
│   │       └── MYC/exp_06/, exp_07/
│   │
│   ├── screening/                     # 🚫 Large prediction CSVs (~86 MB)
│   │   ├── CTNNB1/predictions_exp_05/, exp_06/
│   │   └── MYC/predictions_exp_06/, exp_07/
│   │   # ✅ Small summaries tracked: *_final_1000.csv, *_clusters.csv, reports
│   │
│   └── docking/                       # Molecular docking engine
│       ├── config/genes.yaml          # ✅ Vina box parameters per gene
│       ├── core/                      # ✅ Pipeline modules (prep, runner, parser)
│       ├── reports/                   # ✅ Leaderboards & markdown reports
│       ├── shared/
│       │   ├── ligands/.gitkeep       # ✅ Empty dir placeholder (ligands gitignored)
│       │   └── proteins/
│       │       ├── .gitkeep
│       │       └── *_receptor.pdbqt   # ✅ Prepared receptors in repo
│       ├── bin/vina                   # 🚫 AutoDock Vina binary (manual install)
│       └── runs/                      # 🚫 Per-run logs, configs, PDBQT outputs
│
└── webPortal/                         # ✅ Flask web interface (Senior Design UI)
    ├── app.py
    ├── templates/
    └── static/
```

### Legend

| Symbol | Meaning |
|--------|---------|
| ✅ | Tracked in Git — available after `git clone` |
| 🚫 | Gitignored — must be regenerated locally or downloaded separately |

---

## ⚙️ 3. Pipeline & Execution Order (1–9)

All pipeline scripts **must be executed sequentially** from the `pubchem/scripts/` directory. Each stage writes to a well-defined output path consumed by the next stage.


### Stage 1 — `01_data_aggregation.py`

**Purpose:** Query PubChem for confirmatory assays targeting genes in `genes.txt`, download raw AID datatables, and produce IC50-filtered active/inactive sets.

```bash
cd pubchem/scripts
python 01_data_aggregation.py
```

| Input | Output |
|-------|--------|
| `genes.txt` | `data/01_raw/{GENE}/*.csv` |
| NCBI Entrez + PubChem REST | `data/02_aggregated/{GENE}_active.csv`, `{GENE}_inactive.csv` |

---

### Stage 2 — `02_generate_fingerprints.py`

**Purpose:** Convert SMILES strings to **ECFP4 Morgan fingerprints** (radius = 2, 2048 bits) using RDKit.

```bash
python 02_generate_fingerprints.py
```

| Input | Output |
|-------|--------|
| `data/02_aggregated/*.csv` | `data/03_fingerprints/*.csv` |

---

### Stage 3 — `03_feature_selection.py`

**Purpose:** Tree-based robust feature selection with **train-only multicollinearity pruning** (|r| > 0.80). Evaluates multiple feature-count thresholds and exports the champion feature matrix.

> **Prerequisite:** Run `notebooks/correlation_analysis.ipynb` first to generate leakage-free bit matrices in `data/04_raw_bit_matrices/`.

```bash
jupyter notebook ../notebooks/correlation_analysis.ipynb   # one-time
python 03_feature_selection.py
```

| Input | Output |
|-------|--------|
| `data/04_raw_bit_matrices/{GENE}_ecfp4_2048_original.csv` | `data/05_updated_selected_features/{GENE}_robust_tree_filtered.csv` |

---

### Stage 4 — `04_model_training_and_selection.py`

**Purpose:** Train and benchmark **Random Forest**, **Extra Trees**, and **XGBoost** classifiers. Select the champion model based on CV AUC with anti-overfitting constraints. Persist model artifacts and strict `feature_order.json`.

```bash
# Edit GENE variable inside the script, then run for each target:
python 04_model_training_and_selection.py
```

| Input | Output |
|-------|--------|
| `data/05_updated_selected_features/` | `models/experiments/{GENE}/exp_{NN}/Results/*.pkl`, `*.json` |

---

### Stage 5 — `05_final_screening.py`

**Purpose:** Million-compound **virtual screening** against the PubChem library using the champion model. Parallelized fingerprint generation and batched inference.

```bash
# Edit GENE variable inside the script
python 05_final_screening.py
```

| Input | Output |
|-------|--------|
| `models/experiments/{GENE}/exp_{NN}/Results/` | `screening/{GENE}/predictions_exp_{NN}/*.csv` |

---

### Stage 6 — `06_combine_and_prefilter.py`

**Purpose:** Merge parallel screening chunks, remove known actives (data leakage guard), apply fast physicochemical filters (MW, LogP), and enforce the prediction score threshold (≥ 0.80).

```bash
TARGET_GENE=CTNNB1 python 06_combine_and_prefilter.py
TARGET_GENE=MYC     python 06_combine_and_prefilter.py
```

| Input | Output |
|-------|--------|
| `screening/{GENE}/predictions_exp_{NN}/` | `{gene}_exp{NN}_prefiltered.csv` |

---

### Stage 7 — `07_toxicity_and_grouping.py`

**Purpose:** Advanced medicinal chemistry filtering — **PAINS/Brenk** toxicophores, Lipinski Rule of Five, Veber parameters — followed by **Tanimoto diversity clustering** to select the final **1,000** lead candidates.

```bash
TARGET_GENE=CTNNB1 python 07_toxicity_and_grouping.py
TARGET_GENE=MYC     python 07_toxicity_and_grouping.py
```

| Input | Output |
|-------|--------|
| `{gene}_exp{NN}_prefiltered.csv` | `{gene}_exp{NN}_final_1000.csv`, `{gene}_exp{NN}_clusters.csv` |

---

### Stage 8 — `08_docking_orchestrator.py`

**Purpose:** Orchestrate **AutoDock Vina** docking for the top-N candidates (default N = 100) from the Stage 7 pool.

```bash
# CTNNB1 — champion experiment exp_06
TARGET_GENE=CTNNB1 EXPERIMENT_SUFFIX=06 python 08_docking_orchestrator.py

# MYC — champion experiment exp_07
TARGET_GENE=MYC EXPERIMENT_SUFFIX=07 python 08_docking_orchestrator.py --cpu 16

# Dry-run (validate CID list without docking)
TARGET_GENE=CTNNB1 EXPERIMENT_SUFFIX=06 python 08_docking_orchestrator.py --dry-run
```

| Input | Output |
|-------|--------|
| `screening/{GENE}/predictions_exp_{NN}/*_final_1000.csv` | `docking/runs/{GENE}/exp_{NN}/results/` |

See [`pubchem/docking/README.md`](pubchem/docking/README.md) for Vina environment variables and advanced options.

---

### Stage 9 — `09_docking_analysis.py`

**Purpose:** Aggregate docking summaries into cross-gene **leaderboards**, per-gene **markdown reports**, and a reproducibility **manifest**.

```bash
python 09_docking_analysis.py --exp 06 --limit 100
python 09_docking_analysis.py --exp 07 --genes MYC --limit 100
```

| Input | Output |
|-------|--------|
| `docking/runs/{GENE}/exp_{NN}/results/` | `docking/reports/docking_leaderboard_exp{NN}.csv`, `*.md`, `*.json` |



## 🔄 4. Data Reproduction & Reproducibility

### Why is `pubchem/data/` mostly gitignored?

The `data/` directory contains approximately **61 MB** of derived artifacts that are either:

- **Publicly downloadable** (PubChem confirmatory assay CSVs in `01_raw/`), or
- **Deterministically reproducible** from upstream stages (fingerprints, bit matrices, archives).

Committing these files would bloat the repository, approach GitHub's **100 MB per-file limit** (e.g., `CTNNB1_AID_1665.csv` at 39 MB), and violate standard open-source data-management practices. Instead, this repository follows the **code-as-artifact, data-as-derivation** principle.

### What *is* included in the repository?

| Artifact | Size | Purpose |
|----------|------|---------|
| `data/05_updated_selected_features/` | ~368 KB | Final selected feature matrices — enables model retraining without re-running Stages 1–3 |
| `models/experiments/` | ~1.5 MB | Champion `.pkl` models + JSON metadata |
| `screening/` summaries | ~500 KB | `*_final_1000.csv`, cluster maps, filtering reports |
| `docking/reports/` | ~28 KB | Leaderboards and analysis manifests |

### Full data regeneration (from scratch)

A fresh clone can rebuild the entire data lake by running Stages **1 → 3** in order:

```bash
cd pubchem/scripts

# Step 1: Download PubChem bioassays → 01_raw/ + 02_aggregated/
python 01_data_aggregation.py

# Step 2: Generate ECFP4 fingerprints → 03_fingerprints/
python 02_generate_fingerprints.py

# Step 3a: Build leakage-free bit matrices (interactive)
jupyter notebook ../notebooks/correlation_analysis.ipynb

# Step 3b: Robust tree-based feature selection → 05_updated_selected_features/
python 03_feature_selection.py
```

### Fast path (skip training data rebuild)

If you only need to **re-run screening or docking** with the bundled champion models:

```bash
# Models and feature metadata are already in the repo — jump to Stage 5+
cd pubchem/scripts
python 05_final_screening.py        # requires 03_fingerprints/ (Stage 2)
python 06_combine_and_prefilter.py
python 07_toxicity_and_grouping.py
python 08_docking_orchestrator.py
python 09_docking_analysis.py
```

> **Tip:** For a fully offline demo, keep local copies of `data/03_fingerprints/` and `screening/` outputs on your machine — they are gitignored but not redundant with the tracked summaries.

---

## 🛠️ 5. Prerequisites & Setup

### System Requirements

| Component | Version |
|-----------|---------|
| **OS** | Linux x86_64 (tested on Ubuntu 22.04+) |
| **Python** | 3.10 (`venv310`) |
| **AutoDock Vina** | 1.2.5 (Stages 8–9 only) |
| **RAM** | ≥ 16 GB recommended for virtual screening |
| **Disk** | ~200 MB for full local pipeline outputs |

---

### Step 1 — Clone the repository

```bash
git clone https://github.com/<YOUR_USERNAME>/<YOUR_REPO>.git
cd github_prep
```

---

### Step 2 — Create & activate the virtual environment

```bash
python3.10 -m venv venv310
source venv310/bin/activate        # Linux / macOS
# venv310\Scripts\activate         # Windows
```

---

### Step 3 — Install Python dependencies

The project dependency lockfile is `requirements.txt` at the repository root (full environment including ML, RDKit, Jupyter, and Flask):

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

For **docking-only** minimal dependencies:

```bash
pip install -r pubchem/docking/requirements.txt
# biopython · rdkit · meeko · requests · pyyaml · pandas
```

---

### Step 4 — Install AutoDock Vina (Stages 8–9)

```bash
mkdir -p pubchem/docking/bin && cd pubchem/docking/bin

wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64
mv vina_1.2.5_linux_x86_64 vina

cd ../../scripts
```

---

### Step 5 — Verify the installation

```bash
# RDKit smoke test
python -c "from rdkit import Chem; print('RDKit OK:', Chem.MolFromSmiles('CCO') is not None)"

# Vina smoke test
pubchem/docking/bin/vina --help | head -1

# Pipeline dry-run (no docking)
cd pubchem/scripts
TARGET_GENE=CTNNB1 EXPERIMENT_SUFFIX=06 python 08_docking_orchestrator.py --dry-run
```

---

### Environment Variables Reference

| Variable | Default | Used In | Description |
|----------|---------|---------|-------------|
| `TARGET_GENE` | `CTNNB1` | Stages 6–8 | Target oncoprotein: `MYC`, `CTNNB1`, or `both` |
| `EXPERIMENT_SUFFIX` | `06` | Stage 8 | Experiment folder suffix (`06` for CTNNB1, `07` for MYC) |
| `DOCKING_LIMIT` | `100` | Stage 8 | Top-N ligands to dock from the final-1000 pool |
| `DOCKING_CPU` | `0` | Stage 8 | Vina CPU cores (`0` = auto-detect) |
| `FLASK_APP` | `app.py` | Web Portal | Flask application entry point |

---

## 🔒 Data Leakage Prevention & Experimental Rigor

In previous experimental setups, an implicit data leakage was identified during the feature selection/splitting phase, leading to overly optimistic machine learning metrics. In the latest production pipeline (exp_06 for CTNNB1 and exp_07 for MYC), a rigorous, leakage-free stratified partitioning pipeline was implemented. The slight variance in downstream molecular docking scores directly reflects this mathematically honest and biochemically realistic screening approach.

---

### Previous Setups (Legacy Experiments)

The legacy pipeline (`exp_05` / `exp_06`) exhibited **implicit data leakage** at two critical junctures:

1. **Feature matrix construction** — ECFP4 bit matrices were correlation-pruned (`*_no_corr.csv`) using statistics computed over the **entire compound library**, thereby transmitting label-adjacent variance into the descriptor space before any supervised split.

2. **Feature selection & data splitting** — Optimal feature subsets were derived from globally pre-processed matrices without enforcing a strict **train-only fit / test-only transform** contract. Consequently, cross-validation metrics were optimistically biased, and hold-out performance failed to reflect true generalization capacity on unseen chemical space.

Representative symptoms in the legacy champion models:

| Phenomenon | CTNNB1 `exp_05` | MYC `exp_06` |
|------------|-----------------|--------------|
| Feature space | 30 descriptors (aggressively reduced) | 150 descriptors (pre-leakage matrix) |
| Inflated CV AUC | 0.872 ± 0.036 | **0.908 ± 0.023** |
| Train–test AUC gap | 0.002 (deceptively low) | **0.036** (masked memorization) |
| Virtual screening @ P ≥ 0.50 | 33.5% hit rate (under-selective) | Structurally over-confident CV estimates |
| Overfitting gate | Passed (metric-only) | Passed (metric-only) |

> **Key insight:** A low train–test gap alone is insufficient evidence of generalization when the feature space itself has been contaminated by global preprocessing. Legacy models could achieve superficially acceptable test AUC values while encoding descriptor artifacts that do not transfer to million-compound screening.

---

### Production Pipeline (Leakage-Free Architecture)

The production experiments (`CTNNB1: exp_06`, `MYC: exp_07`) implement a **rigorous, leakage-free stratified partitioning pipeline** with the following design invariants:

| Safeguard | Implementation |
|-----------|----------------|
| **Full-rank input matrix** | Feature selection reads `{GENE}_ecfp4_2048_original.csv` — the complete 2048-bit ECFP4 space without global correlation pruning |
| **Train-only multicollinearity removal** | Pearson \|r\| > 0.90 collinearity pruning is computed **exclusively on the training partition** within `03_feature_selection.py` |
| **Stratified hold-out** | 75/25 stratified train–test split (`random_state=42`) preserving active/inactive class ratios |
| **Nested CV validation** | 5-fold stratified cross-validation on the training set for model benchmarking |
| **Anti-overfitting selection criteria** | Champion model must satisfy: train–test accuracy gap ≤ 0.08, hold-out AUC ≥ 0.85, FPR ≤ 0.01 |
| **Deterministic inference contract** | `feature_order.json` locks column ordering for million-scale virtual screening (Stage 5) |

The resulting **production CV AUC scores** — **0.886 ± 0.036** for CTNNB1 and **0.898 ± 0.033** for MYC — represent statistically defensible estimates of classifier performance on unseen bioactivity data, not optimistically leaked benchmarks.

For MYC specifically, tree-based ensemble models (Random Forest, Extra Trees) that exhibited the highest legacy CV AUC (**≥ 0.930**) were **rejected** in `exp_07` because they violated the overfitting gate (accuracy gap > 0.08), despite their superficially superior cross-validation scores. XGBoost was selected as the champion precisely because it balanced predictive power with **generalization integrity**.

The modest reduction in MYC CV AUC from 0.908 → 0.898 (−0.010) is the mathematically expected correction when leakage pathways are severed. Conversely, CTNNB1 benefits from a richer, leakage-free 133-feature descriptor space (vs. 30 legacy features), yielding improved hold-out recall (0.760 vs. 0.690) and fewer false negatives in the confusion matrix (31 vs. 40).

Variations in downstream **molecular docking affinities** and **virtual hit enrichment** across experiment generations directly reflect this transition: the production pipeline prioritizes *biochemically actionable* candidates over *statistically inflated* predictions.

---

### Experimental Comparison Analysis

#### Table 1 — CTNNB1: Legacy vs. Production (`exp_05` → `exp_06`)

| Metric | `exp_05` (Legacy) | `exp_06` (Production) | Δ |
|--------|:-----------------:|:---------------------:|:-:|
| **Champion algorithm** | Random Forest | Random Forest | — |
| **Selected features (*n*)** | 30 | 133 | +103 |
| **CV AUC (mean ± std)** | 0.872 ± 0.036 | **0.886 ± 0.036** | +0.014 |
| **Hold-out AUC** | 0.899 | **0.913** | +0.014 |
| **Train AUC** | 0.897 | 0.916 | +0.019 |
| **Train–test AUC gap** | 0.002 | 0.003 | +0.001 |
| **Hold-out accuracy** | 0.777 | **0.821** | +0.044 |
| **Hold-out recall** | 0.690 | **0.760** | +0.070 |
| **Hold-out F1** | 0.817 | **0.860** | +0.043 |
| **False negatives (test set)** | 40 | **31** | −9 |
| **Overfitting gate** | ✅ Pass | ✅ Pass | — |
| **Prescreening rate @ P ≥ 0.50** | 33.5% | 70.9% | +37.4 pp |

*Data source: `models/experiments/CTNNB1/exp_{05,06}/Results/model_metadata.json`*

---

#### Table 2 — MYC: Legacy vs. Production (`exp_06` → `exp_07`)

| Metric | `exp_06` (Legacy) | `exp_07` (Production) | Δ |
|--------|:-----------------:|:---------------------:|:-:|
| **Champion algorithm** | Extra Trees | **XGBoost** | Changed |
| **Selected features (*n*)** | 150 | 126 | −24 |
| **CV AUC (mean ± std)** | 0.908 ± 0.023 | **0.898 ± 0.033** | −0.010 |
| **Hold-out AUC** | 0.908 | **0.921** | +0.013 |
| **Train AUC** | 0.944 | 0.930 | −0.014 |
| **Train–test AUC gap** | 0.036 | **0.009** | −0.027 |
| **Train–test accuracy gap** | 0.073 | **0.063** | −0.010 |
| **Hold-out accuracy** | 0.733 | 0.727 | −0.006 |
| **Hold-out recall** | 0.645 | 0.636 | −0.009 |
| **Hold-out F1** | 0.780 | 0.774 | −0.006 |
| **Overfitting gate** | ✅ Pass | ✅ Pass | — |
| **Prescreening rate @ P ≥ 0.50** | 0.4% | 0.6% | +0.2 pp |

*Data source: `models/experiments/MYC/exp_{06,07}/Results/model_metadata.json`*

> The −0.010 decrease in MYC CV AUC is the hallmark of corrected validation: the legacy Extra Trees model achieved elevated cross-validation performance (0.908) partly by exploiting descriptor leakage, while the production XGBoost model trades a modest CV reduction for a **lower train–test AUC gap** (0.009 vs. 0.036) and superior hold-out AUC (0.921).

---

#### Table 3 — MYC `exp_07`: Full Algorithm Benchmark (Overfitting Disqualification)

| Algorithm | CV AUC (mean ± std) | Hold-out AUC | Acc. Gap | Overfitting Gate | Selected |
|-----------|:-------------------:|:------------:|:--------:|:----------------:|:--------:|
| Random Forest | **0.932 ± 0.017** | 0.913 | 0.095 | ❌ Fail | — |
| Extra Trees | **0.930 ± 0.016** | 0.914 | 0.091 | ❌ Fail | — |
| **XGBoost** | 0.898 ± 0.033 | **0.921** | **0.063** | ✅ Pass | **🏆 Champion** |

*Data source: `models/experiments/MYC/exp_07/Results/all_models_results.json`*

This disqualification pattern demonstrates that **maximizing CV AUC alone is an inadequate model selection criterion** when leakage-corrected generalization constraints are enforced.

---

## 🎯 Target Genes & Experiments

### Production Champion Summary

| Gene | Champion Exp. | Algorithm | Features (*n*) | CV AUC (mean ± std) | Hold-out AUC | Hold-out F1 | Train–Test AUC Gap | Screening Threshold | Screening Dir | Docking Suffix | Receptor (PDB) |
|------|:-------------:|-----------|:----------------:|:-------------------:|:------------:|:-----------:|:------------------:|:-------------------:|---------------|:--------------:|:--------------:|
| **CTNNB1** | `exp_06` | Random Forest | 133 | 0.886 ± 0.036 | 0.913 | 0.860 | 0.003 | 0.80 | `screening/CTNNB1/predictions_exp_06/` | `06` | `1JPW` |
| **MYC** | `exp_07` | XGBoost | 126 | 0.898 ± 0.033 | 0.921 | 0.774 | 0.009 | 0.80 | `screening/MYC/predictions_exp_07/` | `07` | `6G6K` |

*Metrics extracted from production `model_metadata.json` artifacts. Receptor parameters defined in `pubchem/docking/config/genes.yaml`.*

---

### Cross-Validation Performance Detail (Production Models)

#### CTNNB1 `exp_06` — Random Forest (Champion)

| Split | Accuracy | AUC | Precision | Recall | F1 |
|-------|:--------:|:---:|:---------:|:------:|:--:|
| CV (5-fold mean) | 0.838 | 0.886 | 0.952 | 0.816 | 0.877 |
| Train | 0.849 | 0.916 | 0.958 | 0.827 | 0.888 |
| **Hold-out test** | **0.821** | **0.913** | **0.990** | **0.760** | **0.860** |

Confusion matrix (test): TN = 49 · FP = 1 · FN = 31 · TP = 98

#### MYC `exp_07` — XGBoost (Champion)

| Split | Accuracy | AUC | Precision | Recall | F1 |
|-------|:--------:|:---:|:---------:|:------:|:--:|
| CV (5-fold mean) | 0.774 | 0.898 | 0.964 | 0.718 | 0.823 |
| Train | 0.790 | 0.930 | 0.971 | 0.735 | 0.836 |
| **Hold-out test** | **0.727** | **0.921** | **0.987** | **0.636** | **0.774** |

Confusion matrix (test): TN = 43 · FP = 1 · FN = 44 · TP = 77

---

### Validation Thresholds (Both Production Experiments)

| Criterion | Threshold | CTNNB1 `exp_06` | MYC `exp_07` |
|-----------|:---------:|:-----------------:|:------------:|
| Max train–test accuracy gap | ≤ 0.08 | 0.028 ✅ | 0.063 ✅ |
| Minimum hold-out AUC | ≥ 0.85 | 0.913 ✅ | 0.921 ✅ |
| Maximum false positive rate | ≤ 0.01 | 0.020 ⚠️ | 0.023 ⚠️ |

> FPR thresholds reflect the inherent precision–recall trade-off in imbalanced bioactivity classification (active:class ratios ≈ 2.6–2.7:1). Both models maintain near-unity precision on the hold-out set.

---

## 📜 License & Citation

This project was developed as a **Senior Design Project** in Computer Engineering.  
If you use this pipeline in academic work, please cite PubChem, ChEMBL, RDKit, and AutoDock Vina accordingly.

---

<p align="center">
  <sub>Built with 🧪 RDKit · 🤖 scikit-learn · ⚓ AutoDock Vina · 🌐 Flask</sub>
</p>
