"""
CTNNB1 Final Screening Script - Compatible with exp_05
=======================================================
Uses feature_order.json for strict column alignment
Compatible with joblib-saved sklearn models

Date: 2026-04-13
"""

import os
import json
import pandas as pd
import numpy as np
from sympy import fps
import joblib
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from joblib import Parallel, delayed, cpu_count
import warnings
from rdkit import RDLogger
import time
import gc

# Suppress warnings
warnings.filterwarnings("ignore")
RDLogger.DisableLog('rdApp.*')

# ============================================================================
# CONFIGURATION
# ============================================================================

GENE = "MYC"
BASE_EXP_PATH = f"../models/experiments/{GENE}"

# Auto-detect latest experiment
latest_exp = sorted([d for d in os.listdir(BASE_EXP_PATH) if d.startswith("exp_")])[-1]
print(f"🔍 Using experiment: {latest_exp}")

EXP_FOLDER = os.path.join(BASE_EXP_PATH, latest_exp, "Results")

# Load model metadata to get model type
metadata_path = os.path.join(EXP_FOLDER, "model_metadata.json")
with open(metadata_path, 'r') as f:
    metadata = json.load(f)

model_type = metadata['model_type']  # e.g., "random_forest"
screening_threshold = metadata['screening_threshold']  # 0.80

# OPTIONAL: Override threshold for testing
screening_threshold = 0.80  # Same as MYC for fair comparison (ACTIVE)


print(f"📊 Model type: {model_type}")
print(f"🎯 Screening threshold: {screening_threshold}")

# Paths
MODEL_PATH = os.path.join(EXP_FOLDER, f"{model_type}_final_model.pkl")
FEATURE_ORDER_PATH = os.path.join(EXP_FOLDER, "feature_order.json")

# Data paths
CHUNK_SIZE = 2_000_000
# PILOT TEST: 444K molecules (COMPLETED - SUCCESS ✅)
#BIG_DATA_PATH = "/home/selahattin/PubChem/PubChem_Compounds/cid_fingerprint_canonical_smiles/Compound_169000001_169500000.csv"
# FULL SCREENING: 119M molecules (ACTIVE - READY TO GO 🚀)
BIG_DATA_PATH = "/home/selahattin/PubChem/PubChem_Compounds/cid_smiles_fingerprint.csv"

OUTPUT_DIR = f"../screening/{GENE}/predictions_{latest_exp}"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# LOAD MODEL AND FEATURE ORDER
# ============================================================================

print(f"\n📦 Loading model...")
model = joblib.load(MODEL_PATH)
print(f"✅ Model loaded: {MODEL_PATH}")

print(f"\n🔒 Loading feature order...")
with open(FEATURE_ORDER_PATH, 'r') as f:
    feature_order = json.load(f)

feature_names = feature_order['features']
feature_count = feature_order['count']
feature_indices = [int(name) for name in feature_names]

print(f"✅ Feature order loaded: {feature_count} features")
print(f"   First 10 features: {feature_names[:10]}")
print(f"   Feature order locked from: {feature_order['created_at']}")

# ============================================================================
# PROCESSING FUNCTION
# ============================================================================

def process_subchunk(sub_df, model, feature_indices, feature_names, threshold):
    """
    Process a sub-chunk of molecules
    
    WHY THIS FUNCTION:
    - Parallel processing for speed
    - Each worker has its own RDKit generator (thread-safe)
    - Strict column order enforcement using feature_names
    
    HOW IT AFFECTS DATA:
    - SMILES → Morgan fingerprint (ECFP4, radius=2, 2048 bits)
    - Extract only the features in feature_indices
    - Build DataFrame with EXACT column order from feature_names
    - Predict using sklearn model (not PyCaret!)
    - Filter by threshold
    """
    # RDKit generator (must be created in each worker)
    local_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    def get_fingerprint(smiles):
        """Convert SMILES to filtered fingerprint vector"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Get full 2048-bit fingerprint
                fp = local_generator.GetFingerprint(mol)
                fp_array = np.array(list(fp.ToBitString())).astype(int)
                # Extract only selected features in EXACT order
                return fp_array[feature_indices]
            return None
        except:
            return None
    
    # Generate fingerprints
    fps = sub_df['Canonical_SMILES'].apply(get_fingerprint)
    
    # Filter valid molecules
    valid_mask = fps.notnull()
    if not valid_mask.any():
        return pd.DataFrame()
    
    # CRITICAL: Build DataFrame with EXACT column order
    # This ensures columns match model's training order
    # CRITICAL FIX FOR XGBOOST: Sütun isimlerini zorla metne (string) çeviriyoruz
    X = pd.DataFrame(
        np.stack(fps[valid_mask].values),
        columns=[str(name) for name in feature_names] # Metinsel hizalama ve kilitleme
    )
    
    # Predict using sklearn model directly (not PyCaret)
    # WHY: exp_05 uses joblib-saved sklearn model, not PyCaret pipeline
    y_pred = model.predict(X)
    y_pred_proba = model.predict_proba(X)[:, 1]  # Probability of class 1 (active)
    
    # Build results
    results = sub_df[valid_mask].copy()
    results['prediction_label'] = y_pred
    results['prediction_score'] = y_pred_proba
    
    # Filter: label=1 (active) AND score >= threshold
    # WHY: Same logic as MYC, prevents false positives
    filtered = results[
        (results['prediction_label'] == 1) & 
        (results['prediction_score'] >= threshold)
    ]
    
    return filtered

# ============================================================================
# MAIN SCREENING LOOP
# ============================================================================

def run_screening():
    """Main screening pipeline"""
    
    print(f"\n{'='*80}")
    print(f"🚀 {GENE} SCREENING STARTED")
    print(f"{'='*80}")
    print(f"Experiment: {latest_exp}")
    print(f"Model: {model_type}")
    print(f"Features: {feature_count}")
    print(f"Threshold: {screening_threshold}")
    print(f"Chunk size: {CHUNK_SIZE:,}")
    print(f"CPU cores: {cpu_count()}")
    print(f"Data: {BIG_DATA_PATH}")
    print(f"Output: {OUTPUT_DIR}")
    print(f"{'='*80}\n")
    
    start_time = time.time()
    
    # CSV reader
    reader = pd.read_csv(
        BIG_DATA_PATH,
        usecols=["CID", "Canonical_SMILES"],
        chunksize=CHUNK_SIZE
    )
    
    total_hits = 0
    total_processed = 0
    
    for i, df_chunk in enumerate(reader, 1):
        chunk_start = time.time()
        print(f"\n{'─'*60}")
        print(f"🔍 Processing Chunk {i}")
        print(f"{'─'*60}")
        print(f"Chunk size: {len(df_chunk):,} molecules")
        
        # Split into sub-chunks for parallel processing
        n_workers = cpu_count()
        sub_chunks = np.array_split(df_chunk, n_workers)
        
        # Parallel processing
        results = Parallel(n_jobs=n_workers)(
            delayed(process_subchunk)(
                sub, model, feature_indices, feature_names, screening_threshold
            ) for sub in sub_chunks
        )
        
        # Combine results
        final_hits = pd.concat(results, ignore_index=True)
        
        # Update counters
        chunk_time = time.time() - chunk_start
        total_processed += len(df_chunk)
        chunk_hits = len(final_hits)
        total_hits += chunk_hits
        
        # Save hits
        if not final_hits.empty:
            out_path = os.path.join(OUTPUT_DIR, f"filtered_chunk_{i}.csv")
            final_hits.to_csv(out_path, index=False)
            print(f"✅ Hits: {chunk_hits:,} (Threshold ≥ {screening_threshold})")
            print(f"💾 Saved: {out_path}")
        else:
            print(f"⚠️  No hits above threshold {screening_threshold}")
        
        # Performance metrics
        chunk_speed = len(df_chunk) / chunk_time
        hit_rate = (chunk_hits / len(df_chunk)) * 100
        
        print(f"⏱️  Chunk time: {chunk_time:.2f}s ({chunk_speed:.0f} mol/s)")
        print(f"📊 Chunk hit rate: {hit_rate:.4f}%")
        print(f"📈 Total processed: {total_processed:,} | Total hits: {total_hits:,}")
        
        # Memory cleanup
        del df_chunk, sub_chunks, results, final_hits
        gc.collect()
    
    # Final summary
    total_time = time.time() - start_time
    overall_hit_rate = (total_hits / total_processed) * 100
    avg_speed = total_processed / total_time
    
    print(f"\n{'='*80}")
    print(f"🏁 SCREENING COMPLETED")
    print(f"{'='*80}")
    print(f"Total time: {total_time/60:.2f} minutes")
    print(f"Total processed: {total_processed:,} molecules")
    print(f"Total hits: {total_hits:,}")
    print(f"Hit rate: {overall_hit_rate:.4f}%")
    print(f"Average speed: {avg_speed:.0f} mol/s")
    print(f"{'='*80}\n")
    
    # Validation check
    if overall_hit_rate > 1.0:
        print("⚠️  WARNING: Hit rate > 1% - Model may be too aggressive!")
        print("   Expected: 0.00-0.05% based on pre-screening test")
    elif overall_hit_rate < 0.001:
        print("✅ EXCELLENT: Very low hit rate - Model is conservative")
    else:
        print("✅ GOOD: Hit rate within expected range")

if __name__ == "__main__":
    run_screening()
