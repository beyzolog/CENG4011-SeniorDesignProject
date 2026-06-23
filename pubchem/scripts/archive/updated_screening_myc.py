import os
import pandas as pd
import numpy as np
from pycaret.classification import load_model, predict_model
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from joblib import Parallel, delayed, cpu_count
import warnings
from rdkit import RDLogger
import time
import gc

# --- 0. UYARILARI SUSTUR ---
# Terminal kirliliğini önlemek ve log dosyasını temiz tutmak için
warnings.filterwarnings("ignore", category=DeprecationWarning)
RDLogger.DisableLog('rdApp.*')

# --- 1. AYARLAR ---
gene = "MYC" #"CTNNB1" 
# En son deneyi bulma mantığı (Daha önce konuştuğumuz gibi)
base_exp_path = f"../models/experiments/{gene}"
latest_exp = sorted([d for d in os.listdir(base_exp_path) if d.startswith("exp_")])[-1]

exp_folder = os.path.join(base_exp_path, latest_exp, "Results")
# Hoca model_2 kullanmış, sen en iyi AUC olan model_1 ile devam edebilirsin
model_path = os.path.join(exp_folder, "model_1/model_1_finalize_model")
feature_list_path = os.path.join(exp_folder, f"{gene}_feature_list.txt")

# HOCANIN CHUNKSIZE DEĞERİ (2 Milyon)
CHUNK_SIZE = 2000000 
big_data_path = "/home/selahattin/PubChem/PubChem_Compounds/cid_smiles_fingerprint.csv" 
output_dir = f"../screening/{gene}/predictions_{latest_exp}"

# THRESHOLD AYARI (False Positive Kontrolü)
CONFIDENCE_THRESHOLD = 0.80  # Sadece %80+ confidence skorlu tahminler

os.makedirs(output_dir, exist_ok=True)

# --- 2. İŞLEME FONKSİYONU ---
def process_subchunk(sub_df, model, feature_indices, feature_names):
    # RDKit Generator her worker içinde tanımlanmalı
    local_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    # SMILES -> Fingerprint -> Filtered Bits
    def get_fp(s):
        try:
            mol = Chem.MolFromSmiles(s)
            if mol:
                return np.array(list(local_generator.GetFingerprint(mol).ToBitString())).astype(int)[feature_indices]
            return None
        except: return None

    fps = sub_df['Canonical_SMILES'].apply(get_fp)
    
    valid_mask = fps.notnull()
    if not valid_mask.any(): return pd.DataFrame()
    
    X = pd.DataFrame(np.stack(fps[valid_mask].values), columns=feature_names)
    
    # Tahmin
    preds = predict_model(model, data=X, verbose=False)
    
    results = sub_df[valid_mask].copy()
    results['prediction_score'] = preds['prediction_score'].values
    results['prediction_label'] = preds['prediction_label'].values
    
    # THRESHOLD-BASED FILTERING (0.80 confidence)
    # Sadece yüksek confidence skorlu tahminleri döndür (False Positive kontrolü)
    # ✅ DOĞRU: Hem label=1 (Active) HEM DE score ≥ 0.80
    return results[(results['prediction_label'] == 1) & 
                   (results['prediction_score'] >= CONFIDENCE_THRESHOLD)]

# --- 3. ANA DÖNGÜ ---
def run_screening():
    print(f"🚀 {gene} için {CHUNK_SIZE} chunksize ile tarama başladı...")
    print(f"    Model: {latest_exp}/model_1")
    print(f"    Confidence Threshold: {CONFIDENCE_THRESHOLD}")
    print(f"    CPU Cores: {cpu_count()}")
    
    start_time = time.time()
    model = load_model(model_path)
    
    with open(feature_list_path, "r") as f:
        feature_names = [line.strip() for line in f if line.strip()]
    feature_indices = [int(name) for name in feature_names]
    print(f"    Feature count: {len(feature_indices)}")

    # Hocanın okuduğu sütunlara uyum
    reader = pd.read_csv(big_data_path, usecols=["CID", "Canonical_SMILES"], chunksize=CHUNK_SIZE)

    total_hits = 0
    total_processed = 0
    
    for i, df_chunk in enumerate(reader):
        chunk_start = time.time()
        print(f"\n🔍 Chunk {i + 1} işleniyor...")
        print(f"    Chunk size: {len(df_chunk)} molecules")
        
        # Paralel bölme
        sub_chunks = np.array_split(df_chunk, cpu_count())
        results = Parallel(n_jobs=cpu_count())(
            delayed(process_subchunk)(sub, model, feature_indices, feature_names) for sub in sub_chunks
        )
        
        final_hits = pd.concat(results, ignore_index=True)
        chunk_time = time.time() - chunk_start
        total_processed += len(df_chunk)
        
        if not final_hits.empty:
            out_path = os.path.join(output_dir, f"filtered_chunk_{i+1}.csv")
            final_hits.to_csv(out_path, index=False)
            total_hits += len(final_hits)
            print(f"    ✅ Hits: {len(final_hits)} (Threshold ≥ {CONFIDENCE_THRESHOLD})")
            print(f"    💾 Kaydedildi: {out_path}")
        else:
            print(f"    ⚠️  No hits above threshold {CONFIDENCE_THRESHOLD}")
        
        print(f"    ⏱️  Chunk time: {chunk_time:.2f}s ({len(df_chunk)/chunk_time:.0f} mol/s)")
        print(f"    📊 Total processed: {total_processed:,} | Total hits: {total_hits:,}")
        
        # Memory cleanup
        del df_chunk, sub_chunks, results, final_hits
        gc.collect()
    
    total_time = time.time() - start_time
    print(f"\n🏁 Tarama tamamlandı!")
    print(f"    Total time: {total_time/60:.2f} minutes")
    print(f"    Total processed: {total_processed:,} molecules")
    print(f"    Total hits: {total_hits:,} (Hit rate: {total_hits/total_processed*100:.4f}%)")
    print(f"    Average speed: {total_processed/total_time:.0f} mol/s")

if __name__ == "__main__":
    run_screening()