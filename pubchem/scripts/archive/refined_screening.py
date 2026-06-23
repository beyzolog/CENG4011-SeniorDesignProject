import os
import pandas as pd
import numpy as np
from pycaret.classification import load_model, predict_model
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from joblib import Parallel, delayed, cpu_count
import warnings
from rdkit import RDLogger

# Uyarıları sustur
warnings.filterwarnings("ignore", category=DeprecationWarning)
RDLogger.DisableLog('rdApp.*')

# --- 1. AYARLAR ---
gene = "CTNNB1" #"MYC"
exp_folder = f"../models/experiments/{gene}/exp_01/Results"
model_path = os.path.join(exp_folder, "model_1/model_1_finalize_model")
feature_list_path = os.path.join(exp_folder, f"{gene}_feature_list.txt")
big_data_path = "/home/selahattin/PubChem/PubChem_Compounds/cid_smiles_fingerprint.csv" 
output_dir = f"../screening/{gene}/predictions"

os.makedirs(output_dir, exist_ok=True)

# --- 2. YARDIMCI FONKSİYONLAR ---

def smiles_to_morgan_filtered(smiles, feature_indices, generator):
    """SMILES'ı Morgan'a çevirir ve sadece modelin beklediği bitleri döndürür."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Fonksiyon içinden gelen generator'ı kullanıyoruz
            fp = generator.GetFingerprint(mol)
            fp_list = np.array(list(fp.ToBitString())).astype(int)
            return fp_list[feature_indices]
        return None
    except:
        return None

def process_subchunk(sub_df, model, feature_indices, feature_names):
    """Küçük bir veri parçasını işler ve tahmin üretir."""
    
    # KRİTİK DÜZELTME: Generator nesnesini burada (worker içinde) tanımlıyoruz.
    # Bu sayede PicklingError hatasının önüne geçiyoruz.
    local_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    
    # Morgan parmak izlerini üret (local_generator'ı gönderiyoruz)
    fps = sub_df['Canonical_SMILES'].apply(lambda x: smiles_to_morgan_filtered(x, feature_indices, local_generator))
    
    valid_mask = fps.notnull()
    if not valid_mask.any():
        return pd.DataFrame()
    
    valid_fps = np.stack(fps[valid_mask].values)
    X = pd.DataFrame(valid_fps, columns=feature_names)
    
    preds = predict_model(model, data=X, verbose=False)
    
    results = sub_df[valid_mask].copy()
    results['Score'] = preds['prediction_score'].values
    results['Label'] = preds['prediction_label'].values
    
    return results[results['Label'] == 1]

# --- 3. ANA TARAMA DÖNGÜSÜ ---
def run_screening():
    print(f"🚀 {gene} için tarama başlatılıyor...")
    
    model = load_model(model_path)
    with open(feature_list_path, "r") as f:
        feature_names = [line.strip() for line in f if line.strip()]
    
    feature_indices = [int(name) for name in feature_names]

    reader = pd.read_csv(big_data_path, usecols=["CID", "Canonical_SMILES"], chunksize=1000000)

    for i, df_chunk in enumerate(reader):
        print(f"📦 Chunk {i+1} işleniyor (1 Milyon molekül)...")
        
        n_cores = cpu_count()
        sub_chunks = np.array_split(df_chunk, n_cores)
        
        # Paralel işlem başlarken model ve feature_indices gönderilir
        results_list = Parallel(n_jobs=n_cores)(
            delayed(process_subchunk)(sub, model, feature_indices, feature_names) for sub in sub_chunks
        )
        
        final_hits = pd.concat(results_list, ignore_index=True)
        if not final_hits.empty:
            out_path = os.path.join(output_dir, f"filtered_chunk_{i+1}.csv")
            final_hits.to_csv(out_path, index=False)
            print(f"✅ {len(final_hits)} potansiyel hit kaydedildi: {out_path}")
        else:
            print(f"⚠️ Chunk {i+1} içinde aktif molekül bulunamadı.")

if __name__ == "__main__":
    run_screening()