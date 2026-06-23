import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

# Ayarlar
BASE_DATA_DIR = "../data"
INPUT_DIR = os.path.join(BASE_DATA_DIR, "02_aggregated")
OUTPUT_DIR = os.path.join(BASE_DATA_DIR, "03_fingerprints")
os.makedirs(OUTPUT_DIR, exist_ok=True)

def smiles_to_fp(smiles):
    """SMILES stringini Morgan Fingerprint (ECFP4) bit dizisine çevirir."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # radius=2 (ECFP4), 2048 bit uzunluk (Standart ML formatı)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
            return "".join([str(b) for b in fp])
        return None
    except:
        return None

def process_file_rdkit(filename):
    if not filename.endswith('.csv'): return
    
    # Label belirleme
    label = 1 if '_active.csv' in filename else 0
    
    print(f"[>] {filename} işleniyor...")
    df = pd.read_csv(os.path.join(INPUT_DIR, filename))
    
    # 1. Label ekle
    df['LABEL'] = label
    
    # 2. Fingerprint hesapla (SMILES sütununu kullanıyoruz)
    # PubChem dosyalarındaki SMILES sütunu ismi: 'PUBCHEM_EXT_DATASOURCE_SMILES'
    smiles_col = 'PUBCHEM_EXT_DATASOURCE_SMILES'
    
    if smiles_col in df.columns:
        df['Morgan_Fingerprint_2048'] = df[smiles_col].apply(smiles_to_fp)
        
        # Parmak izi hesaplanamayan hatalı SMILES'ları temizle
        initial_count = len(df)
        df = df.dropna(subset=['Morgan_Fingerprint_2048'])
        final_count = len(df)
        
        df.to_csv(os.path.join(OUTPUT_DIR, filename), index=False)
        print(f"✓ {filename} bitti. (Kayıt: {final_count}/{initial_count} | Kayıp: {initial_count - final_count})")
    else:
        print(f"x Hata: {filename} içinde SMILES kolonu bulunamadı!")

if __name__ == "__main__":
    files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.csv')]
    for f in files:
        process_file_rdkit(f)
    print("\n[!] Tüm işlemler başarıyla tamamlandı.")