import os
import time
import pandas as pd
import requests
from Bio import Entrez
from io import BytesIO
from multiprocessing import Pool

# --- YAPILANDIRMA ---
Entrez.email = "beyzayoruk@posta.mu.edu.tr" 
BASE_DIR = "../data"
OUTPUT_DIR = os.path.join(BASE_DIR, "02_aggregated")  # Filtrelenmiş/Birleştirilmiş veriler
RAW_DIR = os.path.join(BASE_DIR, "01_raw")       # Ham verilerin tutulduğu yer
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(RAW_DIR, exist_ok=True)

def process_gene(gene_name):
    gene_name = gene_name.strip()
    if not gene_name: return
    
    print(f"\n[>] {gene_name} işleniyor...")
    
    # 1. ADIM: ESearch ile Güncel Confirmatory Assay'leri Bul
    term = f'"{gene_name}"[Gene Symbol] AND (confirmatory[filt])'
    try:
        handle = Entrez.esearch(db="pcassay", retmax=1000, term=term, idtype="acc")
        record = Entrez.read(handle)
        handle.close()
        aid_list = record["IdList"]
        print(f"    - {gene_name}: {len(aid_list)} adet onaylanmış AID bulundu.")
    except Exception as e:
        print(f"    - {gene_name} arama hatası: {e}")
        return

    active_frames = []
    inactive_frames = []
    gene_raw_path = os.path.join(RAW_DIR, gene_name)
    os.makedirs(gene_raw_path, exist_ok=True)

    # 2. ADIM: İndirme, Yedekleme ve Filtreleme
    for aid in aid_list:
        # SMILES içeren datatable URL'si
        url = f"https://pubchem.ncbi.nlm.nih.gov/assay/pcget.cgi?query=download&record_type=datatable&actvty=all&response_type=save&aid={aid}"
        
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 200:
                # --- YAN İŞLEM: Ham Veriyi İstenen Formatla Kaydet ---
                # Format: GEN_AID_NUMARA.csv (Örn: CTNNB1_AID_1904.csv)
                raw_file_name = f"{gene_name}_AID_{aid}.csv"
                raw_file_path = os.path.join(gene_raw_path, raw_file_name)
                
                with open(raw_file_path, 'wb') as f:
                    f.write(r.content)
                
                # --- ANA İŞLEM: Filtreleme ---
                df = pd.read_csv(BytesIO(r.content), low_memory=False)
                
                # Meta-data temizliği
                df['PUBCHEM_RESULT_TAG'] = pd.to_numeric(df['PUBCHEM_RESULT_TAG'], errors='coerce')
                df = df.dropna(subset=['PUBCHEM_RESULT_TAG']).copy()

                # IC50 Filtrelemesi
                if 'Standard Type' in df.columns:
                    df_ic50 = df[df['Standard Type'].str.contains('IC50', case=False, na=False)].copy()
                    
                    if not df_ic50.empty:
                        df_ic50['PubChem Standard Value'] = pd.to_numeric(df_ic50['PubChem Standard Value'], errors='coerce')
                        
                        # Aktifler
                        act = df_ic50[df_ic50['PUBCHEM_ACTIVITY_OUTCOME'] == "Active"]
                        if not act.empty: active_frames.append(act)
                        
                        # İnaktifler (Outcome Active değilse ve > 10 uM)
                        inact = df_ic50[(df_ic50['PUBCHEM_ACTIVITY_OUTCOME'] != "Active") & 
                                        (df_ic50['PubChem Standard Value'] > 10)]
                        if not inact.empty: inactive_frames.append(inact)
            
            time.sleep(0.3) 
        except Exception as e:
            print(f"      ! {aid} işlenirken hata: {e}")
            continue

    # 3. ADIM: Sonuçları Birleştir ve Kaydet
    if active_frames:
        final_active = pd.concat(active_frames, ignore_index=True, sort=False)
        final_active.to_csv(os.path.join(OUTPUT_DIR, f"{gene_name}_active.csv"), index=False)
    
    if inactive_frames:
        final_inactive = pd.concat(inactive_frames, ignore_index=True, sort=False)
        final_inactive.to_csv(os.path.join(OUTPUT_DIR, f"{gene_name}_inactive.csv"), index=False)

    print(f"[!] {gene_name} bitti. Aktif: {len(active_frames)}, İnaktif: {len(inactive_frames)} dosya birleştirildi.")

if __name__ == "__main__":
    # genes.txt dosyasından gen listesini oku
    genes_file = os.path.join(BASE_DIR, "genes.txt")
    if os.path.exists(genes_file):
        with open(genes_file, "r") as f:
            genes = [line.strip() for line in f if line.strip()]
        
        # 4 paralel işlemle çalıştır
        with Pool(4) as p:
            p.map(process_gene, genes)
    else:
        print(f"Hata: {genes_file} bulunamadı!")