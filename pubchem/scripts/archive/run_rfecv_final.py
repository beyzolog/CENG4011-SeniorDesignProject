import pandas as pd
import os
import joblib
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold

# Ayarlar
genes = ["MYC", "CTNNB1"]
input_folder = "../data/04_raw_bit_matrices/"
output_folder = "../data/archive/rfecv_data/"
report_folder = "../data/archive/rfecv_data/results/"

# Klasörleri oluştur
os.makedirs(output_folder, exist_ok=True)
os.makedirs(report_folder, exist_ok=True)

for gene in genes:
    print(f"\n[>>>] {gene} için RFECV (step=1) başlatılıyor...")
    
    # Veriyi yükle
    file_path = os.path.join(input_folder, f"{gene}_ecfp4_2048_no_corr.csv")
    df = pd.read_csv(file_path)
    X = df.drop(columns=["class"])
    y = df["class"]

    # RFECV Kurulumu (Hocanın standardı: step=1)
    estimator = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42, n_jobs=-1)
    cv = StratifiedKFold(n_splits=5)
    
    rfecv = RFECV(estimator=estimator, step=1, cv=cv, scoring='roc_auc', n_jobs=-1, verbose=0)
    rfecv.fit(X, y)

    # Sonuçları Kaydet
    selected_features = X.columns[rfecv.support_]
    filtered_df = pd.concat([X[selected_features], y], axis=1)
    filtered_df.to_csv(os.path.join(output_folder, f"{gene}_rfecv_final.csv"), index=False)

    # Raporu Kaydet
    results = pd.DataFrame({"Feature": X.columns, "Selected": rfecv.support_, "Ranking": rfecv.ranking_})
    results.to_csv(os.path.join(report_folder, f"{gene}_rfecv_results_final.csv"), index=False)
    
    print(f"[✓] {gene} tamamlandı! Optimum özellik sayısı: {rfecv.n_features_}")

print("\n[FINISH] Tüm genler için işlem başarıyla tamamlandı.")