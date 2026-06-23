import pandas as pd
import os
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
from datetime import datetime

# Ayarlar
genes = ["MYC", "CTNNB1"]
input_folder = "../data/04_raw_bit_matrices/"
output_folder = "../data/archive/tree_importance_data/"

# Klasörü oluştur
os.makedirs(output_folder, exist_ok=True)

for gene in genes:
    print(f"\n[>>>] {gene} için Tree-based Importance (Leakage-Free) başlatılıyor...")
    print(f"    Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Veriyi yükle
    file_path = os.path.join(input_folder, f"{gene}_ecfp4_2048_no_corr.csv")
    df = pd.read_csv(file_path)
    X = df.drop(columns=["class"])
    y = df["class"]
    
    print(f"    Dataset boyutu: {len(df)} samples, {len(X.columns)} features")
    print(f"    Class distribution: Active={sum(y==1)}, Inactive={sum(y==0)}")

    # 1. Train-Test Split (80-20, stratified) - LEAKAGE FIX!
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, stratify=y, random_state=42
    )
    print(f"    Train: {len(X_train)} samples, Test: {len(X_test)} samples")

    # 2. Modeli SADECE train set ile eğit
    print(f"    [*] RandomForest eğitimi başladı (train set only)...")
    model = RandomForestClassifier(
        n_estimators=200, 
        class_weight='balanced', 
        random_state=42, 
        n_jobs=-1
    )
    model.fit(X_train, y_train)

    # 3. Önem puanlarını al (train'den)
    importances = model.feature_importances_
    
    # 4. Eşik değer (Threshold) belirle
    threshold = np.mean(importances) 
    selected_indices = importances > threshold
    selected_features = X.columns[selected_indices]
    
    print(f"    Toplam özellik: {len(X.columns)}")
    print(f"    Seçilen özellik (Mean Threshold={threshold:.6f}): {len(selected_features)}")

    # 5. Train ve Test performansını değerlendir (leakage kontrolü)
    y_train_pred = model.predict(X_train)
    y_test_pred = model.predict(X_test)
    
    train_acc = accuracy_score(y_train, y_train_pred)
    test_acc = accuracy_score(y_test, y_test_pred)
    test_precision = precision_score(y_test, y_test_pred, average='binary')
    test_recall = recall_score(y_test, y_test_pred, average='binary')
    test_f1 = f1_score(y_test, y_test_pred, average='binary')
    
    # AUC için probability predictions
    y_train_proba = model.predict_proba(X_train)[:, 1]
    y_test_proba = model.predict_proba(X_test)[:, 1]
    train_auc = roc_auc_score(y_train, y_train_proba)
    test_auc = roc_auc_score(y_test, y_test_proba)
    
    print(f"    [✓] Train Accuracy: {train_acc:.4f}, AUC: {train_auc:.4f}")
    print(f"    [✓] Test Accuracy: {test_acc:.4f}, AUC: {test_auc:.4f}")
    print(f"    [✓] Test Precision: {test_precision:.4f}, Recall: {test_recall:.4f}, F1: {test_f1:.4f}")
    print(f"    [✓] Overfitting Gap: {abs(train_acc - test_acc):.4f}")

    # 6. Sonuçları kaydet (TÜM DATA ile ama sadece seçilen feature'lar)
    filtered_df = pd.concat([X[selected_features], y], axis=1)
    filtered_df.to_csv(os.path.join(output_folder, f"{gene}_new_tree_filtered.csv"), index=False)
    
    # 7. Detaylı log dosyası oluştur
    log_path = os.path.join(output_folder, f"{gene}_feature_selection_log.txt")
    with open(log_path, "w") as f:
        f.write(f"=== {gene} Feature Selection Report (Leakage-Free) ===\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Dataset Info:\n")
        f.write(f"  Total samples: {len(df)}\n")
        f.write(f"  Active (1): {sum(y==1)}\n")
        f.write(f"  Inactive (0): {sum(y==0)}\n")
        f.write(f"  Original features: {len(X.columns)}\n\n")
        f.write(f"Train-Test Split:\n")
        f.write(f"  Train: {len(X_train)} samples ({len(X_train)/len(df)*100:.1f}%)\n")
        f.write(f"  Test: {len(X_test)} samples ({len(X_test)/len(df)*100:.1f}%)\n\n")
        f.write(f"Feature Selection:\n")
        f.write(f"  Method: RandomForest Feature Importance (train only)\n")
        f.write(f"  Threshold: Mean importance = {threshold:.6f}\n")
        f.write(f"  Selected features: {len(selected_features)} ({len(selected_features)/len(X.columns)*100:.1f}%)\n\n")
        f.write(f"Model Performance (Full Feature Set):\n")
        f.write(f"  Train Accuracy: {train_acc:.4f}\n")
        f.write(f"  Train AUC: {train_auc:.4f}\n")
        f.write(f"  Test Accuracy: {test_acc:.4f}\n")
        f.write(f"  Test AUC: {test_auc:.4f}\n")
        f.write(f"  Test Precision: {test_precision:.4f}\n")
        f.write(f"  Test Recall: {test_recall:.4f}\n")
        f.write(f"  Test F1: {test_f1:.4f}\n")
        f.write(f"  Overfitting Gap: {abs(train_acc - test_acc):.4f}\n\n")
        f.write(f"Top 20 Most Important Features:\n")
        top_20_indices = np.argsort(importances)[-20:][::-1]
        for i, idx in enumerate(top_20_indices, 1):
            f.write(f"  {i}. Bit {X.columns[idx]}: {importances[idx]:.6f}\n")
    
    # En önemli 10 özelliği ekrana yazdır
    top_10 = X.columns[np.argsort(importances)[-10:]].tolist()
    print(f"    En önemli 10 bit: {top_10}")
    print(f"    [✓] Log kaydedildi: {log_path}")

print("\n[FINISH] İşlem başarıyla tamamlandı.")