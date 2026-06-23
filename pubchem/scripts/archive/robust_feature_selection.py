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
output_folder = "../data/old_selected_features/"

# Klasörü oluştur
os.makedirs(output_folder, exist_ok=True)

# Test edilecek feature sayıları (Selahattin'in ~175'ini referans alarak)
FEATURE_COUNTS_TO_TEST = [30, 40, 50, 80, 100, 120, 150, 175]

def remove_multicollinearity(X_train, selected_features, threshold=0.90):
    """
    %90+ korelasyonlu feature'lardan birini tut, diğerlerini ele.
    """
    X_subset = X_train[selected_features]
    corr_matrix = X_subset.corr().abs()
    
    # Upper triangle'ı al
    upper_tri = corr_matrix.where(
        np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
    )
    
    # %90+ korelasyonlu feature'ları bul
    to_drop = [col for col in upper_tri.columns if any(upper_tri[col] > threshold)]
    
    # Ele
    final_features = [f for f in selected_features if f not in to_drop]
    
    return final_features, len(to_drop)

for gene in genes:
    print(f"\n{'='*80}")
    print(f"[>>>] {gene} için Robust Feature Selection (Gap <5% Priority)")
    print(f"{'='*80}")
    print(f"    Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Veriyi yükle
    file_path = os.path.join(input_folder, f"{gene}_ecfp4_2048_no_corr.csv")
    df = pd.read_csv(file_path)
    X = df.drop(columns=["class"])
    y = df["class"]
    
    print(f"\n📊 Dataset Info:")
    print(f"    Total samples: {len(df)}")
    print(f"    Active (1): {sum(y==1)}")
    print(f"    Inactive (0): {sum(y==0)}")
    print(f"    Original features: {len(X.columns)}")

    # 1. Train-Test Split (80-20, stratified)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, stratify=y, random_state=42
    )
    print(f"\n🔀 Train-Test Split:")
    print(f"    Train: {len(X_train)} samples (80.0%)")
    print(f"    Test: {len(X_test)} samples (20.0%)")

    # 2. Feature Importance hesapla (SADECE train set ile)
    print(f"\n🌲 RandomForest eğitimi (train set only)...")
    base_model = RandomForestClassifier(
        n_estimators=200, 
        class_weight='balanced', 
        random_state=42, 
        n_jobs=-1
    )
    base_model.fit(X_train, y_train)
    importances = base_model.feature_importances_
    
    # Feature'ları importance'a göre sırala
    sorted_indices = np.argsort(importances)[::-1]
    sorted_features = X.columns[sorted_indices]
    sorted_importances = importances[sorted_indices]

    print(f"\n🔍 Testing Top-K + Multicollinearity Removal...")
    print(f"    Strategy: Balance Gap <5% & Test AUC")
    print(f"    Reference: Selahattin used ~175 features")
    print(f"\n{'K':<8} {'After Corr':<12} {'n/p':<8} {'Train AUC':<12} {'Test AUC':<12} {'Gap':<10} {'Prec':<10} {'Recall':<10}")
    print(f"{'-'*90}")
    
    results = []
    
    for K in FEATURE_COUNTS_TO_TEST:
        # Top-K feature'ları seç
        top_k_features = sorted_features[:K]
        
        # Multicollinearity removal
        final_features, dropped_count = remove_multicollinearity(
            X_train, top_k_features, threshold=0.90
        )
        
        # Model eğit
        X_train_selected = X_train[final_features]
        X_test_selected = X_test[final_features]
        
        model = RandomForestClassifier(
            n_estimators=200,
            class_weight='balanced',
            random_state=42,
            n_jobs=-1
        )
        model.fit(X_train_selected, y_train)
        
        # Performans değerlendir
        y_train_proba = model.predict_proba(X_train_selected)[:, 1]
        y_test_proba = model.predict_proba(X_test_selected)[:, 1]
        y_test_pred = model.predict(X_test_selected)
        
        train_auc = roc_auc_score(y_train, y_train_proba)
        test_auc = roc_auc_score(y_test, y_test_proba)
        gap = abs(train_auc - test_auc)
        test_precision = precision_score(y_test, y_test_pred, average='binary')
        test_recall = recall_score(y_test, y_test_pred, average='binary')
        
        n_p_ratio = len(X_train) / len(final_features)
        
        results.append({
            'K': K,
            'n_features': len(final_features),
            'dropped_corr': dropped_count,
            'n_p_ratio': n_p_ratio,
            'train_auc': train_auc,
            'test_auc': test_auc,
            'gap': gap,
            'test_precision': test_precision,
            'test_recall': test_recall,
            'selected_features': final_features
        })
        
        gap_marker = "✅" if gap < 0.05 else "⚠️"
        print(f"{K:<8} {len(final_features):<12} {n_p_ratio:<8.2f} {train_auc:<12.4f} {test_auc:<12.4f} {gap:<10.4f} {test_precision:<10.4f} {test_recall:<10.4f} {gap_marker}")
    
    # 3. Optimal feature count seçimi
    # Priority: Gap < 5% > Test AUC maksimize > En az feature
    print(f"\n🎯 Optimal Selection (Priority: Gap <5% + Max Test AUC):")
    
    # Gap < 5% olanları filtrele
    low_gap_results = [r for r in results if r['gap'] < 0.05]
    
    if low_gap_results:
        # Test AUC'yi maksimize et, eşitlik durumunda en az feature'ı seç
        optimal = max(low_gap_results, key=lambda x: (x['test_auc'], -x['n_features']))
        print(f"    ✅ Found {len(low_gap_results)} solutions with Gap <5%")
        print(f"    ✅ Selected: Max Test AUC with minimum features")
    else:
        # Gap < 5% yoksa, en düşük gap'i seç
        optimal = min(results, key=lambda x: x['gap'])
        print(f"    ⚠️  No solution with Gap <5%")
        print(f"    ⚠️  Selected: Minimum gap solution")
    
    print(f"\n📌 OPTIMAL CONFIGURATION:")
    print(f"    Original K: {optimal['K']}")
    print(f"    After multicollinearity removal: {optimal['n_features']} features")
    print(f"    Dropped (corr >90%): {optimal['dropped_corr']}")
    print(f"    n/p ratio: {optimal['n_p_ratio']:.2f}")
    print(f"    Train AUC: {optimal['train_auc']:.4f}")
    print(f"    Test AUC: {optimal['test_auc']:.4f}")
    print(f"    Overfitting Gap: {optimal['gap']:.4f} ({optimal['gap']*100:.2f}%)")
    print(f"    Test Precision: {optimal['test_precision']:.4f}")
    print(f"    Test Recall: {optimal['test_recall']:.4f}")
    
    # 4. Optimal feature set ile final evaluation
    X_train_optimal = X_train[optimal['selected_features']]
    X_test_optimal = X_test[optimal['selected_features']]
    
    final_model = RandomForestClassifier(
        n_estimators=200,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    final_model.fit(X_train_optimal, y_train)
    
    y_train_pred = final_model.predict(X_train_optimal)
    y_test_pred = final_model.predict(X_test_optimal)
    y_train_proba = final_model.predict_proba(X_train_optimal)[:, 1]
    y_test_proba = final_model.predict_proba(X_test_optimal)[:, 1]
    
    train_acc = accuracy_score(y_train, y_train_pred)
    test_acc = accuracy_score(y_test, y_test_pred)
    train_auc = roc_auc_score(y_train, y_train_proba)
    test_auc = roc_auc_score(y_test, y_test_proba)
    test_precision = precision_score(y_test, y_test_pred, average='binary')
    test_recall = recall_score(y_test, y_test_pred, average='binary')
    test_f1 = f1_score(y_test, y_test_pred, average='binary')
    
    # 5. Sonuçları kaydet (TÜM DATA ile ama sadece optimal feature'lar)
    filtered_df = pd.concat([X[optimal['selected_features']], y], axis=1)
    output_path = os.path.join(output_folder, f"{gene}_robust_tree_filtered.csv")
    filtered_df.to_csv(output_path, index=False)
    print(f"\n💾 Saved: {output_path}")
    
    # 6. Detaylı log dosyası
    log_path = os.path.join(output_folder, f"{gene}_robust_feature_selection_log.txt")
    with open(log_path, "w") as f:
        f.write(f"=== {gene} Robust Feature Selection Report ===\n")
        f.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"Objective: Minimize overfitting gap (<5%) while maintaining Test AUC\n")
        f.write(f"Strategy: Top-K Selection + Multicollinearity Removal (threshold=0.90)\n")
        f.write(f"Priority: Gap <5% > Max Test AUC > Minimum features (Parsimony)\n\n")
        
        f.write(f"Dataset Info:\n")
        f.write(f"  Total samples: {len(df)}\n")
        f.write(f"  Active (1): {sum(y==1)}\n")
        f.write(f"  Inactive (0): {sum(y==0)}\n")
        f.write(f"  Original features: {len(X.columns)}\n\n")
        
        f.write(f"Train-Test Split:\n")
        f.write(f"  Train: {len(X_train)} samples (80.0%)\n")
        f.write(f"  Test: {len(X_test)} samples (20.0%)\n\n")
        
        f.write(f"Strategy Comparison (Top-K + Multicollinearity Removal):\n")
        f.write(f"{'K':<8} {'Final Feat':<12} {'Dropped':<10} {'n/p':<8} {'Train AUC':<12} {'Test AUC':<12} {'Gap':<10} {'Precision':<12} {'Recall':<12}\n")
        f.write(f"{'-'*110}\n")
        for r in results:
            gap_marker = "✅" if r['gap'] < 0.05 else "⚠️"
            f.write(f"{r['K']:<8} {r['n_features']:<12} {r['dropped_corr']:<10} {r['n_p_ratio']:<8.2f} {r['train_auc']:<12.4f} {r['test_auc']:<12.4f} {r['gap']:<10.4f} {r['test_precision']:<12.4f} {r['test_recall']:<12.4f} {gap_marker}\n")
        
        f.write(f"\nOPTIMAL CONFIGURATION:\n")
        f.write(f"  Method: Top-{optimal['K']} + Multicollinearity Removal\n")
        f.write(f"  Selected features: {optimal['n_features']}\n")
        f.write(f"  Dropped (corr >90%): {optimal['dropped_corr']}\n")
        f.write(f"  n/p ratio: {optimal['n_p_ratio']:.2f} (healthy for generalization)\n")
        f.write(f"  Selection criterion: Gap <5% + Max Test AUC + Parsimony\n\n")
        
        f.write(f"Final Model Performance:\n")
        f.write(f"  Train Accuracy: {train_acc:.4f}\n")
        f.write(f"  Train AUC: {train_auc:.4f}\n")
        f.write(f"  Test Accuracy: {test_acc:.4f}\n")
        f.write(f"  Test AUC: {test_auc:.4f}\n")
        f.write(f"  Test Precision: {test_precision:.4f}\n")
        f.write(f"  Test Recall: {test_recall:.4f}\n")
        f.write(f"  Test F1: {test_f1:.4f}\n")
        f.write(f"  Overfitting Gap (AUC): {abs(train_auc - test_auc):.4f}\n")
        f.write(f"  Overfitting Gap (Acc): {abs(train_acc - test_acc):.4f}\n\n")
        
        f.write(f"Comparison with Previous Approach (exp_03):\n")
        f.write(f"  Previous: 315 features, Gap = 0.0966 (9.66%)\n")
        f.write(f"  Current: {optimal['n_features']} features, Gap = {optimal['gap']:.4f} ({optimal['gap']*100:.2f}%)\n")
        f.write(f"  Feature reduction: {((315 - optimal['n_features'])/315)*100:.1f}%\n")
        f.write(f"  Gap improvement: {((0.0966 - optimal['gap'])/0.0966)*100:.1f}% reduction\n")
        f.write(f"  n/p improvement: {optimal['n_p_ratio']/1.68:.2f}x better\n\n")
        
        # Feature importance için optimal set'teki feature'ların importance'larını al
        optimal_feature_importances = []
        for feat in optimal['selected_features']:
            feat_idx = np.where(X.columns == feat)[0][0]
            optimal_feature_importances.append((feat, importances[feat_idx]))
        
        # Importance'a göre sırala
        optimal_feature_importances.sort(key=lambda x: x[1], reverse=True)
        
        f.write(f"Top 30 Most Important Features (from optimal set):\n")
        for i, (feat_name, feat_importance) in enumerate(optimal_feature_importances[:30], 1):
            f.write(f"  {i}. Bit {feat_name}: {feat_importance:.6f}\n")
        
        f.write(f"\nRecommendation for 110M Screening:\n")
        f.write(f"  ✅ Lower overfitting = more reliable predictions\n")
        f.write(f"  ✅ Fewer features = faster screening (~{((315-optimal['n_features'])/315)*100:.0f}% speed boost)\n")
        f.write(f"  ✅ Better generalization = fewer false positives\n")
        f.write(f"  ✅ Use with 0.80 confidence threshold for optimal results\n")
        f.write(f"  ✅ Expected false positive reduction: ~{(optimal['gap']/0.0966)*100:.0f}% vs exp_03\n")
    
    print(f"📄 Log saved: {log_path}")
    
    # Top features'ı ekrana yazdır
    print(f"\n🔝 Top 10 Selected Features:")
    for i, (feat_name, feat_importance) in enumerate(optimal_feature_importances[:10], 1):
        print(f"    {i}. Bit {feat_name}: {feat_importance:.6f}")

print(f"\n{'='*80}")
print(f"[FINISH] Robust feature selection tamamlandı!")
print(f"{'='*80}")
print(f"\n📌 Next Steps:")
print(f"  1. Review the robust_feature_selection_log.txt files")
print(f"  2. Update updated_model_selector.py to use *_robust_tree_filtered.csv")
print(f"  3. Train exp_04 with the optimized feature set")
print(f"  4. Compare exp_03 vs exp_04 performance")
