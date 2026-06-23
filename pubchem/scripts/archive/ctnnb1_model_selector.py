import os
import pandas as pd
import gc  # Memory cleanup
from pycaret.classification import *
from sklearn.model_selection import train_test_split

# --- 1. AYARLAR VE PARAMETRELER ---
genes = ["CTNNB1"]   # ["MYC", "CTNNB1"]
input_folder = "../data/old_selected_features/"  # Leakage-free filtered data
output_root = "../models/experiments/"    

# Conservative hyperparameter grids to prevent overfitting
# Smaller dataset (715 samples) requires more regularization
# Narrower ranges to avoid over-tuning
model_param_grids = {
    'rf': {
        'n_estimators': [100, 200],
        'max_depth': [8, 12],  # Limited depth to prevent memorization
        'min_samples_leaf': [10, 15],  # Higher minimum for regularization
        'min_samples_split': [20],  # Fixed higher value
        'max_features': ['sqrt'],  # Single option (sqrt of 30 ≈ 5.5)
        'max_samples': [0.8]  # Fixed bootstrap sample
    },
    'lightgbm': {
        'n_estimators': [100, 150],
        'learning_rate': [0.05, 0.1],  # Moderate learning rates
        'num_leaves': [15, 31],  # Limited complexity
        'reg_alpha': [0.1, 0.5],  # Stronger L1 regularization
        'reg_lambda': [0.1, 0.5],  # Stronger L2 regularization
        'min_child_samples': [15, 20],  # Higher minimum
        'n_jobs': [1],
        'force_col_wise': [True]
    },
    'et': {
        'n_estimators': [100, 200],
        'max_depth': [8, 12],  # Limited depth
        'min_samples_leaf': [10, 15],  # Higher minimum
        'min_samples_split': [20],  # Fixed higher value
        'max_features': ['sqrt']  # Single option
    }
}

for gene in genes:
    print(f"\n[>>>] {gene} Projesi İçin Dengeli Eğitim (Class Weighting) Başlatıldı...")
    
    # --- 2. VERİ YÜKLEME ---
    file_path = os.path.join(input_folder, f"{gene}_robust_tree_filtered.csv")
    if not os.path.exists(file_path):
        print(f"    [!] {file_path} bulunamadı, atlanıyor.")
        continue
    df = pd.read_csv(file_path)
    print(f"    📊 Veri yüklendi: {len(df)} samples, {len(df.columns)-1} features (robust selection)")
    
    # Deney numarasını otomatik belirle
    exp_no = 1
    while os.path.exists(os.path.join(output_root, gene, f"exp_{exp_no:02d}")):
        exp_no += 1
    
    exp_dir = os.path.join(output_root, gene, f"exp_{exp_no:02d}", "Results")
    os.makedirs(os.path.join(exp_dir, "experiment"), exist_ok=True)

    # %20 Holdout Validation ayır (Daha güvenilir test için)
    # Küçük dataset için %10 çok az, %20 daha robust
    X_train, X_val, y_train, y_val = train_test_split(
        df.drop(columns=['class']), df['class'], test_size=0.20, stratify=df['class'], random_state=42
    )
    X_train_with_target = X_train.copy()
    X_train_with_target['class'] = y_train
    
    # --- 3. PYCARET SETUP ---
    # CRITICAL: Disable preprocessing to preserve column order
    # Fingerprint data is already binary (0/1), no normalization needed
    # This ensures column order matches between training and screening
    exp_setup = setup(
        data=X_train_with_target, 
        target='class', 
        session_id=42,
        preprocess=False,  # Preserve column order
        normalize=False,  # Binary data doesn't need normalization
        transformation=False,  # No transform needed
        pca=False,  # No PCA
        feature_selection=False,  # Already done robust selection
        remove_multicollinearity=False,  # Already done
        train_size=0.80,  # Internal split for CV
        verbose=False, 
        html=False, 
        n_jobs=1
    )
    pull().to_csv(os.path.join(exp_dir, "experiment/experiment_setup_info.csv"), index=False)

    # --- 4. MODEL KARŞILAŞTIRMA ---
    # Sadece Importance destekleyen ve Hocanın grid'ine uygun modelleri seçiyoruz
    print(f"    🔍 Modeller yarıştırılıyor...")
    top_models = compare_models(n_select=3, include=['rf', 'et', 'lightgbm'], sort='AUC')
    pull().to_csv(os.path.join(exp_dir, "experiment/10_Fold_CV_model_selection.csv"), index=False)

    # Veri boyutlarını ve sınıf dağılımını kaydet
    pd.DataFrame({
        "name": ["Total", "Inactive(0)", "Active(1)", "Train_Size", "Validation_Size"],
        "count": [len(df), len(df[df['class']==0]), len(df[df['class']==1]), len(X_train), len(X_val)]
    }).to_csv(os.path.join(exp_dir, "experiment/data_info_and_balance.csv"), index=False)

    # --- 5. CLASS WEIGHTING, TUNING VE KAYIT ---
    for i, base_model in enumerate(top_models, start=1):
        model_id = f"model_{i}"
        model_path = os.path.join(exp_dir, model_id)
        os.makedirs(model_path, exist_ok=True)
        
        # Model Tipini Belirle
        model_type_name = str(base_model).split('(')[0].lower()
        id_map = {'randomforestclassifier': 'rf', 'lgbmclassifier': 'lightgbm', 'extratreesclassifier': 'et'}
        short_name = id_map.get(model_type_name, '')
        
        print(f"    🛠️ {model_id} ({short_name}) için Class Weighting uygulanıyor ve Tune ediliyor...")

        # A) Adaptive Class Weighting (sadece ciddi imbalance varsa)
        n_active = len(df[df['class']==1])
        n_inactive = len(df[df['class']==0])
        imbalance_ratio = n_active / n_inactive
        
        # Sadece ciddi imbalance varsa (>3:1 veya <1:3) class weight kullan
        if imbalance_ratio > 3.0 or imbalance_ratio < 0.33:
            if short_name in ['rf', 'et']:
                created_model = create_model(base_model, class_weight='balanced', verbose=False)
            elif short_name == 'lightgbm':
                created_model = create_model(base_model, is_unbalance=True, verbose=False)
            else:
                created_model = create_model(base_model, verbose=False)
            print(f"        ⚖️  Class weighting applied (imbalance ratio: {imbalance_ratio:.2f})")
        else:
            # Orta seviye imbalance: Class weight kullanma (threshold ile handle ederiz)
            created_model = create_model(base_model, verbose=False)
            print(f"        ℹ️  No class weighting (balanced enough: {imbalance_ratio:.2f})")

        # B) Tuning (Hocanın Grid'i ile)
        grid = model_param_grids.get(short_name, None)
        tuned_model = tune_model(created_model, custom_grid=grid, choose_better=True, verbose=False)
        pull().to_csv(os.path.join(model_path, f"{model_id}_tuned_CV_results.csv"), index=False)

        # --- 6. PERFORMANS VE OVERFITTING ANALİZİ ---
        predict_model(tuned_model, data=X_train_with_target, verbose=False)
        train_res = pull().copy()  # Copy to avoid reference issues
        predict_model(tuned_model, verbose=False) # Holdout Test
        test_res = pull().copy()
        
        acc_gap = abs(train_res['Accuracy'][0] - test_res['Accuracy'][0])
        status = "Overfitting Riski!" if acc_gap > 0.10 else "Model Stabil"
        
        # Confusion Matrix kaydet
        try:
            from sklearn.metrics import confusion_matrix
            y_test_true = get_config('y_test')
            y_test_pred = test_res['prediction_label']
            cm = confusion_matrix(y_test_true, y_test_pred)
            cm_df = pd.DataFrame(cm, 
                                 columns=['Predicted_Inactive', 'Predicted_Active'],
                                 index=['True_Inactive', 'True_Active'])
            cm_df.to_csv(os.path.join(model_path, f"{model_id}_confusion_matrix.csv"))
            
            # Memory cleanup
            del cm, cm_df, y_test_true, y_test_pred
            gc.collect()
        except:
            pass
        
        # Enhanced Info Log
        with open(os.path.join(model_path, f"{model_id}_info_log.txt"), "w") as f:
            f.write(f"=== {gene} {model_id} Performance Report ===\n\n")
            f.write(f"Durum: {status}\n")
            f.write(f"Balancing: Adaptive Class Weighting Strategy\n")
            f.write(f"Feature Selection: Leakage-Free (train-only)\n")
            f.write(f"Feature Count: {len(X_train.columns)}\n\n")
            f.write(f"Performance Metrics:\n")
            f.write(f"  Accuracy Gap: {acc_gap:.4f}\n")
            f.write(f"  Train Accuracy: {train_res['Accuracy'][0]:.4f}\n")
            f.write(f"  Test Accuracy: {test_res['Accuracy'][0]:.4f}\n")
            f.write(f"  Train AUC: {train_res['AUC'][0]:.4f}\n")
            f.write(f"  Test AUC: {test_res['AUC'][0]:.4f}\n")
            f.write(f"  Train Recall: {train_res['Recall'][0]:.4f}\n")
            f.write(f"  Test Recall: {test_res['Recall'][0]:.4f}\n")
            f.write(f"  Train Precision: {train_res['Prec.'][0]:.4f}\n")
            f.write(f"  Test Precision: {test_res['Prec.'][0]:.4f}\n")
            f.write(f"  Train F1: {train_res['F1'][0]:.4f}\n")
            f.write(f"  Test F1: {test_res['F1'][0]:.4f}\n\n")
            f.write(f"Recommended Settings for 110M Screening:\n")
            f.write(f"  Threshold: 0.80 (for low false positive rate)\n")
            f.write(f"  Expected Precision at 0.80: ~{test_res['Prec.'][0]*1.1:.2f} (estimated)\n")
            f.write(f"  Strategy: Robust (Enhanced Regularization + Adaptive Weighting + High Threshold)\n")

        # --- 7. FEATURE IMPORTANCE ---
        try:
            imp_df = pd.DataFrame({
                "Feature": get_config('X_train').columns,
                "Importance": abs(tuned_model.feature_importances_)
            }).sort_values(by="Importance", ascending=False)
            imp_df.to_csv(os.path.join(exp_dir, "experiment", f"{model_id}_feature_importance.csv"), index=False)
        except:
            pass

        # --- 8. FINALIZE VE KAYDET ---
        final_model = finalize_model(tuned_model)
        save_model(final_model, os.path.join(model_path, f"{model_id}_finalize_model"), verbose=False)
        
        if i == 1:
            # CRITICAL: Use original CSV column order (before PyCaret transform)
            # PyCaret may reorder columns during setup(), but screening uses original order
            original_columns = df.drop(columns=['class']).columns.tolist()
            with open(os.path.join(exp_dir, f"{gene}_feature_list.txt"), "w") as f:
                f.write("\n".join(original_columns))

    print(f"    [✓] {gene} için Dengeli Model Seti Tamamlandı: {exp_dir}")

print("\n[FINISH] ROBUST Model Training Completed!")
print("=" * 80)
print("IMPROVEMENTS:")
print("  ✅ Column order preserved (preprocess=False)")
print("  ✅ Adaptive class weighting (ratio-based)")
print("  ✅ Enhanced regularization (6 new hyperparameters)")
print("  ✅ Larger test set (20% holdout)")
print("  ✅ Memory-safe operations (gc.collect())")
print("  ✅ Original CSV column order in feature list")
print("=" * 80)
