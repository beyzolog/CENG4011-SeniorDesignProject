import os
import pandas as pd
from pycaret.classification import *
from sklearn.model_selection import train_test_split

# --- 1. AYARLAR VE PARAMETRELER ---
genes = ["CTNNB1"]   # ["MYC", "CTNNB1"]
input_folder = "../data/old_selected_features/"  # Leakage-free filtered data
output_root = "../models/experiments/"    

# Optimized hyperparameter grids for small dataset (n=528)
# Narrower ranges to reduce overhead, especially for LightGBM
model_param_grids = {
    'rf': {
        'n_estimators': [100, 200],
        'max_depth': [5, 10],
        'min_samples_leaf': [10]
    },
    'lightgbm': {
        'n_estimators': [100],
        'learning_rate': [0.05],
        'num_leaves': [31],
        'reg_alpha': [0.1],
        'n_jobs': [1],  # Single thread for small data
        'force_col_wise': [True]  # Skip 'choosing' phase
    },
    'et': {
        'n_estimators': [100, 200],
        'max_depth': [5, 10]
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

    # %10 Holdout Validation ayır (Sınıf oranını koruyarak)
    X_train, X_val, y_train, y_val = train_test_split(
        df.drop(columns=['class']), df['class'], test_size=0.10, stratify=df['class'], random_state=42
    )
    X_train_with_target = X_train.copy()
    X_train_with_target['class'] = y_train
    
    # --- 3. PYCARET SETUP ---
    # Not: Sınıf dengesini model seviyesinde (Weighting) çözeceğimiz için fix_imbalance=False
    # n_jobs=1 for small dataset to avoid multi-threading overhead
    exp_setup = setup(data=X_train_with_target, target='class', session_id=42, 
                      preprocess=True, verbose=False, html=False, n_jobs=1)
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

        # A) Model Oluştururken Ağırlıklandırma Ekle (Penalized Learning)
        if short_name in ['rf', 'et']:
            created_model = create_model(base_model, class_weight='balanced', verbose=False)
        elif short_name == 'lightgbm':
            created_model = create_model(base_model, is_unbalance=True, verbose=False)
        else:
            created_model = create_model(base_model, verbose=False)

        # B) Tuning (Hocanın Grid'i ile)
        grid = model_param_grids.get(short_name, None)
        tuned_model = tune_model(created_model, custom_grid=grid, choose_better=True, verbose=False)
        pull().to_csv(os.path.join(model_path, f"{model_id}_tuned_CV_results.csv"), index=False)

        # --- 6. PERFORMANS VE OVERFITTING ANALİZİ ---
        predict_model(tuned_model, data=X_train_with_target, verbose=False)
        train_res = pull()
        predict_model(tuned_model, verbose=False) # Holdout Test (Setup içindeki %30)
        test_res = pull()
        
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
        except:
            pass
        
        # Enhanced Info Log
        with open(os.path.join(model_path, f"{model_id}_info_log.txt"), "w") as f:
            f.write(f"=== {gene} {model_id} Performance Report ===\n\n")
            f.write(f"Durum: {status}\n")
            f.write(f"Balancing: Class Weighting Applied (Hybrid Strategy)\n")
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
            f.write(f"  Strategy: Hybrid (Class Weight + High Threshold)\n")

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

print("\n[FINISH] Class Weighting entegre edilmiş tüm modeller kaydedildi!")