"""
CTNNB1 Final Model Selector - Production Grade Pipeline
========================================================
Strict anti-overfitting and column order preservation
Based on lessons learned from MYC and CTNNB1 experiments

Date: 2026-04-13
"""

import os
import json
import pandas as pd
import numpy as np
import gc
from datetime import datetime
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import (
    accuracy_score, roc_auc_score, precision_score, 
    recall_score, f1_score, confusion_matrix
)
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from xgboost import XGBClassifier
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

GENE = "MYC" # veya MYC
INPUT_FOLDER = "../data/05_updated_selected_features/"
OUTPUT_ROOT = "../models/experiments/"
RANDOM_STATE = 42

# Strict hyperparameters - Anti-overfitting focus
# WHY THESE VALUES:
# These conservative parameters prevent memorization
HYPERPARAMETERS = {
    'extra_trees': {
        'n_estimators': 200,  # More trees = more stable predictions
        'max_depth': 10,  # CRITICAL: Limit tree depth to prevent memorization
                          # - Unlimited depth (None) caused 12% overfitting
                          # - Depth 10 = max 1024 leaf nodes, prevents over-specialization
        'min_samples_split': 25,  # Need 25 samples to split a node
                                  # - Higher = simpler trees = less overfitting
                                  # - 25/536 train samples = ~4.7% minimum
        'min_samples_leaf': 15,  # Each leaf must have 15 samples
                                 # - Prevents tiny leaves that memorize outliers
                                 # - 15/536 = ~2.8% of training data per leaf
        'max_features': 'sqrt',  # Use sqrt(30) ≈ 5.5 features per split
                                 # - Feature subsampling reduces correlation between trees
                                 # - Prevents overfitting to specific feature combinations
        'bootstrap': True,  # Sample with replacement
        'max_samples': 0.8,  # Use 80% of data for each tree
                             # - Increases diversity between trees
                             # - Reduces overfitting
        'class_weight': 'balanced',  # Auto-adjust for 72%/28% imbalance
                                     # - Prevents "predict everything as active"
        'random_state': RANDOM_STATE,
        'n_jobs': -1  # Use all CPU cores
    },
    'random_forest': {
        'n_estimators': 200,
        'max_depth': 10,  # Same rationale as ExtraTrees
        'min_samples_split': 25,
        'min_samples_leaf': 15,
        'max_features': 'sqrt',
        'max_samples': 0.8,
        'class_weight': 'balanced',
        'random_state': RANDOM_STATE,
        'n_jobs': -1
    },
    'xgboost': {
        'n_estimators': 150,  # Fewer trees than RF (XGBoost learns faster)
        'max_depth': 8,  # Even shallower than RF (XGBoost more prone to overfitting)
        'learning_rate': 0.05,  # CRITICAL: Slow learning prevents overfitting
                                # - 0.05 = each tree contributes only 5% to final prediction
                                # - Slower = more robust = better generalization
        'subsample': 0.8,  # Use 80% of samples per tree (like max_samples)
        'colsample_bytree': 0.8,  # Use 80% of features per tree
                                  # - Double randomization (samples + features)
        'min_child_weight': 10,  # Minimum sum of instance weights in a leaf
                                 # - Higher = more conservative = less overfitting
        'gamma': 0.1,  # Minimum loss reduction to make a split
                       # - Higher = fewer splits = simpler trees
        'reg_alpha': 0.5,  # L1 regularization (Lasso)
                           # - Pushes feature weights toward zero
                           # - Feature selection effect
        'reg_lambda': 1.0,  # L2 regularization (Ridge)
                            # - Penalizes large weights
                            # - Smoother predictions
        'scale_pos_weight': 1.0,  # Will be calculated as n_inactive/n_active
                                  # - Handles class imbalance
        'random_state': RANDOM_STATE,
        'n_jobs': -1,
        'eval_metric': 'logloss'  # Optimization metric
    }
}

# Validation thresholds
# WHY THESE SPECIFIC VALUES:
# Based on MYC and CTNNB1 experimental results
MAX_OVERFITTING_GAP = 0.08  # 8% maximum accuracy gap between train and test
                             # - CTNNB1 exp_04/05: 12.65% gap → model failed in screening
                             # - MYC exp_05: ~6% gap → successful screening (0.12% hit rate)
                             # - 8% is conservative threshold for production use

MIN_TEST_AUC = 0.85  # 85% minimum test AUC
                     # - AUC <0.85: Poor discriminative power
                     # - AUC 0.85-0.90: Good
                     # - AUC >0.90: Excellent
                     # - CTNNB1 exp_04: 0.8549 AUC but failed due to overfitting
                     # - Need both good AUC AND low overfitting

MAX_FALSE_POSITIVE_RATE = 0.01  # 1% maximum on passive molecules
                                 # - At 0.80 threshold on 1000 passive molecules
                                 # - 1% = 10 false positives acceptable
                                 # - On 110M screening: 1% = 1.1M false positives (too many!)
                                 # - But combined with 0.80 threshold, actual rate will be lower
                                 # - This catches overly aggressive models BEFORE screening

# ============================================================================
# LOGGING UTILITIES
# ============================================================================

class Logger:
    def __init__(self, log_file):
        self.log_file = log_file
        self.start_time = datetime.now()
        
    def log(self, message, level="INFO"):
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_msg = f"[{timestamp}] [{level}] {message}"
        print(log_msg)
        with open(self.log_file, 'a') as f:
            f.write(log_msg + "\n")
    
    def section(self, title):
        separator = "=" * 80
        self.log(separator)
        self.log(f"  {title}")
        self.log(separator)

# ============================================================================
# COLUMN ORDER MANAGEMENT
# ============================================================================

def save_feature_order(feature_names, output_path):
    """Save feature order as JSON for strict enforcement during screening
    
    WHY THIS EXISTS:
    - PyCaret and sklearn can reorder columns during preprocessing
    - Screening script builds DataFrames from CSV column order
    - Column mismatch = wrong features fed to model = garbage predictions
    - This JSON locks the EXACT order model expects
    
    HOW IT AFFECTS DATA:
    - Training: Saves original CSV column order BEFORE any transformation
    - Screening: Must reorder DataFrame columns to match this exact order
    """
    feature_order = {
        'features': feature_names,  # List of column names in EXACT order
        'count': len(feature_names),  # Validation: screening must have same count
        'created_at': datetime.now().isoformat()  # Audit trail
    }
    with open(output_path, 'w') as f:
        json.dump(feature_order, f, indent=2)
    return feature_order

def load_feature_order(json_path):
    """Load feature order from JSON"""
    with open(json_path, 'r') as f:
        return json.load(f)

# ============================================================================
# STRATIFIED CROSS-VALIDATION
# ============================================================================

def stratified_cv_evaluation(model, X, y, n_folds=5, logger=None):
    """
    Perform stratified k-fold cross-validation
    
    WHY STRATIFIED:
    - Regular K-Fold can create imbalanced folds (e.g., fold with 90% active)
    - Stratified ensures each fold has same class ratio as full dataset
    - Example: If dataset is 70% active, each fold will be ~70% active
    
    HOW IT AFFECTS DATA:
    - Input: 715 samples (72% active, 28% inactive)
    - Each fold: ~143 samples with SAME 72%/28% ratio
    - More reliable performance estimates
    
    Returns: mean metrics and fold-wise results
    """
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=RANDOM_STATE)
    
    fold_results = []
    
    for fold, (train_idx, val_idx) in enumerate(skf.split(X, y), 1):
        X_train_fold, X_val_fold = X.iloc[train_idx], X.iloc[val_idx]
        y_train_fold, y_val_fold = y.iloc[train_idx], y.iloc[val_idx]
        
        # Train model
        model.fit(X_train_fold, y_train_fold)
        
        # Predictions
        y_pred = model.predict(X_val_fold)
        y_pred_proba = model.predict_proba(X_val_fold)[:, 1]
        
        # Metrics
        metrics = {
            'fold': fold,
            'accuracy': accuracy_score(y_val_fold, y_pred),
            'auc': roc_auc_score(y_val_fold, y_pred_proba),
            'precision': precision_score(y_val_fold, y_pred),
            'recall': recall_score(y_val_fold, y_pred),
            'f1': f1_score(y_val_fold, y_pred)
        }
        fold_results.append(metrics)
        
        if logger:
            logger.log(f"  Fold {fold}: AUC={metrics['auc']:.4f}, "
                      f"Acc={metrics['accuracy']:.4f}, "
                      f"Prec={metrics['precision']:.4f}")
    
    # Calculate mean and std
    mean_metrics = {
        'mean_accuracy': np.mean([r['accuracy'] for r in fold_results]),
        'mean_auc': np.mean([r['auc'] for r in fold_results]),
        'mean_precision': np.mean([r['precision'] for r in fold_results]),
        'mean_recall': np.mean([r['recall'] for r in fold_results]),
        'mean_f1': np.mean([r['f1'] for r in fold_results]),
        'std_auc': np.std([r['auc'] for r in fold_results]),
        'std_accuracy': np.std([r['accuracy'] for r in fold_results])
    }
    
    return mean_metrics, fold_results

# ============================================================================
# PRE-SCREENING TEST ON PASSIVE MOLECULES
# ============================================================================

def pre_screening_test(model, feature_order, logger):
    """
    Test model on random passive molecules to detect false positive tendency
    
    WHY THIS EXISTS:
    - CTNNB1 exp_03 had 57% hit rate (68M false positives!)
    - Root cause: Model too aggressive, predicts everything as active
    - This test catches aggressive models BEFORE 110M screening
    
    HOW IT WORKS:
    - Generate 1000 random binary fingerprints (realistic passive molecules)
    - Fingerprint density: 7% (typical for drug-like molecules)
    - If model predicts >1% as active at threshold 0.80 → REJECT MODEL
    
    HOW IT AFFECTS DATA:
    - Model passes: Safe to use for 110M screening
    - Model fails: NOT saved, prevents wasting 3 hours on bad screening
    
    Returns: (passed, hit_rate, details)
    """
    logger.section("PRE-SCREENING TEST: Random Passive Molecules")
    
    # Generate 1000 random binary fingerprints (simulating passive molecules)
    # In production, you would load actual passive molecules from PubChem
    n_test = 1000
    n_features = len(feature_order['features'])
    
    logger.log(f"Generating {n_test} random passive molecule fingerprints...")
    logger.log(f"Features: {n_features}")
    
    # Random binary fingerprints (sparse, like real molecules)
    # Typical fingerprint density: 5-10%
    np.random.seed(RANDOM_STATE)
    X_passive = np.random.binomial(1, 0.07, size=(n_test, n_features))
    X_passive_df = pd.DataFrame(X_passive, columns=feature_order['features'])
    
    # Predict
    y_pred = model.predict(X_passive_df)
    y_pred_proba = model.predict_proba(X_passive_df)[:, 1]
    
    # Count hits at different thresholds
    hits_at_50 = np.sum(y_pred_proba >= 0.50)
    hits_at_70 = np.sum(y_pred_proba >= 0.70)
    hits_at_80 = np.sum(y_pred_proba >= 0.80)
    
    hit_rate_50 = hits_at_50 / n_test
    hit_rate_70 = hits_at_70 / n_test
    hit_rate_80 = hits_at_80 / n_test
    
    logger.log(f"Results on {n_test} passive molecules:")
    logger.log(f"  Hits at threshold 0.50: {hits_at_50} ({hit_rate_50*100:.2f}%)")
    logger.log(f"  Hits at threshold 0.70: {hits_at_70} ({hit_rate_70*100:.2f}%)")
    logger.log(f"  Hits at threshold 0.80: {hits_at_80} ({hit_rate_80*100:.2f}%)")
    
    # Check if model passes (using 0.80 threshold for screening)
    passed = hit_rate_80 <= MAX_FALSE_POSITIVE_RATE
    
    if passed:
        logger.log(f"✅ PASSED: Hit rate {hit_rate_80*100:.2f}% <= {MAX_FALSE_POSITIVE_RATE*100}%", "SUCCESS")
    else:
        logger.log(f"❌ FAILED: Hit rate {hit_rate_80*100:.2f}% > {MAX_FALSE_POSITIVE_RATE*100}%", "ERROR")
        logger.log("Model is too aggressive - will not be saved!", "ERROR")
    
    # CRITICAL: Convert numpy types to Python native for JSON serialization
    # np.sum() returns numpy.int64, which cannot be serialized to JSON
    return passed, hit_rate_80, {
        'hits_50': int(hits_at_50),  # numpy.int64 → Python int
        'hits_70': int(hits_at_70),
        'hits_80': int(hits_at_80),
        'rate_50': float(hit_rate_50),  # numpy.float64 → Python float
        'rate_70': float(hit_rate_70),
        'rate_80': float(hit_rate_80)
    }

# ============================================================================
# MAIN TRAINING PIPELINE
# ============================================================================

def train_final_model():
    """Main training pipeline with strict validation"""
    
    # Setup experiment directory
    exp_no = 1
    while os.path.exists(os.path.join(OUTPUT_ROOT, GENE, f"exp_{exp_no:02d}")):
        exp_no += 1
    
    exp_dir = os.path.join(OUTPUT_ROOT, GENE, f"exp_{exp_no:02d}", "Results")
    os.makedirs(exp_dir, exist_ok=True)
    
    # Initialize logger
    log_file = os.path.join(exp_dir, "training_log.txt")
    logger = Logger(log_file)
    
    logger.section(f"{GENE} FINAL MODEL TRAINING - PRODUCTION PIPELINE")
    logger.log(f"Experiment: exp_{exp_no:02d}")
    logger.log(f"Output directory: {exp_dir}")
    
    # ========================================================================
    # 1. LOAD DATA WITH STRICT COLUMN ORDER
    # ========================================================================
    logger.section("STEP 1: DATA LOADING")
    
    file_path = os.path.join(INPUT_FOLDER, f"{GENE}_robust_tree_filtered.csv")
    if not os.path.exists(file_path):
        logger.log(f"ERROR: File not found: {file_path}", "ERROR")
        return
    
    df = pd.read_csv(file_path)
    logger.log(f"Data loaded: {len(df)} samples, {len(df.columns)-1} features")
    
    # CRITICAL: Lock column order BEFORE any processing
    # WHY: This is the MOST IMPORTANT line for preventing column mismatch
    # - CSV columns: ["1097", "322", "1114", "1722", ...]
    # - If we let sklearn/PyCaret touch this, order might change
    # - Screening builds DataFrame with CSV order
    # - Model expects THIS order
    # HOW IT AFFECTS DATA:
    # - This list becomes the "ground truth" for all future operations
    # - Saved to JSON, screening MUST match this exact order
    feature_columns = [col for col in df.columns if col != 'class']
    logger.log(f"Feature columns locked: {len(feature_columns)} features")
    logger.log(f"First 10 features: {feature_columns[:10]}")
    
    # Save feature order
    feature_order_path = os.path.join(exp_dir, "feature_order.json")
    feature_order = save_feature_order(feature_columns, feature_order_path)
    logger.log(f"✅ Feature order saved: {feature_order_path}")
    
    # Class distribution
    class_dist = df['class'].value_counts()
    n_inactive = class_dist[0]
    n_active = class_dist[1]
    imbalance_ratio = n_active / n_inactive
    
    logger.log(f"Class distribution:")
    logger.log(f"  Inactive (0): {n_inactive} ({n_inactive/len(df)*100:.1f}%)")
    logger.log(f"  Active (1): {n_active} ({n_active/len(df)*100:.1f}%)")
    logger.log(f"  Imbalance ratio: {imbalance_ratio:.2f}")
    
    # ========================================================================
    # 2. TRAIN-TEST SPLIT (Stratified)
    # ========================================================================
    logger.section("STEP 2: TRAIN-TEST SPLIT (Stratified)")
    
    X = df[feature_columns]  # Use locked column order
    y = df['class']
    
    # WHY 25% TEST SIZE:
    # - Previous experiments used 20% (143 samples) → too small, high variance
    # - 25% gives ~179 samples → more reliable test metrics
    # - Stratify ensures class balance in both train and test
    # HOW IT AFFECTS DATA:
    # - Train: ~536 samples (75%) for model learning
    # - Test: ~179 samples (25%) for unbiased evaluation
    # - Both maintain 72%/28% active/inactive ratio
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, stratify=y, random_state=RANDOM_STATE
    )
    
    logger.log(f"Train set: {len(X_train)} samples")
    logger.log(f"Test set: {len(X_test)} samples")
    logger.log(f"Train active: {y_train.sum()} ({y_train.sum()/len(y_train)*100:.1f}%)")
    logger.log(f"Test active: {y_test.sum()} ({y_test.sum()/len(y_test)*100:.1f}%)")
    
    # Save data split info
    split_info = {
        'train_size': len(X_train),
        'test_size': len(X_test),
        'train_active': int(y_train.sum()),
        'test_active': int(y_test.sum()),
        'random_state': RANDOM_STATE
    }
    with open(os.path.join(exp_dir, "data_split.json"), 'w') as f:
        json.dump(split_info, f, indent=2)
    
    # ========================================================================
    # 3. TRAIN MODELS WITH STRATIFIED CV
    # ========================================================================
    logger.section("STEP 3: MODEL TRAINING & VALIDATION")
    
    models_to_train = {
        'extra_trees': ExtraTreesClassifier(**HYPERPARAMETERS['extra_trees']),
        'random_forest': RandomForestClassifier(**HYPERPARAMETERS['random_forest']),
        'xgboost': XGBClassifier(**HYPERPARAMETERS['xgboost'])
    }
    
    # Calculate scale_pos_weight for XGBoost
    # WHY: XGBoost doesn't have 'class_weight' like sklearn
    # - scale_pos_weight tells XGBoost how much to weight positive class
    # - Formula: n_inactive / n_active
    # - Example: 199 inactive / 516 active = 0.39
    # HOW IT AFFECTS DATA:
    # - Without this: Model ignores minority class (inactive)
    # - With this: Model balances learning from both classes
    # - Prevents "predict everything as active" problem
    models_to_train['xgboost'].scale_pos_weight = n_inactive / n_active
    logger.log(f"XGBoost scale_pos_weight: {n_inactive / n_active:.2f}")
    
    best_model = None
    best_model_name = None
    best_score = 0
    all_results = {}
    
    for model_name, model in models_to_train.items():
        logger.log(f"\n{'='*60}")
        logger.log(f"Training: {model_name.upper()}")
        logger.log(f"{'='*60}")
        
        # 5-Fold Stratified CV
        logger.log("Performing 5-fold stratified cross-validation...")
        cv_metrics, fold_results = stratified_cv_evaluation(
            model, X_train, y_train, n_folds=5, logger=logger
        )
        
        logger.log(f"\nCross-Validation Results:")
        logger.log(f"  Mean AUC: {cv_metrics['mean_auc']:.4f} ± {cv_metrics['std_auc']:.4f}")
        logger.log(f"  Mean Accuracy: {cv_metrics['mean_accuracy']:.4f} ± {cv_metrics['std_accuracy']:.4f}")
        logger.log(f"  Mean Precision: {cv_metrics['mean_precision']:.4f}")
        logger.log(f"  Mean Recall: {cv_metrics['mean_recall']:.4f}")
        logger.log(f"  Mean F1: {cv_metrics['mean_f1']:.4f}")
        
        # Train on full training set
        logger.log("\nTraining on full training set...")
        model.fit(X_train, y_train)
        
        # Evaluate on train set
        y_train_pred = model.predict(X_train)
        y_train_proba = model.predict_proba(X_train)[:, 1]
        
        train_metrics = {
            'accuracy': accuracy_score(y_train, y_train_pred),
            'auc': roc_auc_score(y_train, y_train_proba),
            'precision': precision_score(y_train, y_train_pred),
            'recall': recall_score(y_train, y_train_pred),
            'f1': f1_score(y_train, y_train_pred)
        }
        
        # Evaluate on test set
        y_test_pred = model.predict(X_test)
        y_test_proba = model.predict_proba(X_test)[:, 1]
        
        test_metrics = {
            'accuracy': accuracy_score(y_test, y_test_pred),
            'auc': roc_auc_score(y_test, y_test_proba),
            'precision': precision_score(y_test, y_test_pred),
            'recall': recall_score(y_test, y_test_pred),
            'f1': f1_score(y_test, y_test_pred)
        }
        
        # Calculate overfitting gap
        # WHY: Overfitting = model memorizes training data, fails on new data
        # - Train accuracy 95%, Test accuracy 80% → 15% gap → OVERFITTING!
        # - Gap <5%: Excellent generalization
        # - Gap 5-8%: Acceptable
        # - Gap >8%: Reject model
        # HOW IT AFFECTS DATA:
        # - High gap → Model won't perform well on 110M screening
        # - Low gap → Model generalizes well, reliable predictions
        acc_gap = abs(train_metrics['accuracy'] - test_metrics['accuracy'])
        auc_gap = abs(train_metrics['auc'] - test_metrics['auc'])
        
        logger.log(f"\nTrain Metrics:")
        logger.log(f"  Accuracy: {train_metrics['accuracy']:.4f}")
        logger.log(f"  AUC: {train_metrics['auc']:.4f}")
        logger.log(f"  Precision: {train_metrics['precision']:.4f}")
        logger.log(f"  Recall: {train_metrics['recall']:.4f}")
        
        logger.log(f"\nTest Metrics:")
        logger.log(f"  Accuracy: {test_metrics['accuracy']:.4f}")
        logger.log(f"  AUC: {test_metrics['auc']:.4f}")
        logger.log(f"  Precision: {test_metrics['precision']:.4f}")
        logger.log(f"  Recall: {test_metrics['recall']:.4f}")
        
        logger.log(f"\nOverfitting Analysis:")
        logger.log(f"  Accuracy Gap: {acc_gap:.4f} ({acc_gap*100:.2f}%)")
        logger.log(f"  AUC Gap: {auc_gap:.4f} ({auc_gap*100:.2f}%)")
        
        # Validation checks
        # WHY THESE THRESHOLDS:
        # 1. Overfitting gap ≤ 8%: Ensures model generalizes to new molecules
        # 2. Test AUC ≥ 85%: Ensures model has good discriminative power
        # HOW IT AFFECTS DATA:
        # - Both must pass for model to be considered "best"
        # - If all 3 models fail → Training stops, no model saved
        # - Prevents deploying unreliable models to 110M screening
        passed_overfitting = acc_gap <= MAX_OVERFITTING_GAP
        passed_min_auc = test_metrics['auc'] >= MIN_TEST_AUC
        
        logger.log(f"\nValidation Checks:")
        logger.log(f"  Overfitting gap <= {MAX_OVERFITTING_GAP}: {'✅ PASS' if passed_overfitting else '❌ FAIL'}")
        logger.log(f"  Test AUC >= {MIN_TEST_AUC}: {'✅ PASS' if passed_min_auc else '❌ FAIL'}")
        
        # Confusion matrix
        cm = confusion_matrix(y_test, y_test_pred)
        logger.log(f"\nConfusion Matrix:")
        logger.log(f"  TN: {cm[0,0]}, FP: {cm[0,1]}")
        logger.log(f"  FN: {cm[1,0]}, TP: {cm[1,1]}")
        
        # Store results
        # CRITICAL: Convert numpy types to Python native types for JSON serialization
        # NumPy bool_ and float64 cannot be directly serialized to JSON
        all_results[model_name] = {
            'cv_metrics': {k: float(v) for k, v in cv_metrics.items()},  # Convert all to float
            'train_metrics': {k: float(v) for k, v in train_metrics.items()},
            'test_metrics': {k: float(v) for k, v in test_metrics.items()},
            'acc_gap': float(acc_gap),  # numpy.float64 → Python float
            'auc_gap': float(auc_gap),
            'passed_overfitting': bool(passed_overfitting),  # numpy.bool_ → Python bool
            'passed_min_auc': bool(passed_min_auc),
            'confusion_matrix': cm.tolist()  # numpy.ndarray → Python list
        }
        
        # Track best model (by test AUC)
        if test_metrics['auc'] > best_score and passed_overfitting and passed_min_auc:
            best_score = test_metrics['auc']
            best_model = model
            best_model_name = model_name
            logger.log(f"\n🏆 New best model: {model_name} (Test AUC: {best_score:.4f})")
        
        # Memory cleanup
        gc.collect()
    
    # ========================================================================
    # 4. SELECT AND VALIDATE BEST MODEL
    # ========================================================================
    logger.section("STEP 4: BEST MODEL SELECTION")
    
    if best_model is None:
        logger.log("❌ ERROR: No model passed validation criteria!", "ERROR")
        logger.log("All models failed overfitting or minimum AUC requirements.", "ERROR")
        return
    
    logger.log(f"✅ Best model selected: {best_model_name}")
    logger.log(f"Test AUC: {best_score:.4f}")
    
    # ========================================================================
    # 5. PRE-SCREENING TEST
    # ========================================================================
    passed_prescreening, hit_rate, prescreening_details = pre_screening_test(
        best_model, feature_order, logger
    )
    
    if not passed_prescreening:
        logger.log("❌ Model FAILED pre-screening test - will NOT be saved!", "ERROR")
        logger.log(f"Hit rate on passive molecules: {hit_rate*100:.2f}% > {MAX_FALSE_POSITIVE_RATE*100}%", "ERROR")
        return
    
    # ========================================================================
    # 6. SAVE MODEL AND METADATA
    # ========================================================================
    logger.section("STEP 6: SAVING MODEL AND METADATA")
    
    # Save model using joblib (more reliable than pickle for sklearn)
    import joblib
    model_path = os.path.join(exp_dir, f"{best_model_name}_final_model.pkl")
    joblib.dump(best_model, model_path)
    logger.log(f"✅ Model saved: {model_path}")
    
    # Save comprehensive metadata
    # WHY: This metadata file allows screening script to validate model compatibility
    # and understand model's expected performance before running 110M screening
    metadata = {
        'gene': GENE,
        'experiment': f"exp_{exp_no:02d}",
        'model_type': best_model_name,  # Which algorithm was selected (extra_trees/rf/xgboost)
        'training_date': datetime.now().isoformat(),
        'data': {
            'total_samples': int(len(df)),  # Convert to int for JSON
            'train_samples': int(len(X_train)),
            'test_samples': int(len(X_test)),
            'n_features': int(len(feature_columns)),
            'imbalance_ratio': float(imbalance_ratio)  # Convert numpy type to float
        },
        'performance': all_results[best_model_name],  # Already converted to native types
        'prescreening': prescreening_details,  # Pre-screening test results on passive molecules
        'hyperparameters': HYPERPARAMETERS[best_model_name],  # Exact params used for reproducibility
        'validation_thresholds': {
            'max_overfitting_gap': float(MAX_OVERFITTING_GAP),
            'min_test_auc': float(MIN_TEST_AUC),
            'max_false_positive_rate': float(MAX_FALSE_POSITIVE_RATE)
        },
        'feature_order_file': 'feature_order.json',  # CRITICAL: Points to column order file
        'screening_threshold': 0.80  # Recommended threshold for 110M screening
    }
    
    metadata_path = os.path.join(exp_dir, "model_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    logger.log(f"✅ Metadata saved: {metadata_path}")
    
    # Save all results
    results_path = os.path.join(exp_dir, "all_models_results.json")
    with open(results_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    logger.log(f"✅ All results saved: {results_path}")
    
    # ========================================================================
    # 7. FINAL SUMMARY
    # ========================================================================
    logger.section("TRAINING COMPLETED SUCCESSFULLY")
    
    elapsed = datetime.now() - logger.start_time
    logger.log(f"Total time: {elapsed}")
    logger.log(f"Best model: {best_model_name}")
    logger.log(f"Test AUC: {best_score:.4f}")
    logger.log(f"Overfitting gap: {all_results[best_model_name]['acc_gap']*100:.2f}%")
    logger.log(f"Pre-screening hit rate: {hit_rate*100:.2f}%")
    logger.log(f"\n✅ Model is ready for 110M screening!")
    logger.log(f"Feature order: {feature_order_path}")
    logger.log(f"Model file: {model_path}")
    logger.log(f"Metadata: {metadata_path}")
    
    print("\n" + "="*80)
    print("🎉 TRAINING PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*80)
    print(f"Experiment: exp_{exp_no:02d}")
    print(f"Best model: {best_model_name}")
    print(f"Test AUC: {best_score:.4f}")
    print(f"Output directory: {exp_dir}")
    print("="*80)

# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    train_final_model()
