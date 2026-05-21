"""
COVID-19 Variant Classification — Corrected & Reproducible Script
Authors: Arka Ketan Banerjee, Riasha Pal, Anasua Sarkar

Fixes applied vs. original Testing_RandomForest.ipynb:
  1. RandomForestRegressor -> RandomForestClassifier (classification task)
  2. random_state=42 added to ExtraTreesClassifier for reproducibility
  3. Hard-coded absolute paths replaced with relative/configurable paths
  4. Stratified train/test split (ensures class balance in both sets)
  5. Added 10-fold stratified cross-validation for reliable accuracy estimates
  6. Added per-class F1 scores, MCC, FMI, and ROC-AUC
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_validate
from sklearn.metrics import (
    accuracy_score, f1_score, matthews_corrcoef,
    confusion_matrix, ConfusionMatrixDisplay,
    fowlkes_mallows_score, roc_auc_score, RocCurveDisplay,
    classification_report
)
import warnings
warnings.filterwarnings('ignore')

# ── Configuration ──────────────────────────────────────────────────────────────
DATASET_PATH = "Master_Dataset.csv"   # set to your path if running elsewhere
RANDOM_STATE = 42
TEST_SIZE    = 0.20
N_FOLDS      = 10
CLASS_LABELS = {1: "Delta", 2: "Mu", 3: "Omicron"}

# ── Load Data ──────────────────────────────────────────────────────────────────
df = pd.read_csv(DATASET_PATH)
y  = df["Encode"]
X  = df.drop(["Encode"], axis=1)

print(f"Dataset shape : {df.shape}")
print(f"Class counts  :\n{y.value_counts().rename(CLASS_LABELS)}\n")

# ── Single 80/20 Stratified Split ─────────────────────────────────────────────
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y
)

classifiers = {
    "Random Forest": RandomForestClassifier(random_state=RANDOM_STATE),
    "Naive Bayes":   GaussianNB(),
    "Extra Trees":   ExtraTreesClassifier(random_state=RANDOM_STATE),
}

# ── Single-Split Results ───────────────────────────────────────────────────────
print("=" * 65)
print("SINGLE TRAIN/TEST SPLIT RESULTS  (80/20, stratified)")
print("=" * 65)

split_results = {}
fitted_clfs   = {}

for name, clf in classifiers.items():
    clf.fit(X_train, y_train)
    fitted_clfs[name] = clf
    y_pred = clf.predict(X_test)

    acc     = accuracy_score(y_test, y_pred)
    f1_mac  = f1_score(y_test, y_pred, average="macro")
    f1_per  = f1_score(y_test, y_pred, average=None, labels=[1, 2, 3])
    mcc     = matthews_corrcoef(y_test, y_pred)
    fmi     = fowlkes_mallows_score(y_test, y_pred)
    cm      = confusion_matrix(y_test, y_pred, labels=[1, 2, 3])

    split_results[name] = dict(acc=acc, f1_mac=f1_mac, f1_per=f1_per,
                                mcc=mcc, fmi=fmi, cm=cm)
    print(f"\n{name}:")
    print(f"  Accuracy          : {acc:.4f}  ({acc*100:.1f}%)")
    print(f"  Macro F1          : {f1_mac:.4f}")
    print(f"  F1 (Delta/Mu/Omicron): {f1_per.round(4)}")
    print(f"  MCC               : {mcc:.4f}")
    print(f"  FMI               : {fmi:.4f}")
    print(f"  Confusion Matrix  :\n{cm}")

# ── 10-Fold Stratified Cross-Validation ───────────────────────────────────────
print("\n" + "=" * 65)
print(f"{N_FOLDS}-FOLD STRATIFIED CROSS-VALIDATION RESULTS")
print("=" * 65)

cv       = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
scoring  = ["accuracy", "f1_macro", "matthews_corrcoef"]
cv_table = []

for name, clf in classifiers.items():
    fresh = clf.__class__(**clf.get_params())
    scores = cross_validate(fresh, X, y, cv=cv, scoring=scoring)

    row = {
        "Classifier"  : name,
        "Acc Mean"    : scores["test_accuracy"].mean(),
        "Acc Std"     : scores["test_accuracy"].std(),
        "F1 Mean"     : scores["test_f1_macro"].mean(),
        "F1 Std"      : scores["test_f1_macro"].std(),
        "MCC Mean"    : scores["test_matthews_corrcoef"].mean(),
        "MCC Std"     : scores["test_matthews_corrcoef"].std(),
    }
    cv_table.append(row)
    print(f"\n{name}:")
    print(f"  Accuracy  : {row['Acc Mean']:.4f} ± {row['Acc Std']:.4f}  "
          f"({row['Acc Mean']*100:.1f}% ± {row['Acc Std']*100:.1f}%)")
    print(f"  Macro F1  : {row['F1 Mean']:.4f} ± {row['F1 Std']:.4f}")
    print(f"  MCC       : {row['MCC Mean']:.4f} ± {row['MCC Std']:.4f}")

cv_df = pd.DataFrame(cv_table).set_index("Classifier")

# ── Confusion Matrix Plot (Extra Trees — best performer) ──────────────────────
fig, ax = plt.subplots(figsize=(6, 5))
ConfusionMatrixDisplay(
    confusion_matrix=split_results["Extra Trees"]["cm"],
    display_labels=["Delta", "Mu", "Omicron"]
).plot(ax=ax, colorbar=False, cmap="Blues")
ax.set_title("Extra Trees Classifier — Confusion Matrix (80/20 split)")
plt.tight_layout()
plt.savefig("confusion_matrix_extra_trees.png", dpi=150)
plt.close()
print("\nSaved: confusion_matrix_extra_trees.png")

# ── ROC Curve Plot (One-vs-Rest, all three classifiers) ───────────────────────
fig, ax = plt.subplots(figsize=(8, 6))
colors = {"Random Forest": "steelblue", "Naive Bayes": "tomato", "Extra Trees": "seagreen"}

for name, clf in fitted_clfs.items():
    if hasattr(clf, "predict_proba"):
        y_prob = clf.predict_proba(X_test)
        auc = roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro")
        # plot macro OvR curve via micro average proxy label
        RocCurveDisplay.from_predictions(
            (y_test == 1).astype(int),
            y_prob[:, 0],
            name=f"{name} (macro AUC≈{auc:.2f})",
            ax=ax,
            color=colors[name],
            plot_chance_level=(name == "Naive Bayes"),
        )

ax.set_title("ROC Curves — One-vs-Rest (Delta class shown, macro AUC in legend)")
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
plt.tight_layout()
plt.savefig("roc_curves.png", dpi=150)
plt.close()
print("Saved: roc_curves.png")

# ── Summary Table ──────────────────────────────────────────────────────────────
print("\n" + "=" * 65)
print("CROSS-VALIDATION SUMMARY TABLE")
print("=" * 65)
summary = cv_df[["Acc Mean", "Acc Std", "F1 Mean", "F1 Std", "MCC Mean", "MCC Std"]]
summary.columns = ["Acc", "±", "F1", "±", "MCC", "±"]
print(summary.round(4).to_string())
print("\nDone.")
