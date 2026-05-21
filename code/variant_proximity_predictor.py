"""
variant_proximity_predictor.py
==============================
Future SARS-CoV-2 Variant Characterisation via Feature-Space Proximity

Authors : Arka Ketan Banerjee, Riasha Pal, Anasua Sarkar
Purpose : Given a new variant's genomic sequences (.fasta), this script
          (a) projects it into the trained 5-class feature space,
          (b) reports a distance-weighted probability over known variants,
          (c) flags genuine novelty (far from ALL known clusters), and
          (d) raises a recombinant warning (equidistant between two clusters).

Usage
-----
  # 1 -- Train once and save models:
      python variant_proximity_predictor.py --train

  # 2 -- Predict a new variant from a FASTA file:
      python variant_proximity_predictor.py --predict path/to/new_variant.fasta --label "XYZ"

  # 3 -- Train + immediately run a demo prediction on held-out samples:
      python variant_proximity_predictor.py --train --demo

Notes
-----
  * Does NOT modify COVID19_Classifier_Fixed.py or any existing file.
  * Requires: numpy, pandas, scikit-learn, umap-learn, matplotlib, joblib
  * Install missing deps:  pip install umap-learn --break-system-packages
"""

import argparse
import os
import re
import warnings
from itertools import product
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import mahalanobis
from sklearn.ensemble import ExtraTreesClassifier, IsolationForest, RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import umap

warnings.filterwarnings("ignore")

# ---- Paths -------------------------------------------------------------------
HERE      = Path(__file__).parent          # .../code/
ROOT      = HERE.parent                    # .../E:/claude/
DATA_DIR  = ROOT / "datasets"
MODEL_DIR = ROOT / "models"
FIG_DIR   = ROOT / "figures"

CODON_CSV = DATA_DIR / "Combined_5class_Codon.csv"
TR_CSV    = DATA_DIR / "NCBI_TR_Fixed_Dataset.csv"   # cols: mono di tri quad pent hex

MODEL_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

# ---- Class metadata ----------------------------------------------------------
CLASS_LABELS = {1: "Alpha", 2: "Mu", 3: "Delta", 4: "Omicron", 5: "Beta"}
CLASS_COLORS = {
    "Alpha":   "#4a9eda",
    "Mu":      "#e74c3c",
    "Delta":   "#27ae60",
    "Omicron": "#e67e22",
    "Beta":    "#8e44ad",
    "NEW":     "#000000",
}
N_CLASSES    = 5
RANDOM_STATE = 42

# All 64 codons in alphabetical order (index 0-63 matches dataset columns '0'-'63')
_BASES = ["A", "C", "G", "T"]
CODONS = ["".join(t) for t in product(_BASES, repeat=3)]   # AAA ... TTT
TR_COLS = ["mono", "di", "tri", "quad", "pent", "hex"]

# ---- Novelty / recombinant thresholds ----------------------------------------
# IsolationForest.decision_function is 0-centred: >0 normal, <0 anomalous
NOVELTY_SCORE_THRESHOLD = 0.0
# If ratio of 2nd-closest / closest Mahalanobis distance is below this,
# the variant may be a recombinant of the two nearest clusters.
RECOMBINANT_RATIO_MAX   = 1.40


# ==============================================================================
# 1. DATA LOADING
# ==============================================================================

def load_combined_features():
    """
    Merge codon-usage (64 features) and tandem-repeat (6 features) datasets.

    Returns
    -------
    df : full DataFrame with Encode + 70 feature columns
    X  : (716, 70) float array
    y  : (716,)    int array  (class labels 1-5)
    """
    codon = pd.read_csv(CODON_CSV)
    tr    = pd.read_csv(TR_CSV)

    # Rename codon feature columns to 'c_0' ... 'c_63'
    codon_feat_cols = [str(i) for i in range(64)]
    codon = codon.rename(columns={c: f"c_{c}" for c in codon_feat_cols})

    feature_cols = [f"c_{i}" for i in range(64)] + TR_COLS
    df = pd.concat(
        [codon.reset_index(drop=True),
         tr[TR_COLS].reset_index(drop=True)],
        axis=1
    )

    X = df[feature_cols].values.astype(float)
    y = df["Encode"].values.astype(int)
    return df, X, y


# ==============================================================================
# 2. FEATURE EXTRACTION FROM FASTA
# ==============================================================================

def _parse_fasta(fasta_path):
    """Minimal FASTA parser -- returns list of uppercase sequence strings."""
    sequences, current = [], []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current).upper())
                current = []
            else:
                current.append(line)
    if current:
        sequences.append("".join(current).upper())
    return sequences


def extract_codon_features(fasta_path):
    """
    Count codon usage across all sequences in a FASTA file (reading-frame 0).
    Returns a (64,) int array aligned with CODONS list (alphabetical order).
    """
    seqs        = _parse_fasta(fasta_path)
    counts      = np.zeros(64, dtype=int)
    codon_index = {c: i for i, c in enumerate(CODONS)}

    for seq in seqs:
        for start in range(0, len(seq) - 2, 3):
            codon = seq[start:start + 3]
            if codon in codon_index:
                counts[codon_index[codon]] += 1
    return counts


def _count_tandem_repeats(seq, unit_len):
    """Count non-overlapping tandem repeats of a given unit length in seq."""
    pattern = r"(.{" + str(unit_len) + r"})\1+"
    return sum(len(m.group()) // unit_len - 1
               for m in re.finditer(pattern, seq))


def extract_tr_features(fasta_path):
    """
    Compute tandem-repeat features (mono ... hex) across all sequences.
    Returns a (6,) int array aligned with TR_COLS.
    """
    seqs   = _parse_fasta(fasta_path)
    counts = np.zeros(6, dtype=int)
    for seq in seqs:
        for i, unit_len in enumerate(range(1, 7)):   # 1=mono ... 6=hex
            counts[i] += _count_tandem_repeats(seq, unit_len)
    return counts


def extract_all_features(fasta_path):
    """
    Full 70-feature vector (64 codon + 6 TR) from a FASTA file.
    Returns a (70,) float array ready for inference.
    """
    codon = extract_codon_features(fasta_path).astype(float)
    tr    = extract_tr_features(fasta_path).astype(float)
    return np.concatenate([codon, tr])


# ==============================================================================
# 3. TRAINING
# ==============================================================================

def train_and_save(verbose=True):
    """
    Train all components on the full 5-class dataset and persist to MODEL_DIR.

    Saved artefacts
    ---------------
    scaler.pkl          StandardScaler
    umap_model.pkl      Supervised UMAP (2-D); supports .transform() on new data
    rf_clf.pkl          RandomForestClassifier
    et_clf.pkl          ExtraTreesClassifier
    iso_forest.pkl      IsolationForest for novelty detection
    cluster_stats.pkl   {class_id: {centroid, inv_cov, label, n}}
    label_map.pkl       CLASS_LABELS dict
    embedding.pkl       (716, 2) training embedding for background plots
    y_train.pkl         (716,) labels for background plots
    """
    if verbose:
        print("Loading combined feature dataset ...")
    _, X, y = load_combined_features()

    # Scale
    scaler   = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Supervised UMAP
    if verbose:
        print("Fitting supervised UMAP (this may take ~30 s) ...")
    umap_model = umap.UMAP(
        n_components=2,
        n_neighbors=20,
        min_dist=0.1,
        metric="euclidean",
        target_metric="categorical",
        target_weight=0.5,
        random_state=RANDOM_STATE,
        transform_seed=RANDOM_STATE,
    )
    umap_model.fit(X_scaled, y=y)
    embedding = umap_model.transform(X_scaled)

    # Classifiers (80/20 split for held-out evaluation)
    if verbose:
        print("Training RandomForest and ExtraTrees classifiers ...")
    X_tr, X_te, y_tr, y_te = train_test_split(
        X_scaled, y, test_size=0.2, stratify=y, random_state=RANDOM_STATE
    )
    rf_clf = RandomForestClassifier(n_estimators=300, random_state=RANDOM_STATE)
    et_clf = ExtraTreesClassifier(n_estimators=300, random_state=RANDOM_STATE)
    rf_clf.fit(X_tr, y_tr)
    et_clf.fit(X_tr, y_tr)

    if verbose:
        names = [CLASS_LABELS[c] for c in sorted(np.unique(y_te))]
        print("\n-- ExtraTrees (80/20 split) --")
        print(classification_report(y_te, et_clf.predict(X_te), target_names=names))

    # Isolation Forest
    if verbose:
        print("Fitting Isolation Forest ...")
    iso_forest = IsolationForest(
        contamination=0.05,
        random_state=RANDOM_STATE,
        n_estimators=200,
    )
    iso_forest.fit(X_scaled)

    # Per-class cluster statistics for Mahalanobis distance
    if verbose:
        print("Computing per-class cluster statistics ...")
    cluster_stats = {}
    for cls in sorted(np.unique(y)):
        mask     = y == cls
        X_cls    = X_scaled[mask]
        centroid = X_cls.mean(axis=0)
        cov      = np.cov(X_cls, rowvar=False)
        cov     += np.eye(cov.shape[0]) * 1e-6   # regularise
        try:
            inv_cov = np.linalg.inv(cov)
        except np.linalg.LinAlgError:
            inv_cov = np.linalg.pinv(cov)
        cluster_stats[cls] = {
            "centroid": centroid,
            "inv_cov":  inv_cov,
            "label":    CLASS_LABELS[cls],
            "n":        int(mask.sum()),
        }

    # Save everything
    artefacts = {
        "scaler":        scaler,
        "umap_model":    umap_model,
        "rf_clf":        rf_clf,
        "et_clf":        et_clf,
        "iso_forest":    iso_forest,
        "cluster_stats": cluster_stats,
        "label_map":     CLASS_LABELS,
        "embedding":     embedding,
        "y_train":       y,
    }
    for name, obj in artefacts.items():
        joblib.dump(obj, MODEL_DIR / f"{name}.pkl")

    if verbose:
        print(f"\nAll models saved to {MODEL_DIR}\n")
    return artefacts


def load_models():
    """Load all saved model artefacts from MODEL_DIR."""
    keys = ["scaler", "umap_model", "rf_clf", "et_clf",
            "iso_forest", "cluster_stats", "label_map", "embedding", "y_train"]
    missing = [k for k in keys if not (MODEL_DIR / f"{k}.pkl").exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing model files: {missing}\nRun:  python variant_proximity_predictor.py --train"
        )
    return {k: joblib.load(MODEL_DIR / f"{k}.pkl") for k in keys}


# ==============================================================================
# 4. INFERENCE
# ==============================================================================

def predict_from_features(X_new, models, variant_label="NEW"):
    """
    Full characterisation of a new variant from its 70-feature vector.

    Parameters
    ----------
    X_new         : (70,) float array -- codon (64) + TR (6) features
    models        : artefact dict from load_models() or train_and_save()
    variant_label : display name for the new variant

    Returns
    -------
    result dict with keys:
      variant           : display label
      rf_proba          : {class_label: probability}   RF classifier
      et_proba          : {class_label: probability}   ExtraTrees
      ensemble_proba    : average of rf and et probas
      top_match         : name of closest variant by ensemble probability
      mahal_distances   : {class_label: Mahalanobis distance to centroid}
      closest_by_dist   : name of closest variant by Mahalanobis distance
      novelty_score     : IsolationForest decision_function score (>0=normal)
      is_novel          : True if score < 0 (anomalous / genuinely new)
      recombinant_flag  : True if top-2 Mahalanobis distances are similar
      recombinant_pair  : (name1, name2) if recombinant_flag else None
      umap_coords       : (x, y) 2-D UMAP projection coordinates
    """
    scaler        = models["scaler"]
    umap_model    = models["umap_model"]
    rf_clf        = models["rf_clf"]
    et_clf        = models["et_clf"]
    iso_forest    = models["iso_forest"]
    cluster_stats = models["cluster_stats"]
    label_map     = models["label_map"]

    X_new = np.asarray(X_new, dtype=float).reshape(1, -1)
    X_sc  = scaler.transform(X_new)

    # Classifier probabilities
    classes  = rf_clf.classes_
    rf_prob  = rf_clf.predict_proba(X_sc)[0]
    et_prob  = et_clf.predict_proba(X_sc)[0]
    ens_prob = (rf_prob + et_prob) / 2

    rf_proba  = {label_map[c]: round(float(p), 4) for c, p in zip(classes, rf_prob)}
    et_proba  = {label_map[c]: round(float(p), 4) for c, p in zip(classes, et_prob)}
    ens_proba = {label_map[c]: round(float(p), 4) for c, p in zip(classes, ens_prob)}
    top_match = max(ens_proba, key=ens_proba.get)

    # Mahalanobis distances to each cluster centroid
    mahal = {}
    for cls, stats in cluster_stats.items():
        try:
            dist = mahalanobis(X_sc[0], stats["centroid"], stats["inv_cov"])
        except Exception:
            diff = X_sc[0] - stats["centroid"]
            dist = float(np.sqrt(diff @ diff))
        mahal[stats["label"]] = round(float(dist), 4)

    sorted_mahal    = sorted(mahal.items(), key=lambda x: x[1])
    closest_by_dist = sorted_mahal[0][0]

    # Recombinant check: are the two nearest clusters suspiciously close?
    d1, d2 = sorted_mahal[0][1], sorted_mahal[1][1]
    ratio   = (d2 / d1) if d1 > 0 else float("inf")
    recombinant_flag = ratio < RECOMBINANT_RATIO_MAX
    recombinant_pair = (sorted_mahal[0][0], sorted_mahal[1][0]) if recombinant_flag else None

    # Novelty score (0-centred; negative = anomalous)
    novelty_score = float(iso_forest.decision_function(X_sc)[0])
    is_novel      = novelty_score < NOVELTY_SCORE_THRESHOLD

    # UMAP projection into the trained 2-D embedding
    umap_coords = tuple(float(v) for v in umap_model.transform(X_sc)[0])

    return {
        "variant":          variant_label,
        "label_map":        label_map,
        "rf_proba":         rf_proba,
        "et_proba":         et_proba,
        "ensemble_proba":   ens_proba,
        "top_match":        top_match,
        "mahal_distances":  mahal,
        "closest_by_dist":  closest_by_dist,
        "novelty_score":    novelty_score,
        "is_novel":         is_novel,
        "recombinant_flag": recombinant_flag,
        "recombinant_pair": recombinant_pair,
        "umap_coords":      umap_coords,
    }


def predict_from_fasta(fasta_path, models, variant_label=None):
    """
    End-to-end prediction: FASTA file -> characterisation result dict.
    Extracts codon + TR features then calls predict_from_features().
    """
    fasta_path = Path(fasta_path)
    if variant_label is None:
        variant_label = fasta_path.stem
    print(f"Extracting features from {fasta_path.name} ...")
    X_new = extract_all_features(fasta_path)
    return predict_from_features(X_new, models, variant_label=variant_label)


# ==============================================================================
# 5. REPORTING & VISUALISATION
# ==============================================================================

def print_report(result):
    """Print a human-readable prediction summary to stdout."""
    v = result["variant"]

    print("\n" + "=" * 60)
    print(f"  PREDICTION REPORT  --  Variant: {v}")
    print("=" * 60)

    print("\n  Ensemble probability (RF + ExtraTrees average):")
    for name, p in sorted(result["ensemble_proba"].items(), key=lambda x: -x[1]):
        bar = "#" * int(p * 30)
        print(f"    {name:<10}  {p:.4f}  {bar}")

    print(f"\n  Top probabilistic match  : {result['top_match']}")
    print(f"  Closest by Mahalanobis   : {result['closest_by_dist']}")

    print("\n  Mahalanobis distance to each cluster centroid:")
    for name, d in sorted(result["mahal_distances"].items(), key=lambda x: x[1]):
        print(f"    {name:<10}  {d:.4f}")

    nv = result["novelty_score"]
    nv_tag = "!! NOVEL -- far from all known clusters !!" if result["is_novel"] else "within known space"
    print(f"\n  Novelty score  : {nv:.4f}  ({nv_tag})")

    if result["recombinant_flag"]:
        p1, p2 = result["recombinant_pair"]
        print(f"  !! RECOMBINANT WARNING: {v} sits close to BOTH {p1} and {p2} clusters")
    else:
        print("  Recombinant check : no signal")

    x, y = result["umap_coords"]
    print(f"\n  UMAP 2-D coords  : x={x:.3f}, y={y:.3f}")
    print("=" * 60 + "\n")


def plot_projection(result, models, save=True):
    """
    Plot the trained UMAP embedding (background = all 716 training points)
    with the new variant projected as a large labelled star marker.
    """
    embedding = models["embedding"]   # (716, 2)
    y_train   = models["y_train"]
    label_map = models["label_map"]
    ux, uy    = result["umap_coords"]
    v_label   = result["variant"]

    fig, ax = plt.subplots(figsize=(9, 7))

    for cls in sorted(np.unique(y_train)):
        mask  = y_train == cls
        cname = label_map[cls]
        ax.scatter(
            embedding[mask, 0], embedding[mask, 1],
            c=CLASS_COLORS[cname], label=cname,
            s=20, alpha=0.55, linewidths=0,
        )

    ax.scatter(
        ux, uy,
        c=CLASS_COLORS["NEW"], marker="*",
        s=420, zorder=10, label=f"NEW: {v_label}",
        edgecolors="white", linewidths=0.8,
    )
    ax.annotate(
        f"  {v_label}", (ux, uy),
        fontsize=11, fontweight="bold", color="black",
        xytext=(8, 8), textcoords="offset points",
    )

    ax.set_title(
        f"UMAP Projection -- Variant '{v_label}' in 5-class Feature Space",
        fontsize=13, fontweight="bold",
    )
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")
    ax.legend(loc="best", fontsize=9, framealpha=0.85)
    plt.tight_layout()

    if save:
        out = FIG_DIR / f"umap_projection_{v_label}.png"
        plt.savefig(out, dpi=150)
        print(f"Plot saved -> {out}")
    plt.show()
    plt.close()


# ==============================================================================
# 6. DEMO -- sanity-check on held-out samples
# ==============================================================================

def run_demo(models):
    """
    Take one held-out sample per variant and predict it.
    All 5 should return the correct variant as the top match.
    """
    print("\n" + "=" * 60)
    print("  DEMO -- predicting held-out samples from training data")
    print("=" * 60)
    _, X, y = load_combined_features()

    for cls in sorted(np.unique(y)):
        idx    = np.where(y == cls)[0][-1]   # last sample of each class
        X_new  = X[idx]
        result = predict_from_features(
            X_new, models, variant_label=f"{CLASS_LABELS[cls]}_held_out"
        )
        top     = result["top_match"]
        correct = "PASS" if top == CLASS_LABELS[cls] else "FAIL"
        print(
            f"  [{correct}] True: {CLASS_LABELS[cls]:<10}  Predicted: {top:<10}"
            f"  (p={result['ensemble_proba'][top]:.3f}, "
            f"novelty={result['novelty_score']:.3f}, "
            f"novel_flag={result['is_novel']})"
        )
    print()


# ==============================================================================
# 7. CLI ENTRY POINT
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Predict SARS-CoV-2 variant characteristics "
            "by feature-space proximity."
        )
    )
    parser.add_argument(
        "--train", action="store_true",
        help="Train and save all models (run once before predicting).",
    )
    parser.add_argument(
        "--predict", type=str, default=None,
        help="Path to a new variant .fasta file to characterise.",
    )
    parser.add_argument(
        "--label", type=str, default=None,
        help="Display label for the new variant (default: filename stem).",
    )
    parser.add_argument(
        "--demo", action="store_true",
        help="Run a quick sanity-check on held-out training samples.",
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip UMAP plot generation.",
    )
    args = parser.parse_args()

    if args.train:
        models = train_and_save(verbose=True)
    else:
        models = load_models()

    if args.demo:
        run_demo(models)

    if args.predict:
        result = predict_from_fasta(args.predict, models, variant_label=args.label)
        print_report(result)
        if not args.no_plot:
            plot_projection(result, models, save=True)


if __name__ == "__main__":
    main()
