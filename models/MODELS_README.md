# Trained Models

Model `.pkl` files are **not tracked in Git** (listed in `.gitignore`) because they are large binary files that are fully reproducible from the data.

## Generate models

Run the following command from the repository root to train and save all nine model files:

```bash
python code/variant_proximity_predictor.py --train
```

This will create the following files in this directory:

| File | Description |
|------|-------------|
| `scaler.pkl` | StandardScaler fitted on the 70-D training features |
| `umap_model.pkl` | Supervised UMAP embedding (target_weight=0.5, n_components=2) |
| `rf_clf.pkl` | Random Forest classifier (n_estimators=300, random_state=42) |
| `et_clf.pkl` | Extra Trees classifier (n_estimators=300, random_state=42) |
| `iso_forest.pkl` | Isolation Forest novelty detector (contamination=0.05, n_estimators=200) |
| `cluster_stats.pkl` | Per-class centroid and regularised inverse covariance matrix |
| `label_map.pkl` | Integer-to-variant-name mapping |
| `embedding.pkl` | 2-D UMAP coordinates for all 716 training sequences |
| `y_train.pkl` | Training labels (for UMAP background plot) |

Training takes approximately 2–5 minutes on a standard laptop CPU.
All random seeds are fixed; output is fully deterministic.
