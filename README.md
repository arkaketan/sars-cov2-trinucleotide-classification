# Alignment-Free SARS-CoV-2 Variant Classification Using Overlapping Trinucleotide Frequencies

**Authors:** Arka Ketan Banerjee, Anasua Sarkar

---

## Overview

This repository contains all code, data, and figures for the manuscript:

> **Alignment-Free SARS-CoV-2 Variant Classification Using Overlapping Trinucleotide Frequencies: A Machine Learning Study with Statistical Validation**

We present a systematic evaluation of overlapping trinucleotide (3-mer) frequency analysis for five-class SARS-CoV-2 variant classification (Alpha, Beta, Delta, Mu, Omicron) using 716 complete genome sequences from NCBI. The Extra Trees classifier achieves **94.6% ± 3.1%** five-class and **98.4% ± 1.9%** three-class accuracy under 10-fold stratified cross-validation.

Key contributions:
- 64-dimensional overlapping trinucleotide frequency feature extraction (no alignment, no ORF annotation)
- SHAP-based identification of CpG-containing trinucleotides (CCG, GCG, CGG) as the dominant discriminative signal, consistent with ZAP-mediated CpG suppression biology
- Bootstrap confidence intervals, McNemar pairwise tests, per-class precision-recall curves, and calibration analysis
- Proximity-based novelty detection framework (Mahalanobis distance + Isolation Forest) for flagging sequences outside the known variant distribution

---

## Repository Structure

```
sars-cov2-trinucleotide-classification/
├── README.md
├── LICENSE
├── .gitignore
├── requirements.txt
├── environment.yml
│
├── manuscript/
│   └── COVID19_Full_Analysis_v8.docx      # Full manuscript (submission version)
│
├── code/
│   ├── COVID19_Classifier_Fixed.py         # Main classification pipeline
│   └── variant_proximity_predictor.py      # Proximity-based surveillance framework
│
├── data/
│   ├── Combined_5class_Codon.csv           # 716 × 65 (64 trinucleotide features + class label)
│   ├── NCBI_TR_Fixed_Dataset.csv           # 716 × 7  (6 tandem-repeat features + class label)
│   ├── Combined_3class_Codon.csv           # 558 × 65 (Alpha, Delta, Omicron only)
│   └── Combined_3class_TR.csv             # 558 × 7  (3-class tandem-repeat)
│
├── figures/
│   ├── fig07_umap_proximity.png            # UMAP 5-class embedding with proximity validation
│   ├── fig08_ensemble_proba.png            # Ensemble probability bars (held-out samples)
│   ├── fig09_feature_importance.png        # ET top-30 feature importances
│   ├── fig10_shap_per_class.png            # Per-class SHAP feature importance
│   ├── fig11_pr_curves.png                 # Per-class precision-recall curves
│   ├── fig12_calibration.png              # Calibration reliability diagrams
│   └── fig13_workflow_diagram.png          # Proximity-based characterisation workflow
│
└── models/                                 # Generated at runtime (not tracked in Git)
```

---

## Installation

### Option 1 — pip

```bash
pip install -r requirements.txt
```

### Option 2 — conda (recommended)

```bash
conda env create -f environment.yml
conda activate sars-trinuc
```

Python 3.9+ is required.

---

## Usage

### 1. Classification pipeline

Trains and evaluates Random Forest, Naive Bayes, and Extra Trees classifiers on the trinucleotide and tandem-repeat feature sets. Produces cross-validation results, confusion matrices, and ROC curves.

```bash
# From the repo root
python code/COVID19_Classifier_Fixed.py
```

The script expects `Master_Dataset.csv` (combined 5-class dataset) in the working directory. You can point it at `data/Combined_5class_Codon.csv` by editing the `DATASET_PATH` constant at the top of the file.

**Output files generated:**
- `confusion_matrix_extra_trees.png`
- `roc_curves.png`
- Console summary of 10-fold CV accuracy, macro F1, MCC, FMI

### 2. Proximity-based variant surveillance framework

Trains a UMAP embedding, Mahalanobis distance scorer, and Isolation Forest novelty detector on the 5-class feature set, then evaluates held-out samples.

```bash
# Train all models (saved to models/)
python code/variant_proximity_predictor.py --train

# Run held-out validation demo (5 samples, one per variant)
python code/variant_proximity_predictor.py --demo

# Predict from a new FASTA file
python code/variant_proximity_predictor.py --predict my_sequence.fasta --label "Unknown"

# Predict without generating the UMAP plot
python code/variant_proximity_predictor.py --predict my_sequence.fasta --no-plot
```

**Data paths** are configured at the top of `variant_proximity_predictor.py`:
```python
CODON_CSV = "data/Combined_5class_Codon.csv"
TR_CSV    = "data/NCBI_TR_Fixed_Dataset.csv"
MODEL_DIR = "models/"
```

**Demo output (held-out single-sequence validation):**

| True Variant | Top Match | Ensemble p | Novelty Score | Novel / Recombinant |
|---|---|---|---|---|
| Alpha   | Alpha   | 0.855 | −0.001 | Yes* / No |
| Mu      | Mu      | 0.848 | −0.246 | Yes** / No |
| Delta   | Delta   | 0.975 | +0.146 | No / No |
| Omicron | Omicron | 0.983 | +0.114 | No / No |
| Beta    | Beta    | 0.892 | +0.154 | No / No |

---

## Data

### Feature CSV format

**`Combined_5class_Codon.csv`** — overlapping trinucleotide frequencies

| Column | Description |
|--------|-------------|
| `Encode` | Class label: 1=Alpha, 2=Mu, 3=Delta, 4=Omicron, 5=Beta |
| `c_0` … `c_63` | Relative frequency of each of the 64 trinucleotides (AAA … TTT) |

**`NCBI_TR_Fixed_Dataset.csv`** — tandem-repeat counts

| Column | Description |
|--------|-------------|
| `Encode` | Class label (same encoding as above) |
| `mono` … `hex` | Count of tandem repeats with unit length 1–6 |

### Class distribution

| Variant | Label | Sequences |
|---------|-------|-----------|
| Alpha (B.1.1.7)     | 1 | 200 |
| Mu (B.1.621)        | 2 | 8   |
| Delta (B.1.617.2)   | 3 | 200 |
| Omicron (BA.1)      | 4 | 200 |
| Beta (B.1.351)      | 5 | 108 |

All 716 sequences were retrieved from NCBI Nucleotide in May 2025. NCBI accession numbers are listed in Supplementary Table S1 of the manuscript.

### Reproducing the dataset from NCBI

```python
from Bio import Entrez
Entrez.email = "your@email.com"
# Use organism-specific search terms per variant (see manuscript Section 2.1)
```

---

## Figures

| Figure | File | Description |
|--------|------|-------------|
| Fig. 7  | `fig07_umap_proximity.png` | Supervised UMAP embedding; X = cluster centroid; ★ = held-out validation samples |
| Fig. 8  | `fig08_ensemble_proba.png` | RF + ET ensemble probability bars for 5 held-out samples |
| Fig. 9  | `fig09_feature_importance.png` | Extra Trees top-30 feature importances (blue = trinucleotide, red = TR) |
| Fig. 10 | `fig10_shap_per_class.png` | Per-class mean |SHAP| for top 15 features |
| Fig. 11 | `fig11_pr_curves.png` | Per-class precision-recall curves with average precision scores |
| Fig. 12 | `fig12_calibration.png` | Reliability diagrams for RF and ET |
| Fig. 13 | `fig13_workflow_diagram.png` | Proximity-based variant characterisation workflow |

Figures 1–6 (topological polar maps, UMAP/t-SNE projections, confusion matrices, bootstrap CI chart, McNemar table) are generated directly by `COVID19_Classifier_Fixed.py`.

---

## Reproducing All Results

```bash
# 1. Install dependencies
conda env create -f environment.yml && conda activate sars-trinuc

# 2. Run main classification pipeline
python code/COVID19_Classifier_Fixed.py

# 3. Train proximity models
python code/variant_proximity_predictor.py --train

# 4. Run held-out validation demo
python code/variant_proximity_predictor.py --demo
```

All random seeds are fixed (`random_state=42`) for full reproducibility.

---

## Citation

If you use this code or data, please cite:

```bibtex
@article{banerjee2025sars,
  title   = {Alignment-Free {SARS-CoV-2} Variant Classification Using Overlapping
             Trinucleotide Frequencies: A Machine Learning Study with Statistical Validation},
  author  = {Banerjee, Arka Ketan and Sarkar, Anasua},
  journal = {[Journal name — to be updated upon acceptance]},
  year    = {2025},
  doi     = {10.5281/zenodo.20319358}
}
```

The archived code release is also available on Zenodo:
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20319358.svg)](https://doi.org/10.5281/zenodo.20319358)

---

## Licence

This project is released under t