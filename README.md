# Alignment-Free SARS-CoV-2 Variant Classification Using Overlapping Trinucleotide Frequencies

**Arka Ketan Banerjee**¹ · **Anasua Sarkar**²

¹ Department of Bio-science and Engineering, Jadavpur University, Kolkata-700032, India  
² Department of Computer Science and Engineering, Jadavpur University, Kolkata-700032, India  

Corresponding author: arkaketan1@gmail.com

---

## Overview

This repository contains the full analysis code, datasets, figures, and manuscript for the study:

> *Alignment-Free SARS-CoV-2 Variant Classification Using Overlapping Trinucleotide Frequencies: A Machine Learning Approach*

We classify SARS-CoV-2 variants (Alpha, Delta, Omicron, Beta, Mu) using 64-dimensional overlapping trinucleotide (3-mer) frequency vectors extracted directly from complete genome sequences — no alignment required.

**Key results:**
- Extra Trees: 94.6% ± 3.1% five-class accuracy; 98.4% ± 1.9% three-class accuracy
- Random Forest: 93.4% ± 3.1% five-class; 97.8% ± 1.9% three-class
- Top discriminating features (SHAP): CCG, GCG, CGG trinucleotides
- CpG O/E ANOVA: F = 14.19, p = 3.82e-11 across five variants

---

## Repository Structure

```
.
├── manuscript/
│   └── COVID19_Full_Analysis_v15_Final.docx   # Submission-ready manuscript
├── code/
│   ├── COVID19_Classifier_Fixed.py            # Main classification pipeline
│   ├── variant_proximity_predictor.py         # Mahalanobis proximity & ensemble predictor
│   └── build_manuscript.py                    # Manuscript builder (python-docx)
├── datasets/
│   ├── NCBI_TR_Dataset.csv                    # Five-class trinucleotide features (NCBI, n=716)
│   ├── NCBI_Codon_Dataset.csv                 # Five-class codon features (NCBI, n=716)
│   ├── ml_results.json                        # ML classification results
│   ├── stats_5class.json                      # Five-class statistical summaries
│   └── stats_3class.json                      # Three-class statistical summaries
├── figures/
│   └── final/                                 # All manuscript figures (Fig 1-13)
├── accessions/
│   ├── NCBI_accessions.txt                    # NCBI sequence accession numbers
│   └── GISAID_EPI_ISL_accessions.txt          # GISAID EPI_ISL accession numbers
├── requirements.txt
└── README.md
```

---

## Datasets

### Five-Class Dataset (NCBI, n=716)
Complete SARS-CoV-2 genome sequences retrieved from NCBI Nucleotide:
- Alpha (n=150), Delta (n=200), Omicron (n=200), Beta (n=108), Mu (n=58)
- Accession numbers: `accessions/NCBI_accessions.txt`

### Three-Class Dataset (NCBI + GISAID, n=558)
- Delta (n=250), Mu (n=58), Omicron (n=250)
- NCBI base sequences augmented with 50 additional sequences per class from GISAID
- GISAID EPI_ISL accession numbers: `accessions/GISAID_EPI_ISL_accessions.txt`
- Raw GISAID sequences cannot be redistributed per the GISAID Data Access Agreement

---

## Methods Summary

1. **Feature extraction**: 64-dimensional overlapping trinucleotide frequency vectors from complete genome FASTA sequences
2. **Classifiers**: Random Forest, Gaussian Naive Bayes, Extra Trees (scikit-learn, random_state=42)
3. **Validation**: 10-fold stratified cross-validation; bootstrap 95% CI (1,000 resamples)
4. **Statistics**: McNemar's test with Bonferroni correction (Family A: 18 tests, alpha_adj=0.0028; Family B: 10 CpG t-tests, alpha_adj=0.005)
5. **Interpretability**: SHAP values; CpG O/E ratio ANOVA
6. **Proximity**: Mahalanobis distance with Ledoit-Wolf shrinkage (alpha=0.617 for Mu)
7. **Visualisation**: t-SNE (max_iter=500, perplexity=30) and UMAP (n_neighbors=15, min_dist=0.1)

---

## Requirements

```
pip install -r requirements.txt
```

---

## Usage

```python
# Run the main classification pipeline
python code/COVID19_Classifier_Fixed.py

# Run proximity predictor
python code/variant_proximity_predictor.py
```

---

## Citation

If you use this code or data, please cite:

> Banerjee AK, Sarkar A. Alignment-Free SARS-CoV-2 Variant Classification Using Overlapping Trinucleotide Frequencies: A Machine Learning Approach. 2025.

Archived release: Zenodo DOI: 10.5281/zenodo.20319358

---

## GISAID Acknowledgement

We gratefully acknowledge the originating laboratories responsible for obtaining SARS-CoV-2 specimens and the submitting laboratories where genetic sequence data were generated and shared via the GISAID EpiCoV database.

---

## Licence

This project is licensed under the MIT Licence. See LICENSE for details.
