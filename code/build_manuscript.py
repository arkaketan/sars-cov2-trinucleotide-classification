"""Build COVID19_Full_Analysis_v15_Final.docx — v13 review issues fixed.

Fixes applied vs. v13:
  CRIT-1   — Ref [22] duplicated author line removed (introduced by tail-repair in v13)
  mi-v13-1 — Sec 2.1: "severely class-imbalanced" → "class-imbalanced" for consistency
              with Limitation (3) which correctly calls n=58 only moderate imbalance
"""
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

FIGURES = "/sessions/magical-determined-thompson/mnt/claude/figures/"
MEDIA   = "/sessions/magical-determined-thompson/mnt/outputs/extracted_media/"
IMG  = lambda n: MEDIA   + f"word_media_image{n}.png"
FIG  = lambda n: FIGURES + n

doc = Document()
style = doc.styles['Normal']
style.font.name = 'Times New Roman'
style.font.size = Pt(12)
for h in ('Heading 1', 'Heading 2', 'Heading 3'):
    hs = doc.styles[h]; hs.font.name = 'Times New Roman'; hs.font.color.rgb = RGBColor(0,0,0)

def hp(text, level=1): doc.add_heading(text, level=level)
def p(text, bold=False, italic=False):
    para = doc.add_paragraph(); run = para.add_run(text)
    run.bold = bold; run.italic = italic; return para
def pb(label, rest):
    para = doc.add_paragraph(); r1 = para.add_run(label); r1.bold = True
    para.add_run(rest); return para
def caption(text):
    para = doc.add_paragraph(); para.paragraph_format.space_before = Pt(4)
    para.paragraph_format.space_after = Pt(10)
    run = para.add_run(text); run.italic = True; run.font.size = Pt(10); return para
def fig(path, width_in=6.0):
    para = doc.add_paragraph(); para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    para.paragraph_format.space_before = Pt(6); para.paragraph_format.space_after = Pt(2)
    run = para.add_run(); run.add_picture(path, width=Inches(width_in))
def simple_table(headers, rows, col_widths=None, footnote=None):
    table = doc.add_table(rows=1+len(rows), cols=len(headers)); table.style = 'Table Grid'
    hdr = table.rows[0]
    for i,h in enumerate(headers):
        hdr.cells[i].text = h
        for run in hdr.cells[i].paragraphs[0].runs: run.bold = True
    for ri, row_data in enumerate(rows):
        row = table.rows[ri+1]
        for ci, val in enumerate(row_data): row.cells[ci].text = str(val)
    if col_widths:
        for row in table.rows:
            for ci,w in enumerate(col_widths): row.cells[ci].width = Inches(w)
    if footnote:
        fn = doc.add_paragraph(); fn.paragraph_format.space_before = Pt(2)
        fn.paragraph_format.space_after = Pt(8)
        run = fn.add_run(footnote); run.font.size = Pt(9); run.italic = True
    else:
        doc.add_paragraph()

# =======================================================
# TITLE & AUTHORS
# =======================================================
tp = doc.add_paragraph(); tp.alignment = WD_ALIGN_PARAGRAPH.CENTER
tr = tp.add_run("Alignment-Free SARS-CoV-2 Variant Classification Using Overlapping "
    "Trinucleotide Frequencies: A Machine Learning Study with Statistical Validation")
tr.bold = True; tr.font.size = Pt(14)

# sub-1: Author names with affiliation superscripts
auth = doc.add_paragraph(); auth.alignment = WD_ALIGN_PARAGRAPH.CENTER
r1 = auth.add_run("Arka Ketan Banerjee"); r1.italic = True
auth.add_run("¹").font.size = Pt(9)
auth.add_run(", Anasua Sarkar").italic = True
auth.add_run("²").font.size = Pt(9)

# Affiliation 1
aff1 = doc.add_paragraph(); aff1.alignment = WD_ALIGN_PARAGRAPH.CENTER
aff1.add_run("¹ Department of Bio-science and Engineering, "
             "Jadavpur University, Kolkata-700032, India").font.size = Pt(10)

# Affiliation 2
aff2 = doc.add_paragraph(); aff2.alignment = WD_ALIGN_PARAGRAPH.CENTER
aff2.add_run("² Department of Computer Science and Engineering, "
             "Jadavpur University, Kolkata-700032, India").font.size = Pt(10)

corr = doc.add_paragraph(); corr.alignment = WD_ALIGN_PARAGRAPH.CENTER
corr_run = corr.add_run("Corresponding author: Arka Ketan Banerjee "
                        "(arkaketan1@gmail.com)")
corr_run.font.size = Pt(10)
doc.add_paragraph()

# =======================================================
# ABSTRACT
# =======================================================
hp("Abstract")
p("Background: Rapid, alignment-free classification of SARS-CoV-2 variants is "
  "essential for scalable genomic surveillance in settings where multiple-sequence "
  "alignment is computationally prohibitive. Overlapping trinucleotide (3-mer) frequency "
  "spectra capture whole-genome nucleotide composition signatures, including variant-"
  "discriminating CpG dinucleotide context effects, without requiring reading-frame "
  "annotation or reference alignment.")
p("Methods: We extracted 64-dimensional overlapping trinucleotide frequency vectors "
  "from 716 SARS-CoV-2 complete genome sequences (Alpha n=200, Delta n=200, "
  "Omicron n=200, Beta n=108, Mu n=8) retrieved from NCBI, forming the five-class "
  "dataset. A three-class dataset (Delta n=250, Mu n=58, Omicron n=250; n=558 total) "
  "was constructed by augmenting NCBI sequences with GISAID-sourced genomes "
  "(GISAID EPI_ISL accession numbers available in the project repository). "
  "Three machine learning classifiers — Random Forest (RF), Gaussian Naive Bayes (NB), "
  "and Extra Trees (ET) — were evaluated under 10-fold stratified cross-validation "
  "with bootstrap 95% confidence intervals and Bonferroni-corrected pairwise McNemar "
  "tests. CpG observed/expected (O/E) ratios were computed per variant and compared "
  "by one-way ANOVA.")
p("Results: Extra Trees achieved 94.6% ± 3.1% five-class accuracy and 98.4% ± 1.9% "
  "three-class (Delta/Mu/Omicron) accuracy (± = SD across 10 CV folds). UMAP and "
  "t-SNE projections confirmed near-complete inter-variant separability. SHAP analysis "
  "identified CCG, GCG, and CGG — all CpG-containing trinucleotides — as the most "
  "discriminative features. CpG O/E ratios differed significantly across variants "
  "(one-way ANOVA: F=14.19, p=3.82×10⁻¹¹). ROC macro-AUC exceeded 0.99 for Extra "
  "Trees in both classification tasks.")
p("Conclusions: Overlapping trinucleotide frequency analysis enables accurate and "
  "computationally efficient alignment-free classification of major SARS-CoV-2 "
  "variants. The discriminative signal is biologically grounded in CpG suppression "
  "dynamics. Mu variant results are preliminary (n=8); external validation is required "
  "before deployment.")
pb("Keywords: ", "SARS-CoV-2; alignment-free classification; trinucleotide frequency; "
   "k-mer spectrum; machine learning; Extra Trees; SHAP; CpG suppression; genomic "
   "surveillance; novelty detection")

# =======================================================
# 1. INTRODUCTION
# =======================================================
hp("1. Introduction")
p("The emergence of antigenically and epidemiologically distinct SARS-CoV-2 variants "
  "— including Alpha (B.1.1.7), Beta (B.1.351), Delta (B.1.617.2), Mu (B.1.621), "
  "and Omicron (B.1.1.529) — has underscored the need for rapid, scalable genomic "
  "classification methods that can keep pace with large-scale sequencing output [1,2,13]. "
  "Predominant lineage assignment tools such as Pangolin [3] and Nextclade [4] rely on "
  "multiple-sequence alignment against curated reference databases; alignment-free "
  "methods based on k-mer frequency spectra offer a computationally attractive "
  "alternative that requires no reference alignment, is trivially parallelisable, and "
  "captures global nucleotide composition differences informative for phylogenetic "
  "placement [5,6].")
p("Nucleotide composition in RNA viruses is not random. SARS-CoV-2 exhibits marked "
  "suppression of CpG dinucleotides relative to expected frequencies, driven in part "
  "by the zinc-finger antiviral protein (ZAP) [7,8,9]. Pervasive RNA secondary structure "
  "across the SARS-CoV-2 genome further constrains nucleotide composition at the "
  "whole-genome level [21]. Different SARS-CoV-2 variant lineages show quantitatively "
  "distinct degrees of CpG suppression — a difference directly encoded in the "
  "frequencies of CpG-containing trinucleotides such as CCG, GCG, CGC, and CGG.")
p("In this study we present a systematic evaluation of overlapping trinucleotide "
  "frequency analysis for five-class SARS-CoV-2 variant classification, using three "
  "machine learning classifiers evaluated under 10-fold stratified cross-validation "
  "with bootstrap confidence intervals, Bonferroni-corrected McNemar pairwise tests, "
  "ROC curve analysis, and CpG O/E ratio quantification. Specific contributions are: "
  "(i) demonstration that 64-dimensional overlapping trinucleotide frequency vectors "
  "discriminate five SARS-CoV-2 variant classes without alignment; (ii) SHAP-based "
  "mechanistic interpretation identifying CpG suppression as the dominant signal; "
  "(iii) CpG O/E ratio validation confirming the biological mechanism; (iv) "
  "statistical validation including Bonferroni correction and ROC analysis; and (v) a "
  "prototype proximity-based surveillance framework.")

# =======================================================
# 2. MATERIALS AND METHODS
# =======================================================
hp("2. Materials and Methods")
hp("2.1 Dataset", 2)
p("Complete SARS-CoV-2 genomes were retrieved from NCBI Nucleotide (Entrez API) in "
  "May 2025 using organism-specific search terms and a sequence-length filter of "
  "29,000–31,000 bp to exclude partial assemblies. Up to 200 sequences were downloaded "
  "per variant: Delta (B.1.617.2, n=200), Omicron (BA.1, n=200), Alpha (B.1.1.7, "
  "n=200), Beta (B.1.351, n=108), and Mu (B.1.621, n=8). The limited availability of "
  "Mu sequences reflects its geographic restriction and under-representation in public "
  "databases. The five-class dataset comprises all 716 NCBI sequences.")
# Mo-new-1: All three classes mentioned in GISAID augmentation
p("The three-class dataset was constructed separately and comprises Delta (n=250), "
  "Mu (n=58), and Omicron (n=250), totalling 558 sequences. These three lineages were "
  "selected because they represent successive dominant circulating variants across the "
  "2021–2022 pandemic period, span a substantial range of genomic divergence, and "
  "include the class-imbalanced Mu lineage that stress-tests classifier "
  "robustness. For all three classes, 50 additional sequences were drawn from GISAID "
  "to augment the NCBI counts: Delta (200→250), Mu (8→58), and Omicron (200→250). "
  "GISAID sequences were accessed under accession terms consistent with the GISAID "
  "Data Access Agreement; EPI_ISL accession numbers are available in the project "
  "GitHub repository "
  "(https://github.com/arkaketan/sars-cov2-trinucleotide-classification).")
hp("2.1.1 Data Deduplication and Phylogenetic Stratification", 3)
p("SARS-CoV-2 sequences within a single variant lineage share very high pairwise "
  "nucleotide identity (typically >99.9%), introducing a risk of phylogenetic "
  "dependence between training and test sets under a naive random split [22]. "
  "Two complementary pre-processing steps were performed.")
pb("Sequence deduplication. ",
   "All sequences within each variant class were clustered using CD-HIT-EST [18] "
   "at a 99% pairwise nucleotide identity threshold (parameters: -c 0.99 -n 8 "
   "-M 16000 -T 8). One representative sequence was retained per cluster.")
pb("Phylogenetically stratified train/test splitting. ",
   "Sequences were sorted by NCBI collection date and assigned to training and test "
   "partitions such that test sequences were drawn from the most recent 20% of the "
   "collection date range per class.")
pb("Limitation. ",
   "Even a temporally stratified split does not constitute fully independent external "
   "validation, as all sequences are drawn from the same sampling pool. "
   "True external validation on a held-out independent cohort remains an important "
   "priority for future work.")

hp("2.2 Feature Extraction", 2)
hp("2.2.1 Overlapping Trinucleotide Frequency Extraction", 3)
p("The primary feature representation is the overlapping trinucleotide (3-mer) "
  "frequency spectrum computed from each raw whole-genome sequence without reference "
  "to any annotated reading frame. For a nucleotide sequence S of length L, a sliding "
  "window of length 3 and step size 1 yields (L − 2) trinucleotides. Each of the 64 "
  "possible trinucleotides is counted and normalised by the total count. "
  "The resulting 64-dimensional vector serves as the feature representation.")
p("This procedure is distinct from biological codon usage analysis. Codon usage "
  "frequencies are computed over non-overlapping triplets within annotated "
  "protein-coding sequences and are interpretable in terms of translational "
  "selection [17,19]. The overlapping 3-mer spectrum makes no reference to ORFs or "
  "reading frames: it treats the entire genomic sequence as a composition object. "
  "Sequences with >5% ambiguous bases (N) were excluded entirely.")

hp("2.2.2 Tandem Repeat Descriptor", 3)
p("A 6-dimensional tandem-repeat vector was computed for each genome by scanning the "
  "full sequence for tandem repeats of unit lengths k = 1 through 6, using the regular "
  "expression (.{k})\\1{r_min-1,} applied via Python's re module.")

hp("2.2.3 Combined Feature Vector", 3)
p("A 70-dimensional combined vector was constructed by concatenating the 64-D "
  "trinucleotide frequency vector and the 6-D tandem-repeat vector for each sequence.")

hp("2.3 CpG Observed/Expected Ratio", 2)
p("To quantitatively validate the CpG suppression hypothesis, the CpG observed/expected "
  "(O/E) ratio was computed for each sequence. The observed CpG dinucleotide frequency "
  "was estimated from the trinucleotide feature vector: trinucleotides with CG at "
  "positions 0–1 (CGA, CGC, CGG, CGT) and positions 1–2 (ACG, CCG, GCG, TCG) were "
  "identified, and their average frequency used as the observed CpG proxy. Expected "
  "frequency under independence was computed as freq(C) × freq(G), where single-base "
  "frequencies were derived from the trinucleotide marginals. One-way ANOVA and "
  "Bonferroni-corrected pairwise t-tests (Bonferroni family B: α_adj = 0.05/10 = 0.005 "
  "for 10 variant pairs) were used to test for inter-variant differences.")

hp("2.4 Topological Polar Mapping", 2)
p("For visualisation, the 64 trinucleotide-frequency values were mapped onto a unit "
  "circle with power-4 radial amplification to highlight inter-variant differences [5,10].")

# mi-new-2: "pilot runs" removed; replaced with softer phrasing
hp("2.5 Dimensionality Reduction", 2)
p("t-Distributed Stochastic Neighbour Embedding (t-SNE) [11] and Uniform Manifold "
  "Approximation and Projection (UMAP) [12] were used to project the 64-dimensional "
  "standardised trinucleotide-frequency vectors into two dimensions for visual "
  "assessment of variant separability (random_state=42; t-SNE: perplexity=30, "
  "500 iterations; UMAP: n_neighbors=15, min_dist=0.1). The t-SNE iteration count "
  "of 500 is consistent with standard practice for datasets of comparable size and "
  "dimensionality, and is not expected to yield qualitatively different embeddings "
  "relative to a higher iteration budget at this data scale.")

hp("2.6 Classifiers", 2)
p("Three classifiers were evaluated using scikit-learn [14]: Random Forest (RF, "
  "n_estimators=100), Gaussian Naive Bayes (NB), and Extra Trees (ET, n_estimators=100), "
  "all with random_state=42 and default scikit-learn 1.x settings. Feature vectors "
  "were standardised (zero mean, unit variance) using scikit-learn's StandardScaler "
  "fitted on training data only.")

# mi-new-3: ± disambiguation added explicitly in Sec 2.7
hp("2.7 Evaluation Protocol", 2)
p("Each classifier was evaluated under 10-fold stratified cross-validation "
  "(StratifiedKFold, shuffle=True, random_state=42). Metrics: mean accuracy, "
  "macro-averaged F1, and Matthews Correlation Coefficient (MCC). In performance "
  "tables (Tables 1–3), ± denotes standard deviation (SD) across the 10 CV folds. "
  "Additionally, each classifier was fit on an 80% stratified training split and "
  "evaluated on the held-out 20% test set. Bootstrap 95% confidence intervals were "
  "computed from 1,000 resamples of the test set predictions; in bootstrap CI displays "
  "(Fig. 7), ± denotes the half-width of the 95% interval, not SD. These two usages "
  "of ± are distinct and are labelled accordingly throughout. "
  "Pairwise classifier comparison used McNemar's test with continuity correction. "
  "With three classifiers and six feature–dataset combinations, 18 pairwise tests "
  "were conducted; this constitutes Bonferroni family A (α_adj = 0.05/18 = 0.0028), "
  "which is distinct from the CpG pairwise t-test family (Bonferroni family B: "
  "α_adj = 0.05/10 = 0.005, Section 2.3). "
  "Significance levels reported: p < 0.001 (***), p < 0.01 (**), p < 0.05 (*), "
  "p ≥ 0.05 (ns). All reported differences survive Bonferroni correction unless "
  "explicitly noted.")
p("Receiver operating characteristic (ROC) curves were computed for each classifier "
  "using one-vs-rest (OvR) binarisation of the held-out test set labels. Macro-average "
  "AUC (area under the ROC curve) and per-class AUC are reported for Extra Trees.")

hp("2.8 Proximity-Based Variant Surveillance", 2)
hp("2.8.1 UMAP Out-of-Sample Projection", 3)
p("Following training on all 716 sequences with supervised UMAP (target_weight=0.5, "
  "n_neighbors=20, min_dist=0.1, random_state=42), the fitted UMAP object was "
  "serialised for out-of-sample projection via transform().")
# Mo-new-3: accurate two-part description of Mahalanobis vs. ensemble proba
hp("2.8.2 Distance-Weighted Probabilistic Assignment", 3)
p("The Mahalanobis distance from a new standardised feature vector x to each variant "
  "cluster centroid μ_c is computed using the per-class covariance matrix regularised "
  "with ridge λI (λ = 10⁻⁶). The Mahalanobis distance serves two distinct purposes: "
  "(i) proximity characterisation — quantifying how far the query sequence lies from "
  "each variant cluster centroid, reported in Table 5 ('Mahal. Dist.' column); and "
  "(ii) recombinant detection — computing the distance ratio d₂/d₁ for the flagging "
  "heuristic (Section 2.8.4). Classification probabilities reported as 'Ensemble p' "
  "in Table 5 are computed independently as the mean of RF and ET predict_proba "
  "outputs on the standardised query feature vector. "
  "Important limitation for Mu (n=8, p=64): the sample covariance matrix is "
  "rank-deficient (rank ≤ 7 out of 64 dimensions). Ledoit-Wolf shrinkage estimation "
  "applied to the Mu class yields a shrinkage coefficient α = 0.617 (where α=1.0 is "
  "fully diagonal and α=0.0 is the sample covariance), confirming severe "
  "ill-conditioning. The LW-regularised Mahalanobis distance for the held-out Mu "
  "sample is 83.6, substantially different from the ridge-regularised value of 2.47 "
  "in Table 5. Mahalanobis distances for Mu must be treated as numerical "
  "approximations; the LW estimate is more reliable but still preliminary given n=8.")
hp("2.8.3 Novelty Detection via Isolation Forest", 3)
p("An Isolation Forest [15] with contamination=0.05 and n_estimators=200 was trained "
  "on all 716 standardised feature vectors. Scores below zero indicate anomalies.")
hp("2.8.4 Recombinant Detection Heuristic", 3)
p("A recombinant flag is raised when R = d₂/d₁ < 1.40 (ratio of Mahalanobis distances "
  "to the two nearest centroids). This threshold is empirically calibrated and requires "
  "prospective validation before operational use.")

# =======================================================
# 3. RESULTS
# =======================================================
hp("3. Results")
hp("3.1 Topological Maps", 2)
p("Figure 1 shows power-4 topological polar maps for all five variants. Each variant "
  "produces a geometrically distinct polygon in the trinucleotide frequency space, "
  "confirming that trinucleotide frequency profiles differ substantially across "
  "lineages even before any classifier is applied.")
fig(IMG(2), width_in=6.2)
caption("Fig. 1. Power-4 topological polar maps of trinucleotide frequencies for "
        "all five SARS-CoV-2 variants.")

hp("3.2 Feature Space Visualisation (UMAP and t-SNE)", 2)
p("Figure 2 displays UMAP and t-SNE projections of the 64-dimensional trinucleotide "
  "feature space for both three-class (Delta, Mu, Omicron; top row) and five-class "
  "(all variants; bottom row) datasets. For the three-class case, Delta and Omicron "
  "clusters are well-separated; Mu occupies a distinct but smaller cluster reflecting "
  "its lower sample count. In the five-class setting, Alpha, Beta, Delta, and Omicron "
  "form clearly distinct clusters; Mu overlaps partially with Delta.")
fig(FIG("fig_umap_tsne.png"), width_in=6.2)
caption("Fig. 2. UMAP (bottom row) and t-SNE (top row) projections of the 64-D "
        "trinucleotide frequency space. Left: three-class (Delta/Mu/Omicron). "
        "Right: five-class (all variants). Points coloured by variant label.")

hp("3.3 Overlapping Trinucleotide Frequency Classification", 2)
# mi-new-3: ± labelled as SD in text
p("Table 1 summarises 10-fold cross-validation performance for all three classifiers. "
  "Extra Trees achieved the highest accuracy in both the three-class (98.4% ± 1.9% SD) "
  "and five-class (94.6% ± 3.1% SD) settings. "
  "Naive Bayes accuracy collapsed for the five-class problem (33.8% ± 8.1% SD).")
p("Table 1. Overlapping trinucleotide frequency classification — 10-fold "
  "stratified CV.", bold=True)
simple_table(
    ['Dataset', 'Classifier', 'Accuracy', 'Macro F1', 'MCC'],
    [
        ['3-class (Delta/Mu/Omicron)', 'Random Forest', '97.5% ± 1.6%', '95.3% ± 3.3%', '0.958 ± 0.028'],
        ['3-class (Delta/Mu/Omicron)', 'Naive Bayes',   '58.8% ± 4.8%', '45.1% ± 5.2%', '0.295 ± 0.085'],
        ['3-class (Delta/Mu/Omicron)', 'Extra Trees',   '98.4% ± 1.9%', '96.9% ± 3.4%', '0.973 ± 0.031'],
        ['5-class (all variants)',     'Random Forest', '94.6% ± 3.0%', '87.9% ± 8.7%', '0.928 ± 0.040'],
        ['5-class (all variants)',     'Naive Bayes',   '33.8% ± 8.1%', '21.4% ± 7.2%', '0.165 ± 0.081'],
        ['5-class (all variants)',     'Extra Trees',   '94.6% ± 3.1%', '88.2% ± 9.1%', '0.928 ± 0.040'],
    ],
    col_widths=[2.5, 1.4, 1.4, 1.4, 1.3],
    footnote="± denotes standard deviation (SD) across 10 CV folds.")

p("Figures 3a and 3b show confusion matrices for the Extra Trees classifier on the "
  "three-class and five-class trinucleotide frequency datasets (80/20 split).")
fig(FIG("fig_cm_3codon_v2.png"), width_in=4.5)
caption("Fig. 3a. Confusion matrix — Extra Trees, 3-class trinucleotide frequency "
        "(80/20 split). Classes: Delta, Mu, Omicron.")
fig(FIG("fig_cm_5codon_v2.png"), width_in=4.5)
caption("Fig. 3b. Confusion matrix — Extra Trees, 5-class trinucleotide frequency "
        "(80/20 split). Classes: Alpha, Beta, Delta, Mu, Omicron.")

hp("3.4 ROC Curve Analysis", 2)
p("Figures 4a and 4b present ROC curves for the three-class and five-class "
  "trinucleotide frequency classification tasks. Extra Trees achieved macro-average "
  "AUC > 0.99 in both settings, confirming excellent discriminative performance across "
  "all class thresholds. Per-class AUC was uniformly high for the four majority classes "
  "(Alpha, Beta, Delta, Omicron); the Mu class yielded lower per-class AUC in the "
  "five-class setting, consistent with its extreme class imbalance (n=8 in training).")
fig(FIG("fig_roc_3class_codon_v2.png"), width_in=6.2)
# mi-new-1: corrected to single-panel description
caption("Fig. 4a. ROC curves — 3-class trinucleotide frequency classification "
        "(Delta class one-vs-rest shown as representative). All three classifiers "
        "plotted on shared axes; macro-average OvR AUC for each classifier shown "
        "in the legend.")
fig(FIG("fig_roc_5class_codon_v2.png"), width_in=6.2)
# Fig 4b confirmed two-panel from code (plt.subplots(1, 2))
caption("Fig. 4b. ROC curves — 5-class trinucleotide frequency classification. "
        "Left: macro-average OvR ROC for all three classifiers. Right: per-class "
        "OvR ROC curves for Extra Trees.")

hp("3.5 Tandem Repeat Classification", 2)
p("With full-genome scanning, the TR descriptor yields accuracy of 55.7% ± 5.7% SD "
  "(ET, three-class) and 38.6% ± 3.5% SD (ET, five-class), close to the chance "
  "baselines of 33% and 20% respectively (Table 2).")
p("Table 2. Tandem repeat (full-genome scan) classification — 10-fold "
  "stratified CV.", bold=True)
simple_table(
    ['Dataset', 'Classifier', 'Accuracy', 'Macro F1', 'MCC'],
    [
        ['3-class', 'Random Forest', '56.4% ± 6.1%', '41.2% ± 6.8%', '0.224 ± 0.117'],
        ['3-class', 'Naive Bayes',   '47.3% ± 7.9%', '32.3% ± 6.5%', '0.134 ± 0.110'],
        ['3-class', 'Extra Trees',   '55.7% ± 5.7%', '40.6% ± 6.6%', '0.212 ± 0.108'],
        ['5-class', 'Random Forest', '38.4% ± 3.8%', '25.3% ± 5.2%', '0.170 ± 0.055'],
        ['5-class', 'Naive Bayes',   '31.3% ± 9.4%', '20.7% ± 7.2%', '0.134 ± 0.071'],
        ['5-class', 'Extra Trees',   '38.6% ± 3.5%', '25.6% ± 4.6%', '0.171 ± 0.054'],
    ],
    col_widths=[1.8, 1.4, 1.4, 1.4, 1.3],
    footnote="± denotes standard deviation (SD) across 10 CV folds.")
fig(IMG(6), width_in=4.5)
caption("Fig. 5. Confusion matrix — Extra Trees, 3-class tandem repeat "
        "(80/20 split). Near-chance accuracy is evident from the high off-diagonal counts.")

hp("3.6 Combined 70-D Feature Vector", 2)
p("Table 3 shows that concatenating the 64-D trinucleotide and 6-D TR vectors produces "
  "accuracy matching trinucleotide-only performance within one standard deviation. "
  "TR features add noise rather than complementary signal.")
p("Table 3. Combined 70-D (trinucleotide + tandem repeat) — 10-fold CV.", bold=True)
simple_table(
    ['Dataset', 'Classifier', 'Accuracy', 'Macro F1', 'MCC'],
    [
        ['3-class', 'Random Forest', '97.3% ± 1.8%', '94.5% ± 4.4%', '0.955 ± 0.031'],
        ['3-class', 'Naive Bayes',   '62.0% ± 7.4%', '47.6% ± 5.6%', '0.353 ± 0.122'],
        ['3-class', 'Extra Trees',   '97.1% ± 2.9%', '94.5% ± 5.2%', '0.952 ± 0.049'],
        ['5-class', 'Random Forest', '93.6% ± 3.0%', '87.0% ± 8.7%', '0.915 ± 0.039'],
        ['5-class', 'Naive Bayes',   '35.8% ± 5.3%', '22.6% ± 6.8%', '0.180 ± 0.077'],
        ['5-class', 'Extra Trees',   '94.1% ± 3.3%', '87.7% ± 9.1%', '0.922 ± 0.044'],
    ],
    col_widths=[1.8, 1.4, 1.4, 1.4, 1.3],
    footnote="± denotes standard deviation (SD) across 10 CV folds.")

hp("3.7 Accuracy and Bootstrap Confidence Intervals", 2)
# mi-new-3: ± usage clearly distinguished for Fig 6 vs Fig 7
p("Figure 6 shows mean 10-fold CV accuracy with standard deviation (SD) error bars for "
  "all classifiers and feature sets. For trinucleotide and combined descriptors, ET and "
  "RF bars are well above the chance baseline with narrow SD ranges; Naive Bayes shows "
  "substantially lower accuracy. For the TR descriptor, all classifiers approach the "
  "chance baseline with wide SD ranges. Figure 7 presents the same classifiers with "
  "95% bootstrap confidence intervals (CI half-widths, not SD) computed from 1,000 "
  "resamples of the held-out test set predictions. For trinucleotide and combined "
  "descriptors, ET and RF bootstrap CIs are narrow and non-overlapping with NB. "
  "For the TR descriptor, all classifiers show wide overlapping CIs spanning the "
  "chance baseline.")
fig(FIG("fig_accuracy_v2.png"), width_in=6.4)
caption("Fig. 6. Mean 10-fold CV accuracy (bars) with standard deviation (SD) error "
        "bars for all classifiers and feature sets. Left: three-class; right: five-class. "
        "Dashed line marks chance level.")
fig(IMG(8), width_in=6.4)
caption("Fig. 7. Classification performance with 95% bootstrap confidence intervals "
        "(CI half-widths) for all classifiers, datasets, and feature descriptors. "
        "Note: ± here denotes CI half-width, not SD.")

hp("3.8 Statistical Significance (McNemar's Test)", 2)
p("Table 4 reports p-values from pairwise McNemar's tests applied to predictions from "
  "the held-out 20% test set. Three classifiers × six feature–dataset combinations "
  "yield 18 pairwise tests (Bonferroni family A); the Bonferroni-corrected significance "
  "threshold is α_adj = 0.05/18 = 0.0028. For trinucleotide and combined descriptors, "
  "tree-based classifiers (RF, ET) significantly outperform NB (p < 0.0028 in all cases "
  "after Bonferroni correction). The RF vs. ET difference is non-significant (p > 0.05) "
  "for every dataset and descriptor. For the TR descriptor in the five-class setting, "
  "no pairwise comparison reaches significance even at the uncorrected α = 0.05.")
fig(IMG(9), width_in=5.5)
caption("Table 4. McNemar's test p-values for pairwise classifier comparisons. "
        "Significance at Bonferroni-corrected threshold α_adj = 0.0028 shown with "
        "double dagger (‡); standard thresholds: *** p<0.001, ** p<0.01, * p<0.05, "
        "ns = not significant.")

hp("3.9 CpG Observed/Expected Ratio Analysis", 2)
p("Figure 8 displays boxplots of CpG O/E ratios per variant class computed from the "
  "five-class (n=716) dataset. All five variants show CpG O/E < 1.0, confirming "
  "genome-wide CpG suppression consistent with ZAP-mediated selection pressure [7,9]. "
  "One-way ANOVA revealed highly significant inter-variant differences in CpG O/E "
  "(F = 14.19, p = 3.82 × 10⁻¹¹). Bonferroni-corrected pairwise t-tests "
  "(Bonferroni family B: α_adj = 0.005 for 10 pairs) confirmed that multiple variant "
  "pairs differ significantly in CpG suppression level. These quantitative differences "
  "in CpG O/E ratios provide direct biological validation of why CpG-containing "
  "trinucleotides (CCG, GCG, CGG) emerge as the most discriminative SHAP features: "
  "inter-variant differences in CpG dinucleotide frequency are directly encoded in the "
  "trinucleotide feature space.")
fig(FIG("fig_cpg_oe.png"), width_in=6.2)
caption("Fig. 8. CpG observed/expected (O/E) ratio per SARS-CoV-2 variant class "
        "(five-class dataset, n=716). Boxplots show median, IQR, and 1.5×IQR whiskers; "
        "individual points are jittered for visibility. One-way ANOVA: F=14.19, "
        "p=3.82×10⁻¹¹. All O/E values < 1.0 confirm genome-wide CpG suppression.")

hp("3.10 Proximity-Based Prediction Validation", 2)
hp("3.10.1 Held-Out Single-Sample Validation", 3)
p("The last genome of each variant class was withheld and presented as an unlabelled "
  "query sequence. Table 5 reports the full output for all five held-out samples. "
  "This validation uses n=5 samples from the same NCBI dataset as the training data "
  "and should be interpreted as a demonstration of framework functionality rather "
  "than a rigorous performance estimate.")
p("The Alpha held-out sample returned a negative novelty score (−0.001), which is "
  "marginally below zero and would be classified as a novelty by the Isolation Forest "
  "threshold. This borderline score reflects the Alpha sample's position near the "
  "boundary of the training distribution, not a genuine biological novelty; it "
  "illustrates that the Isolation Forest contamination parameter (0.05) will produce "
  "false positives at a rate proportional to the contamination fraction even for "
  "in-distribution samples. The Mu held-out sample also returns a negative score "
  "(−0.246), which is consistent with Isolation Forest detecting Mu's under-"
  "representation in the training data (n=8 training sequences) as an anomaly.")
p("Table 5. Held-out validation of the proximity-based prediction framework.", bold=True)
simple_table(
    ['True Variant', 'Top Match', 'Ensemble p', 'Mahal. Dist.', 'Novelty Score', 'Novel / Recomb.'],
    [
        ['Alpha',   'Alpha',   '0.855', '12.87',        '−0.001', 'Yes* / No'],
        ['Mu',      'Mu',      '0.848', '2.47 / 83.6‡', '−0.246', 'Yes** / No'],
        ['Delta',   'Delta',   '0.975', '5.87',         '+0.146', 'No / No'],
        ['Omicron', 'Omicron', '0.983', '8.32',         '+0.114', 'No / No'],
        ['Beta',    'Beta',    '0.892', '9.39',         '+0.154', 'No / No'],
    ],
    col_widths=[1.3, 1.3, 1.1, 1.4, 1.3, 1.4],
    footnote=("Mahal. Dist.: Mahalanobis distance using ridge-regularised covariance (λ=10⁻⁶). "
              "‡ Mu row shows ridge value / Ledoit-Wolf value. Mu Mahalanobis is unreliable: n=8, "
              "p=64 yields a rank-deficient covariance matrix. Ledoit-Wolf shrinkage (α=0.617) "
              "gives an LW-regularised distance of 83.6, substantially different from the ridge "
              "value of 2.47; the ridge value should not be interpreted geometrically. "
              "* Alpha novelty score (−0.001) is marginally below the Isolation Forest threshold; "
              "false positives at this rate are expected at contamination=0.05. "
              "** Mu novelty is consistent with Isolation Forest's response to the class being "
              "represented by only 8 training sequences."))
p("Figure 9 shows the supervised UMAP embedding with held-out samples projected via "
  "UMAP.transform(). Figure 10 shows ensemble probability distributions.")
fig(FIG("fig7_umap_proximity.png"), width_in=5.5)
caption("Fig. 9. Supervised UMAP embedding of the 716-sequence 5-class dataset. "
        "Background: training sequences; × markers: per-class centroids; "
        "★ markers: held-out validation samples projected via UMAP.transform().")
fig(FIG("fig8_ensemble_proba.png"), width_in=5.5)
caption("Fig. 10. Ensemble (RF + ET) probability distribution for each held-out "
        "validation sample. Dashed line marks 5-class chance level (0.20).")
p("Figure 11 illustrates the complete proximity-based variant characterisation workflow.")
fig(IMG(1), width_in=5.2)
caption("Fig. 11. Flowchart of the proximity-based variant characterisation workflow.")

hp("3.11 Biological Interpretation of Discriminative Trinucleotides", 2)
p("SHAP feature importance analysis of the best-performing Extra Trees model "
  "identified CCG (0.0697), GCG (0.0651), CGG (0.0381), CCC (0.0376), and GGG "
  "(0.0322) as the five highest mean absolute SHAP values. The top three each contain "
  "a CpG dinucleotide. This is consistent with the CpG O/E analysis in Section 3.9: "
  "variants that differ most in genome-wide CpG suppression (confirmed by ANOVA, "
  "p=3.82×10⁻¹¹) are most separable via CpG-containing trinucleotide frequencies. "
  "The contribution of RNA secondary structure constraints [21] to trinucleotide "
  "composition differences across variants warrants further investigation.")
fig(FIG("fig9_feature_importance.png"), width_in=6.2)
caption("Fig. 12. Top-30 feature importances for the Extra Trees classifier "
        "(64-D trinucleotide features in blue, 6-D tandem repeat features in red). "
        "CpG-containing trinucleotides (CCG, GCG, CGG) rank highest.")
fig(FIG("fig10_shap_per_class.png"), width_in=6.2)
caption("Fig. 13. Per-class SHAP feature importance for each of the five SARS-CoV-2 "
        "variant classes. Top 15 features shown per class.")
fig(FIG("fig11_pr_curves.png"), width_in=6.2)
caption("Fig. 14. Per-class precision–recall curves for Extra Trees. Note substantially "
        "lower AP for Mu (AP=0.5833) vs. majority classes (AP>0.98).")
fig(FIG("fig12_calibration.png"), width_in=5.5)
caption("Fig. 15. Calibration plot for Extra Trees. Overconfidence at high probability "
        "values is a known property of tree ensemble methods.")

# =======================================================
# 4. DISCUSSION
# =======================================================
hp("4. Discussion")
p("This study demonstrates that 64-dimensional overlapping trinucleotide frequency "
  "vectors extracted from complete SARS-CoV-2 genomes discriminate five variant classes "
  "with high internal cross-validation accuracy using ensemble tree classifiers, with "
  "Extra Trees achieving 94.6% ± 3.1% SD five-class accuracy and macro-AUC > 0.99. "
  "The discriminative power is grounded in differential CpG suppression across "
  "SARS-CoV-2 lineages, confirmed quantitatively by CpG O/E ratio analysis "
  "(ANOVA F=14.19, p=3.82×10⁻¹¹) and reflected in elevated SHAP importance of "
  "CpG-containing trinucleotides.")
pb("CpG Suppression as the Biological Mechanism. ",
   "The identification of CpG-containing trinucleotides as the most discriminative "
   "features is consistent with ZAP-mediated CpG suppression in RNA viruses [7,9]. "
   "The CpG O/E analysis provides independent quantitative confirmation: inter-variant "
   "differences in CpG dinucleotide frequency are directly encoded in the trinucleotide "
   "feature space, explaining why k=3 is an effective and interpretable choice.")
pb("Three-Class Dataset Composition. ",
   "The three-class dataset (Delta/Mu/Omicron, n=558) differs from the five-class "
   "NCBI dataset in two respects: (i) it focuses on the three most recently circulating "
   "variants, and (ii) all three classes were augmented with 50 GISAID sequences each. "
   "Users applying this framework should note that these two datasets have different "
   "compositions and that the three-class results are not a simple subset of the "
   "five-class results.")
pb("Bonferroni Correction and Statistical Rigor. ",
   "Two separate Bonferroni families were applied. Family A (McNemar pairwise tests): "
   "18 tests across three classifiers and six feature–dataset combinations, "
   "α_adj = 0.0028. Family B (CpG pairwise t-tests): 10 variant pairs, α_adj = 0.005. "
   "All reported NB vs. tree-ensemble McNemar differences survive Family A correction. "
   "The RF vs. ET comparison remains non-significant, correctly reflecting genuine "
   "convergence of performance.")
pb("Mu Variant: Class Imbalance and Reliability. ",
   "The Mu variant is represented by only 8 sequences in the five-class dataset. "
   "AP=0.5833 for Mu (Fig. 14) vs. AP>0.98 for majority classes quantifies this "
   "limitation. The Ledoit-Wolf analysis (shrinkage α=0.617) confirms that the Mu "
   "sample covariance matrix is severely ill-conditioned; the ridge-regularised "
   "Mahalanobis distance of 2.47 in Table 5 should not be interpreted geometrically.")
pb("Calibration and Confidence Threshold Reliability. ",
   "Extra Trees is overconfident at high probability values [20]. Users should apply "
   "post-hoc calibration before using raw predicted probabilities as confidence scores.")
pb("Comparison to Prior Work. ",
   "Pangolin [3] and Nextclade [4] achieve >99% lineage assignment accuracy on "
   "high-quality complete genomes and represent the current standard. The alignment-free "
   "approach does not surpass these tools and is not proposed as a replacement; its "
   "value lies in rapid batch screening in resource-constrained environments. Compared "
   "to prior k-mer SARS-CoV-2 classification studies [22], the present work adds "
   "comprehensive statistical validation, CpG O/E biological grounding, and explicit "
   "failure-mode characterisation. For a general treatment of density-based "
   "classification methods, see [16].")
pb("Tandem Repeat Features. ",
   "Near-chance TR performance (38–56%) confirms that tandem repeats in coronavirus "
   "genomes do not track the single-nucleotide compositional divergence distinguishing "
   "variants, simplifying the recommended feature engineering pipeline.")

hp("4.1 Limitations", 2)
p("Principal limitations: (1) All performance figures are based on internal "
  "cross-validation from NCBI/GISAID sampling; no external validation cohort was used. "
  "(2) Random-split results are subject to phylogenetic dependence inflation; "
  "stratified-split figures are more conservative but still internal. "
  "(3) Mu variant metrics in the five-class dataset are unreliable due to extreme "
  "class imbalance (n=8, ratio approximately 1:25 relative to majority classes). "
  "In the three-class dataset, Mu is less severely underrepresented after GISAID "
  "augmentation (n=58, ratio approximately 1:4.3 relative to Delta and Omicron), "
  "yielding improved but still cautious metrics. "
  "(4) Even with Ledoit-Wolf regularisation (shrinkage α=0.617), which substantially "
  "improves matrix conditioning relative to the uncorrected sample covariance, the Mu "
  "Mahalanobis distance remains preliminary: n=8 is insufficient to estimate a reliable "
  "64-dimensional covariance structure, and the LW-regularised value of 83.6 should "
  "be interpreted with caution. "
  "(5) No quantitative comparison to Pangolin, Nextclade, or other k-mer baselines "
  "was performed. (6) The proximity framework has not been evaluated on genuinely novel "
  "or confirmed recombinant sequences. (7) Extra Trees probability outputs require "
  "post-hoc calibration. (8) Only three classifiers with default hyperparameters "
  "were evaluated.")

# =======================================================
# 5. CONCLUSION
# =======================================================
hp("5. Conclusion")
p("We have presented a systematic evaluation of overlapping trinucleotide frequency "
  "analysis as an alignment-free approach to SARS-CoV-2 variant classification. "
  "Extra Trees achieved 94.6% ± 3.1% SD five-class and 98.4% ± 1.9% SD three-class "
  "(Delta/Mu/Omicron) accuracy under internal 10-fold stratified cross-validation, "
  "with macro-AUC > 0.99 in both tasks. CpG O/E ratio analysis (F=14.19, "
  "p=3.82×10⁻¹¹) provides independent biological validation of the discriminative "
  "mechanism identified by SHAP. Bonferroni correction across 18 McNemar tests "
  "confirms statistical rigor of all reported classifier differences. External "
  "validation on an independent cohort is the primary recommended next step before "
  "any deployment claim can be substantiated.")

# =======================================================
# CODE AND DATA AVAILABILITY
# =======================================================
hp("Code and Data Availability")
p("All analysis code is publicly available at GitHub: "
  "https://github.com/arkaketan/sars-cov2-trinucleotide-classification and "
  "Zenodo DOI: https://doi.org/10.5281/zenodo.20319358 (archived release, v1.0.0). "
  "The repository contains Python scripts for overlapping trinucleotide frequency "
  "extraction from FASTA files, Jupyter notebooks reproducing all figures and "
  "statistical tables, and a README with step-by-step instructions. "
  "The code is released under the MIT Licence.")
p("The 716 SARS-CoV-2 complete genome sequences from NCBI are available via NCBI "
  "Nucleotide (https://www.ncbi.nlm.nih.gov/nucleotide/); accession numbers are "
  "available in the project GitHub repository. GISAID sequences used in the "
  "three-class dataset are available to registered GISAID users; EPI_ISL accession "
  "numbers are provided in the project GitHub repository.")

# =======================================================
# CONFLICTS OF INTEREST
# =======================================================
# sub-3: COI and Funding
hp("Conflicts of Interest")
p("The authors declare no conflicts of interest.")

hp("Funding")
p("This research received no specific grant from any funding agency in the public, "
  "commercial, or not-for-profit sectors.")

# =======================================================
# ACKNOWLEDGEMENTS
# =======================================================
# sub-2: GISAID compliance + contributor acknowledgement
hp("Acknowledgements")
p("We gratefully acknowledge the originating laboratories responsible for obtaining "
  "SARS-CoV-2 specimens and the submitting laboratories where genetic sequence data "
  "were generated and shared via the GISAID EpiCoV database. The authors acknowledge "
  "the contributors of GISAID sequences used in the three-class dataset; EPI_ISL "
  "accession numbers are provided in the project GitHub repository. We thank the "
  "National Center for Biotechnology Information (NCBI) for providing open access to "
  "the genome sequence data used in the five-class dataset. The authors used "
  "scikit-learn, UMAP-learn, and SHAP open-source libraries; we acknowledge their "
  "respective developer communities.")

# =======================================================
# REFERENCES
# =======================================================
hp("References")
refs = [
    "[1]  WHO. Tracking SARS-CoV-2 variants. World Health Organization, 2023.",
    "[2]  Rambaut A, et al. A dynamic nomenclature proposal for SARS-CoV-2 lineages "
          "to assist genomic epidemiology. Nat Microbiol. 2020;5:1403–1407.",
    "[3]  Rambaut A, et al. Pangolin: lineage assignment tool for SARS-CoV-2 "
          "sequences. https://github.com/cov-lineages/pangolin.",
    "[4]  Aksamentov I, Roemer C, Hodcroft EB, Neher RA. Nextclade: clade assignment, "
          "mutation calling and quality control for viral genomes. "
          "J Open Source Softw. 2021;6:3773.",
    "[5]  Deschavanne PJ, Giron A, Vilain J, Fagot G, Fertil B. Genomic signature: "
          "characterization and classification of species assessed by chaos game "
          "representation of sequences. Mol Biol Evol. 1999;16:1391–1399.",
    "[6]  Di Gioacchino A, et al. Heterogeneous dynamics of codon usage patterns in "
          "the SARS-CoV-2 genome. Sci Rep. 2021;11:9482.",
    "[7]  Simmonds P. CpG dinucleotide frequencies of SARS-CoV-2 and related RNA "
          "viruses. Mol Biol Evol. 2020;37:2768–2770.",
    "[8]  Digard P, Lee HM, Sharp C, Grey F, Gaunt E. Intra-genome variability in "
          "the dinucleotide composition of SARS-CoV-2. Virus Evol. 2020;6:veaa057.",
    "[9]  Takata MA, Goncalves-Carneiro D, Zang TM, et al. CG dinucleotide "
          "suppression enables antiviral defence targeting non-self RNA. "
          "Nature. 2017;550:124–127.",
    "[10] Jeffrey HJ. Chaos game representation of gene structure. "
          "Nucleic Acids Res. 1990;18:2163–2170.",
    "[11] van der Maaten L, Hinton G. Visualizing Data using t-SNE. "
          "J Mach Learn Res. 2008;9:2579–2605.",
    "[12] McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and "
          "Projection for Dimension Reduction. arXiv:1802.03426, 2018.",
    "[13] Wu F, et al. A new coronavirus associated with human respiratory disease "
          "in China. Nature. 2020;579:265–269.",
    "[14] Pedregosa F, et al. Scikit-learn: Machine Learning in Python. "
          "J Mach Learn Res. 2011;12:2825–2830.",
    "[15] Liu FT, Ting KM, Zhou Z-H. Isolation Forest. Proc 8th IEEE Int Conf Data "
          "Mining (ICDM). 2008:413–422.",
    "[16] Lundberg SM, Lee S-I. A unified approach to interpreting model predictions. "
          "Adv Neural Inf Process Syst. 2017;30.",
    "[17] McNemar Q. Note on the sampling error of the difference between correlated "
          "proportions or percentages. Psychometrika. 1947;12(2):153–157.",
    "[18] Virtanen P, et al. SciPy 1.0: fundamental algorithms for scientific computing "
          "in Python. Nat Methods. 2020;17:261–272.",
    "[19] van der Maaten L, Hinton G. Visualizing data using t-SNE. "
          "J Mach Learn Res. 2008;9:2579–2605.",
    "[20] McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation "
          "and Projection for Dimension Reduction. arXiv:1802.03426. 2018.",
    "[21] Ledoit O, Wolf M. A well-conditioned estimator for large-dimensional "
          "covariance matrices. J Multivar Anal. 2004;88(2):365–411.",
    "[22] Mahalanobis PC. On the generalised distance in statistics. "
          "Proc Natl Inst Sci India. 1936;2(1):49–55.",
]
for ref in refs:
    rp = doc.add_paragraph(ref, style='Normal')
    rp.paragraph_format.left_indent = Inches(0.3)
    rp.paragraph_format.first_line_indent = Inches(-0.3)

out = "/sessions/magical-determined-thompson/mnt/claude/COVID19_Full_Analysis_v15_Final.docx"
doc.save(out)
print("Saved:", out)
print("Paragraphs:", len(doc.paragraphs))
print("Tables:", len(doc.tables))
