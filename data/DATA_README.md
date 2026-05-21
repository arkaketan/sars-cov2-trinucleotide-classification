# Dataset Notes

## Included files

| File | Rows | Columns | Description |
|------|------|---------|-------------|
| `Combined_5class_Codon.csv` | 716 | 65 | 64-D overlapping trinucleotide frequencies + class label |
| `NCBI_TR_Fixed_Dataset.csv` | 716 | 7 | 6-D tandem-repeat counts + class label |
| `Combined_3class_Codon.csv` | 558 | 65 | 3-class subset (Alpha, Delta, Omicron) — trinucleotide |
| `Combined_3class_TR.csv` | 558 | 7 | 3-class subset — tandem-repeat |

## Class encoding

| Encode value | Variant | Lineage | Sequences |
|---|---|---|---|
| 1 | Alpha   | B.1.1.7   | 200 |
| 2 | Mu      | B.1.621   | 8   |
| 3 | Delta   | B.1.617.2 | 200 |
| 4 | Omicron | BA.1      | 200 |
| 5 | Beta    | B.1.351   | 108 |

## Feature extraction

Trinucleotide frequencies were computed using a sliding window of length 3 and step 1 over each complete genome sequence (no reading-frame annotation). Each of the 64 possible trinucleotides {A,C,G,T}^3 is counted and normalised by total trinucleotides (L−2). Column names are `c_0` through `c_63`, ordered AAA, AAC, AAG, ... TTT.

Tandem-repeat counts were extracted using IMEx (Mudunuri & Nagarajaram, 2007) scanning unit lengths 1–6 across the full ~30,000 bp genome.

## Provenance

All sequences were retrieved from NCBI Nucleotide (https://www.ncbi.nlm.nih.gov/nucleotide/) in May 2025. NCBI accession numbers are provided in Supplementary Table S1 of the manuscript. No sequences were obtained from GISAID.

## Reproducing from raw sequences

```python
from Bio import Entrez, SeqIO
Entrez.email = "your@email.com"

# Example: fetch Alpha genomes
handle = Entrez.esearch(db="nucleotide",
                        term="SARS-CoV-2[Organism] B.1.1.7[All Fields] complete genome",
                        retmax=200)
```

See `code/variant_proximity_predictor.py` (`extract_codon_features` and `extract_tr_features`) for the exact extraction logic used in this study.
