# Missense Variant Pathogenicity Prediction with ESM-2 and Gene-Held-Out Validation

Predicting whether a missense variant is likely disease-causing, using protein language model scores, evaluated on genes the model has never seen, and applied to prioritize Variants of Uncertain Significance (VUS) from ClinVar.

**Live demo:**
https://variant-pathogenicity-esm2-m7to79upw5p7sanscfhbgh.streamlit.app/


> **Research prototype. Not a clinical tool.** Predictions must not be used for diagnosis or treatment decisions.

---

## Why this problem matters

Genetic tests often return variants of uncertain significance (VUS): changes where nobody knows whether the variant causes disease. These results leave patients and clinicians without an answer. Computational evidence can help curators decide which VUS deserve follow-up first.

## What this project does

1. Builds a labeled dataset of missense variants from ClinVar (pathogenic vs. benign).
2. Scores every variant with **ESM-2** (35M parameters), a protein language model, using a zero-shot log-likelihood ratio: log P(mutant amino acid) minus log P(original amino acid) at the masked position.
3. Adds simple biochemical features: BLOSUM62, change in hydrophobicity, volume, and charge, and relative position in the protein.
4. Trains and compares models using **two validation designs**: a random split and a **gene-held-out split** (entire genes are kept out of training).
5. Checks calibration and per-gene performance, and ranks VUS by predicted pathogenicity.

## Data

| Source | Use |
|---|---|
| ClinVar `variant_summary.txt.gz` (NCBI FTP, downloaded November 2026) | Variants and clinical labels |
| UniProt (reviewed human proteins, REST API) | Protein sequences |
| ESM-2 `facebook/esm2_t12_35M_UR50D` (Hugging Face) | Language model scores |
| Biopython | BLOSUM62 matrix |

**Filtering steps (with counts):**

| Step | Rows |
|---|---|
| GRCh38 single-nucleotide variants with a protein change | 3,465,432 |
| True missense (20 standard amino acids, no stop codons) | 2,572,693 |
| Review status with assertion criteria or better | 2,472,860 |
| After removing duplicates | 2,338,514 |
| Genes with at least 10 pathogenic and 10 benign variants, top 80 by size | 80 genes |
| Genes dropped: no UniProt sequence found (10) or ClinVar/UniProt residue mismatch (8) | **62 genes** |
| **Final labeled variants** | **26,827** (14,055 pathogenic, 12,772 benign) |
| **VUS scored** | **94,048** |

Labels: Pathogenic and Likely pathogenic = 1; Benign and Likely benign = 0. Conflicting classifications were excluded.

**Reference check:** for every variant, the wild-type amino acid reported by ClinVar must match the UniProt sequence at that position. Genes where this failed (isoform or transcript mismatch) were removed.

## Methods

- **ESM-2 scoring:** the residue at each position is masked, and the model's probabilities for all 20 amino acids are recorded once per unique (gene, position) pair (77,985 positions) and reused across variants. Long proteins use a 1,000-residue window centered on the variant.
- **Models:** ESM score alone (no training), logistic regression, and LightGBM (200 trees, learning rate 0.05).
- **Validation:** 5-fold cross-validation, either random (stratified) or **grouped by gene** so test genes are never seen during training.
- **Metrics:** AUROC, AUPRC, MCC (mean over folds).

## Results

| Split | Model | AUROC | AUPRC | MCC |
|---|---|---|---|---|
| Random | ESM score only (no training) | 0.813 | 0.871 | 0.474 |
| Random | Simple features + LogReg | 0.733 | 0.749 | 0.316 |
| Random | Simple features + LightGBM | 0.814 | 0.832 | 0.496 |
| Random | ESM + simple + LightGBM | 0.883 | 0.913 | 0.622 |
| **Gene-held-out** | ESM score only (no training) | 0.808 | 0.856 | 0.460 |
| **Gene-held-out** | Simple features + LogReg | 0.731 | 0.733 | 0.317 |
| **Gene-held-out** | Simple features + LightGBM | 0.800 | 0.797 | 0.469 |
| **Gene-held-out** | **ESM + simple + LightGBM** | **0.868** (SD 0.044) | **0.893** | **0.582** |

**Key findings**

1. **The language model adds real value.** On held-out genes, adding ESM features raised AUROC from 0.800 to 0.868. The untrained ESM score alone (0.808) already matched the simple-feature LightGBM model.
2. **The gap between random and gene-held-out splits was small** (0.883 vs. 0.868). In this dataset, random splitting inflated performance only modestly, but gene-held-out validation is still the more honest estimate for unseen genes.
3. **Performance varies a lot by gene.** AUROC is about 0.99 for collagen genes (COL1A2, COL3A1, COL4A1, COL4A3, COL4A4) and 0.61 to 0.71 for KMT2C, SETBP1, RYR1, ANKRD11, TP53, KAT6A, DNAH11, and RYR2. A possible explanation is that pathogenicity in the harder genes depends on mechanisms or domains that a sequence-only model does not capture. This was not tested here.
4. **Calibration (gene-held-out predictions)** is close to the ideal diagonal, with slight underconfidence at low probabilities and slight overconfidence at higher ones.

Figures: `figures/fig_roc.png`, `figures/fig_calibration.png`, `figures/fig_shap.png`. Per-gene results: `results/per_gene_auroc.csv`.

## VUS prioritization

The final model was trained on all labeled variants and applied to 94,048 VUS: 5,042 scored above 0.9 and 3,025 below 0.1. The top-ranked VUS were mostly cysteine substitutions (for example to Y, F, or W) in genes such as FBN1, FBN2, NOTCH1, and COL3A1, which is consistent with the known structural importance of cysteines in these proteins. The top 1,000 are in `results/top1000_vus_predictions.csv`.

These are **computational prioritizations for follow-up, not classifications**. Many top scores are exactly 1.0, and VUS may differ in distribution from the labeled set, so the calibration seen on labeled variants may not transfer to them.

## Limitations

- ClinVar labels contain errors, and submitters sometimes disagree.
- Some ClinVar labels were influenced by older prediction tools, so there is partial **circularity**.
- Only genes with at least 10 variants in each class were used, so results may not generalize to rare genes. Gene selection itself may introduce bias.
- The labeled set is nearly balanced, which is unlike real clinical data where benign variants dominate. AUPRC and MCC values would be lower under realistic class ratios.
- Only missense variants are modeled. No structural features (for example AlphaFold) were used.
- Hyperparameters were not tuned, and only one ESM-2 size (35M) was tested.
- The 62 genes were chosen before splitting, using label counts from the full dataset.

## Future work

- Compare ESM-2 model sizes (8M, 150M).
- Add AlphaFold confidence (pLDDT) and domain annotations.
- Time-based validation: test on VUS that were later reclassified in newer ClinVar releases.
- Calibrate probabilities explicitly (for example isotonic regression) and evaluate under realistic class ratios.

## How to reproduce

1. Open the notebook in `notebooks/` in Google Colab with a T4 GPU.
2. Download `variant_summary.txt.gz` from `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/`.
3. Run the cells in order: filter ClinVar in chunks, extract protein changes, select genes, fetch UniProt sequences, validate reference residues, score positions with ESM-2, build features, then train and evaluate.
4. Install dependencies: `pip install -r requirements.txt`.

Large intermediate files (ClinVar, ESM position scores) are not included because of size.

## Repository structure

```
.
├── README.md
├── requirements.txt
├── notebooks/      # Colab notebook
├── results/        # metrics tables, per-gene AUROC, top VUS predictions
└── figures/        # ROC, calibration, SHAP plots
```

## Author

Tehreem Raouf · https://www.linkedin.com/in/tehreem-raouf
