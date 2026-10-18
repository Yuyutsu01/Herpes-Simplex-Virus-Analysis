# Herpes Simplex Virus Proteome Analysis (HSV-1 & HSV-2)
### Curated Reference Engineering, Representation Learning, and Functional Classification Benchmark

[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyTorch CUDA](https://img.shields.io/badge/PyTorch-CUDA%2012.1-red.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioinformatics](https://img.shields.io/badge/Bioinformatics-Proteomics-green.svg)](https://github.com/)

---

## Overview

This repository provides an end-to-end, publication-grade bioinformatics and machine learning framework for the **Herpes Simplex Virus (HSV-1 and HSV-2) proteome**. Starting from 22,710 heterogeneous public database entries, the pipeline establishes a gold-standard 149-protein canonical reference set, propagates verified labels to 9,149 near-duplicate clinical variants (yielding a 9,298-protein corpus), extracts three distinct feature modalities (physicochemical, $k$-mer TF-IDF, and ProtBERT transformer embeddings), and conducts comprehensive within-species and cross-species benchmarking across 5 core viral functional ontologies:
- **Regulatory** (e.g., ICP0, ICP4, ICP27, ICP22)
- **Replication** (e.g., DNA polymerase, helicase-primase, thymidine kinase)
- **Capsid_Assembly** (e.g., VP5 major capsid, portal, terminase, scaffolding proteins)
- **Tegument** (e.g., VP16, VP22, large tegument UL36)
- **Glycoprotein_Envelope** (e.g., gB, gC, gD, gH, gL)

---

## Pipeline Architecture

```text
========================================================================================
                                 PIPELINE FLOWCHART
========================================================================================

  +-----------------------+     +-------------------------------+
  |  non_redundant.fasta  |     |  HSV2_non_redundant (2).fasta |
  |   (HSV-1, 11,811 seq) |     |      (HSV-2, 10,899 seq)      |
  +-----------+-----------+     +---------------+---------------+
              |                                 |
              +----------------+----------------+
                               |
                               v
               [STEP 1: CURATION & VERIFICATION]
               - run_full_curation.py -> hsv_curation_records.tsv (22,710 records)
               - step1_verify_lengths.py -> RefSeq validation (NC_001806.2, NC_001798.2)
               - step2_build_annotation.py -> 149 canonical KEEP reference proteins
               - step3_export_fasta.py -> hsv1_hsv2_curated_reference.fasta
                               |
                               v
                 [STEP 2: TRAINING CORPUS EXPANSION]
               - step1_merge_mapping.py -> pairwise alignment (>=90% id, <=10% len diff)
               - step2_build_training_set.py -> hsv_training_set.tsv (9,298 clean sequences)
               - step3_generate_splits.py -> group-aware 100-fold CV & transfer splits
                               |
                               v
                     [STEP 3: FEATURE EXTRACTION]
               - extract_physchem.py -> (9298 x 26) macro-descriptors + molar fractions
               - extract_kmer_tfidf.py -> (9298 x 128) 3/5-mer TF-IDF + TruncatedSVD
               - extract_protbert.py -> (9298 x 1024) GPU-accelerated ProtBERT embeddings
                               |
                               v
               [STEP 4: BENCHMARKING & ERROR ANALYSIS]
               - run_benchmarking.py -> Grouped 100-fold CV, transfer, ablation
               - run_visualizations_and_analysis.py -> UMAP/t-SNE, RF Gini, paired stats
                               |
                               v
                    [STEP 5: PUBLICATION FIGURES]
               - results/figures/fig1_confusion_matrix.png
               - results/figures/fig2_roc_pr_curves.png
               - results/figures/fig3_feature_comparison.png
               - results/figures/fig4_cross_species_transfer.png
               - results/figures/fig5_umap_tsne_protbert.png
               - results/figures/fig6_physchem_importance.png
========================================================================================
```

---

## Key Experimental Results

### 1. Within-Species Grouped Cross-Validation (100 Folds)
Cross-validation strictly holds out entire ortholog groups (zero homology leakage). ProtBERT representations dramatically outperform traditional sequence representations:

| Feature Modality | Optimal Model | Macro-F1 (Mean ± Std) | Balanced Accuracy | MCC |
|:---|:---|:---:|:---:|:---:|
| **ProtBERT Embeddings (1,024-d)** | **Logistic Regression** | **0.6331 ± 0.1376** | **0.6695 ± 0.1255** | **0.6341 ± 0.1499** |
| **Combined Representation (1,178-d)** | Logistic Regression | 0.6322 ± 0.1372 | 0.6685 ± 0.1250 | 0.6330 ± 0.1492 |
| **Physicochemical Descriptors (26-d)** | Random Forest | 0.3731 ± 0.1434 | 0.4129 ± 0.1456 | 0.2763 ± 0.1940 |
| **K-mer TF-IDF + SVD (128-d)** | Random Forest | 0.3134 ± 0.1432 | 0.3562 ± 0.1432 | 0.2073 ± 0.1877 |
| **Majority Baseline** | Dummy Classifier | 0.0935 ± 0.0000 | 0.2000 ± 0.0000 | 0.0000 ± 0.0000 |

* **Statistical Significance:**
  - ProtBERT vs. Physicochemical: $\Delta = +0.2600$, Paired $t = 14.10$ ($p = 2.03 \times 10^{-25}$), Wilcoxon $W = 77.0$ (**$p = 3.86 \times 10^{-17}$**).
  - ProtBERT vs. K-mer TF-IDF: $\Delta = +0.3197$, Paired $t = 20.68$ ($p = 1.02 \times 10^{-37}$), Wilcoxon $W = 6.0$ (**$p = 4.67 \times 10^{-18}$**).

### 2. Cross-Species Zero-Shot Transfer
Evaluating generalization across herpesvirus species boundaries:
- **HSV-1 $\to$ HSV-2 Transfer:** **0.8903 Macro-F1** (Balanced Accuracy: 0.8911, MCC: 0.8660)
- **HSV-2 $\to$ HSV-1 Transfer:** **0.9404 Macro-F1** (Balanced Accuracy: 0.9409, MCC: 0.9254)

---

## Repository Structure

```text
Herpes-Simplex-Virus-Analysis/
├── README.md                           # Main GitHub documentation
├── generate_curation_report.py         # Curation report compiler
├── data/                               # Raw FASTA sequences & reference genomes
│   ├── non_redundant.fasta             # Raw HSV-1 input sequences (11,811 headers)
│   ├── HSV2_non_redundant (2).fasta    # Raw HSV-2 input sequences (10,899 headers)
│   └── refseq_cache/                   # Cached NCBI RefSeq genomes (NC_001806.2, NC_001798.2)
├── curation/                           # Curation pipeline & canonical reference
│   ├── run_full_curation.py            # Master parsing, deduplication, & status labeling
│   ├── step1_verify_lengths.py         # RefSeq verification of ambiguous lengths
│   ├── verify_lengths.json             # Verification machine-readable output
│   ├── verify_lengths.md               # UL31 and UL27 verification narrative
│   ├── step2_build_annotation.py       # Canonical 149-protein annotation builder
│   ├── hsv_annotation.tsv              # 149 canonical reference annotations (12 columns)
│   ├── hsv_annotation_audit.md         # Audit breakdown of canonical proteome
│   ├── step3_export_fasta.py           # Standardized reference FASTA exporter
│   ├── hsv1_hsv2_curated_reference.fasta           # 149 canonical protein sequences
│   ├── hsv1_hsv2_curated_reference.manifest.tsv   # Cryptographic SHA-256 hashes
│   ├── hsv1_hsv2_curated_reference.excluded_duplicates.tsv
│   ├── hsv_curation_records.tsv        # Master 22,710-record curation table
│   ├── hsv_curation_table.md           # Master Markdown curation table
│   └── hsv_curation_report.md          # Curator summary report
├── training_set/                       # Expanded training corpus & split manifests
│   ├── step1_merge_mapping.py          # Near-duplicate alignment mapping
│   ├── merge_to_canonical.tsv          # Alignment metadata for 18,455 MERGE candidates
│   ├── merge_unresolved.md             # Discarded record audit log
│   ├── step2_build_training_set.py     # Training corpus builder (16 columns)
│   ├── hsv_training_set.tsv            # Master 9,298-protein training dataset
│   ├── hsv_training_set_validation.log # Training set sanity check logs
│   ├── step3_generate_splits.py        # Precomputed split generator
│   ├── within_species_groupsplit.json  # 100-fold grouped CV splits
│   ├── leave_species_out.json          # Zero-shot cross-species transfer splits
│   ├── step4_training_summary.py       # Comprehensive training corpus summary
│   ├── training_set_summary.md         # Statistical matrix report
│   ├── step4_validate_all.py           # End-to-end dataset integrity audit
│   ├── step5_verify_manifest.py        # Matrix shape & NaN validation
│   ├── VALIDATION_REPORT.md            # Master validation certificate
│   └── record_ids.json                 # 9,298 aligned sequence index manifest
├── features/                           # Feature extraction & matrix artifacts
│   ├── extract_physchem.py             # 26-d physicochemical feature extractor
│   ├── physchem_features.npy           # (9298 x 26) float32 matrix
│   ├── physchem_feature_names.json     # Feature name list
│   ├── extract_kmer_tfidf.py           # 128-d k-mer TF-IDF + SVD extractor
│   ├── kmer_tfidf_features.npy         # (9298 x 128) float32 matrix
│   ├── extract_protbert.py             # 1024-d ProtBERT GPU extractor
│   ├── protbert_embeddings.npy         # (9298 x 1024) float32 matrix
│   └── feature_extraction_report.md    # Extraction telemetry & audit report
├── benchmarking/                       # ML models, benchmarking, & evaluations
│   ├── run_benchmarking.py             # Full CV & cross-species benchmark
│   ├── within_species_results.csv      # 100-fold CV results per model/feature
│   ├── cross_species_results.csv       # Cross-species transfer results
│   ├── feature_comparison.md           # Model ablation table
│   ├── benchmark_report.md             # Benchmark summary report
│   ├── run_visualizations_and_analysis.py # Manifold & statistical testing pipeline
│   ├── error_analysis.csv              # Top 20 most confident prediction errors
│   ├── statistical_tests.csv           # Paired t-test and Wilcoxon p-values
│   ├── curation_summary_stats.json     # High-level curation counts
│   └── curation_details_stats.json     # Flag and warning breakdowns
├── results/figures/                    # 300 DPI publication figures
│   ├── fig1_confusion_matrix.png       # Normalized confusion matrix (grouped CV)
│   ├── fig2_roc_pr_curves.png          # Multi-class ROC and PR curves
│   ├── fig3_feature_comparison.png     # Feature comparison bar chart
│   ├── fig4_cross_species_transfer.png # Cross-species transfer gap
│   ├── fig5_umap_tsne_protbert.png     # Side-by-side UMAP & t-SNE projections
│   └── fig6_physchem_importance.png    # Top 15 Random Forest Gini feature importances
└── docs/                               # Comprehensive project documentation
    └── DOCUMENTATION.md                # Detailed repository manifest for every file
```

---

## Installation & Prerequisites

### Environment Setup
- Python: **3.10** or **3.11** (tested on 3.11.1)
- CUDA: **12.1+** with compatible NVIDIA GPU (tested on NVIDIA GeForce RTX 4050, 6 GB VRAM)

```bash
# Clone the repository
git clone https://github.com/Yuyutsu01/Herpes-Simplex-Virus-Analysis.git
cd Herpes-Simplex-Virus-Analysis

# Create and activate a virtual environment
python -m venv venv
venv\Scripts\activate  # On Windows
# source venv/bin/activate  # On Linux/macOS

# Install PyTorch with CUDA 12.1 acceleration
pip install --upgrade torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121

# Install required dependencies
pip install biopython pandas numpy scipy scikit-learn transformers umap-learn matplotlib seaborn tabulate
```

---

## How to Reproduce

All scripts are completely idempotent and configured to run directly from the repository root:

```bash
# ====================================================================
# PHASE 1: Curation & Canonical Reference Construction
# ====================================================================
python curation/run_full_curation.py
python curation/step1_verify_lengths.py
python curation/step2_build_annotation.py
python curation/step3_export_fasta.py

# ====================================================================
# PHASE 2: Training Corpus Expansion & Split Manifests
# ====================================================================
python training_set/step1_merge_mapping.py
python training_set/step2_build_training_set.py
python training_set/step3_generate_splits.py
python training_set/step4_training_summary.py
python training_set/step4_validate_all.py

# ====================================================================
# PHASE 3: Feature Extraction (GPU Recommended for ProtBERT)
# ====================================================================
python features/extract_physchem.py
python features/extract_kmer_tfidf.py
python features/extract_protbert.py
python training_set/step5_verify_manifest.py

# ====================================================================
# PHASE 4: Benchmarking & Manuscript Visualizations
# ====================================================================
python benchmarking/run_benchmarking.py
python benchmarking/run_visualizations_and_analysis.py
```

---

## Data Availability Statement

All data generated by this study are freely accessible within this repository:
1. **Canonical Reference Proteome:** 149 non-redundant reference proteins (75 HSV-1, 74 HSV-2) are provided in [`curation/hsv1_hsv2_curated_reference.fasta`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv1_hsv2_curated_reference.fasta) and annotated in [`curation/hsv_annotation.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv_annotation.tsv).
2. **Expanded Training Corpus:** The full 9,298-protein labeled dataset is provided in [`training_set/hsv_training_set.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/hsv_training_set.tsv), with strict one-to-one index alignment against [`training_set/record_ids.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/record_ids.json) and all feature matrices in [`features/`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/).
3. **Reference Genomes:** Canonical RefSeq GenBank records for HSV-1 strain 17 (`NC_001806.2`) and HSV-2 strain HG52 (`NC_001798.2`) are cached in [`data/refseq_cache/`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/data/refseq_cache/).