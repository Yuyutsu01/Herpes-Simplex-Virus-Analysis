# Herpes Simplex Virus Proteome Analysis (HSV-1 & HSV-2)
## Complete Repository Manifest & Technical Documentation

This document provides a comprehensive, file-by-file specification and architectural guide for every script, dataset, manifest, report, and figure in the repository.

---

## Table of Contents
1. [`data/` — Raw Sequence Inputs & Reference Genomes](#1-data--raw-sequence-inputs--reference-genomes)
2. [`curation/` — Curation Pipeline & Canonical Reference Standard](#2-curation--curation-pipeline--canonical-reference-standard)
3. [`training_set/` — Expanded Training Corpus & Evaluation Splits](#3-training_set--expanded-training-corpus--evaluation-splits)
4. [`features/` — Feature Extraction & Representation Matrices](#4-features--feature-extraction--representation-matrices)
5. [`benchmarking/` — Machine Learning Benchmarking & Error Analyses](#5-benchmarking--machine-learning-benchmarking--error-analyses)
6. [`results/figures/` — Publication Visualizations](#6-resultsfigures--publication-visualizations)
7. [`root/` — Root Level Files & Global Compilers](#7-root--root-level-files--global-compilers)

---

## 1. `data/` — Raw Sequence Inputs & Reference Genomes

Contains raw, uncurated input FASTA files and downloaded reference genome caches:

* **[`data/non_redundant.fasta`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/data/non_redundant.fasta)**:
  - **Description:** Original source FASTA file for HSV-1 proteins collected from diverse public repositories (GenBank, UniProt, PRF, PDB, patents).
  - **Content:** Contains 11,811 raw sequence records spanning canonical proteins, clinical isolate variants, synthetic constructs, truncations, and non-HSV contaminants.
  - **File Size:** ~7.95 MB.
* **[`data/HSV2_non_redundant (2).fasta`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/data/HSV2_non_redundant%20(2).fasta)**:
  - **Description:** Original source FASTA file for HSV-2 proteins.
  - **Content:** Contains 10,899 raw sequence records spanning HSV-2 clinical isolates, laboratory strains (HG52, 333, MS), fragments, and patents.
  - **File Size:** ~7.10 MB.
* **[`data/refseq_cache/`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/data/refseq_cache/)**:
  - **Description:** Local cache folder for full NCBI RefSeq GenBank flatfiles downloaded via Entrez `efetch` (`rettype=gbwithparts`). Enables completely offline, reproducible genomic cross-verification.
  - **Files:**
    - `NC_001806.2.gbk`: Complete annotated reference genome for **Human alphaherpesvirus 1 strain 17** (152,222 bp double-stranded DNA).
    - `NC_001798.2.gbk`: Complete annotated reference genome for **Human alphaherpesvirus 2 strain HG52** (154,746 bp double-stranded DNA).

---

## 2. `curation/` — Curation Pipeline & Canonical Reference Standard

Contains the multi-phase bioinformatics curation pipeline, RefSeq verification tools, and the gold-standard 149-protein canonical reference set:

* **[`curation/run_full_curation.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/run_full_curation.py)**:
  - **Purpose:** Primary curation engine. Parses all 22,710 raw FASTA headers, cleans descriptions, infers viral taxonomy, flags contaminations, indexes exact sequence duplicates (via SHA-256), clusters near-duplicates ($\ge 90\%$ identity), and assigns recommended curation statuses (`KEEP`, `MERGE`, `EXCLUDE`, `REVIEW`).
* **[`curation/step1_verify_lengths.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/step1_verify_lengths.py)**:
  - **Purpose:** Genomic cross-verification script resolving apparent length discrepancies in candidate canonical proteins against canonical RefSeq genomes using Biopython's `PairwiseAligner`. Specifically validates HSV-2 UL31 (306 aa) and HSV-2 UL27 (904 aa).
* **[`curation/verify_lengths.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/verify_lengths.md)**:
  - **Purpose:** Narrative verification report confirming that accession `WZW62192.1` (UL31, 306 aa) and accession `QAU10510.1` (UL27 gB, 904 aa) are authentic, full-length coding sequences with authentic Met1 starts and clean stop codons.
* **[`curation/verify_lengths.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/verify_lengths.json)**:
  - **Purpose:** Machine-readable output containing alignment identities, coverages, coordinates, and `CONFIRMED_FULL_LENGTH` verdicts.
* **[`curation/step2_build_annotation.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/step2_build_annotation.py)**:
  - **Purpose:** Compiles the 12-column annotation table for the 149 canonical reference proteins (75 HSV-1, 74 HSV-2), harmonizing gene nomenclature, temporal expression classes (`IE`, `E`, `L`, `E-L`), functional classes, orthology groups (`OG_*`), and membrane localization.
* **[`curation/hsv_annotation.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv_annotation.tsv)**:
  - **Purpose:** Canonical reference annotation table (149 rows, 12 columns). Acts as the authoritative source of truth for viral gene orthology and biological properties.
* **[`curation/hsv_annotation_audit.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv_annotation_audit.md)**:
  - **Purpose:** Detailed auditing report cross-checking all 149 canonical proteins, verifying 100% proteome completeness and zero missing core genes.
* **[`curation/step3_export_fasta.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/step3_export_fasta.py)**:
  - **Purpose:** Extracts and wraps canonical sequences to standard 80-character line lengths, formats standardized FASTA headers (`>{protein_id}|species={sp}|gene={g}|og={og}|{orig_header}`), and fingerprints each sequence with SHA-256.
* **[`curation/hsv1_hsv2_curated_reference.fasta`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv1_hsv2_curated_reference.fasta)**:
  - **Purpose:** Curated reference FASTA containing the 149 canonical protein sequences.
* **[`curation/hsv1_hsv2_curated_reference.manifest.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv1_hsv2_curated_reference.manifest.tsv)**:
  - **Purpose:** Cryptographic manifest listing protein IDs, lengths, and SHA-256 hashes for all 149 canonical sequences.
* **[`curation/hsv1_hsv2_curated_reference.excluded_duplicates.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv1_hsv2_curated_reference.excluded_duplicates.tsv)**:
  - **Purpose:** Audit tracking file confirming that zero redundant sequences contaminated the canonical reference standard.
* **[`curation/hsv_curation_records.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv_curation_records.tsv)**:
  - **Purpose:** Master curation database table containing all 22,710 records with 20 structured columns.
* **[`curation/hsv_curation_table.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv_curation_table.md)**:
  - **Purpose:** Complete Markdown-formatted curation table for GitHub inspection (~7.4 MB).
* **[`curation/hsv_curation_report.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/curation/hsv_curation_report.md)**:
  - **Purpose:** Curator executive summary reporting status breakdowns, contamination flags, and the explicit KEEP list.

---

## 3. `training_set/` — Expanded Training Corpus & Evaluation Splits

Contains the label propagation pipeline that expanded the dataset to 9,298 clean sequences and precomputed cross-validation splits:

* **[`training_set/step1_merge_mapping.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/step1_merge_mapping.py)**:
  - **Purpose:** Maps 18,455 candidate `MERGE` records against their canonical reference using pairwise Needleman-Wunsch global alignment, requiring $\ge 90.0\%$ sequence identity, $\le 10.0\%$ length difference, and standard 20-amino-acid alphabet.
* **[`training_set/merge_to_canonical.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/merge_to_canonical.tsv)**:
  - **Purpose:** Alignment mapping table recording identity, coverage, and `ASSIGNED` vs. `UNRESOLVED` status for all 18,455 merge candidates.
* **[`training_set/merge_unresolved.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/merge_unresolved.md)**:
  - **Purpose:** Audit log detailing why 9,235 merge candidates were dropped (length ratio violations, $<90\%$ identity, or non-HSV contaminants).
* **[`training_set/step2_build_training_set.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/step2_build_training_set.py)**:
  - **Purpose:** Merges the 149 canonical references with 9,149 verified clinical isolate variants, propagating functional and temporal labels into `hsv_training_set.tsv`.
* **[`training_set/hsv_training_set.tsv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/hsv_training_set.tsv)**:
  - **Purpose:** Master machine learning dataset containing 9,298 high-quality viral protein sequences across 16 standardized columns.
* **[`training_set/hsv_training_set_validation.log`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/hsv_training_set_validation.log)**:
  - **Purpose:** Execution log confirming that 100% of rows satisfy length, amino acid alphabet, and functional category constraints.
* **[`training_set/step3_generate_splits.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/step3_generate_splits.py)**:
  - **Purpose:** Precomputes group-aware, stratified cross-validation and cross-species evaluation splits using `StratifiedGroupKFold`.
* **[`training_set/within_species_groupsplit.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/within_species_groupsplit.json)**:
  - **Purpose:** Precomputed 100-fold evaluation manifest ($2 \text{ species} \times 10 \text{ repeats} \times 5 \text{ folds}$) grouped by `ortholog_group` with zero group leakage.
* **[`training_set/leave_species_out.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/leave_species_out.json)**:
  - **Purpose:** Zero-shot cross-species transfer split manifest (HSV-1 $\to$ HSV-2 and HSV-2 $\to$ HSV-1).
* **[`training_set/step4_training_summary.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/step4_training_summary.py)**:
  - **Purpose:** Computes contingency tables, class distributions, and sanity checks for the training corpus.
* **[`training_set/training_set_summary.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/training_set_summary.md)**:
  - **Purpose:** Summary report documenting class distributions and confirming `READY FOR FEATURE EXTRACTION: YES`.
* **[`training_set/step4_validate_all.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/step4_validate_all.py)**:
  - **Purpose:** Automated audit verifying 100% data integrity across all curation and training files.
* **[`training_set/VALIDATION_REPORT.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/VALIDATION_REPORT.md)**:
  - **Purpose:** Formal audit certificate verifying complete compliance of all datasets.
* **[`training_set/step5_verify_manifest.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/step5_verify_manifest.py)**:
  - **Purpose:** Validates numerical feature matrices (dimensions, dtypes, non-finite values) and exports `record_ids.json`.
* **[`training_set/record_ids.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/training_set/record_ids.json)**:
  - **Purpose:** Alignment manifest containing 9,298 record identifiers matching the exact row index of all numerical matrices.

---

## 4. `features/` — Feature Extraction & Representation Matrices

Contains feature extraction pipelines and persisted float32 NumPy matrices:

* **[`features/extract_physchem.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/extract_physchem.py)**:
  - **Purpose:** Extracts 26 physicochemical features using Biopython's `ProteinAnalysis`: 6 macro-descriptors (length, molecular weight, aromaticity, instability index, hydropathy GRAVY, isoelectric point) and 20 amino acid molar fractions.
* **[`features/physchem_features.npy`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/physchem_features.npy)**:
  - **Shape:** `(9298, 26)` | **Dtype:** `float32` | **Size:** 0.92 MB | **NaNs/Infs:** 0.
* **[`features/physchem_feature_names.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/physchem_feature_names.json)**:
  - **Purpose:** JSON list containing the 26 feature column names.
* **[`features/extract_kmer_tfidf.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/extract_kmer_tfidf.py)**:
  - **Purpose:** Generates 3-mer and 5-mer peptide tokens, computes TF-IDF representations (7,000 sparse terms), reduces dimensionality to 128 components via `TruncatedSVD` (capturing 96.21% explained variance), and standardizes using `StandardScaler`.
* **[`features/kmer_tfidf_features.npy`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/kmer_tfidf_features.npy)**:
  - **Shape:** `(9298, 128)` | **Dtype:** `float32` | **Size:** 4.54 MB | **NaNs/Infs:** 0.
* **[`features/extract_protbert.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/extract_protbert.py)**:
  - **Purpose:** GPU-accelerated deep feature extractor using `Rostlab/prot_bert` (420M parameters) in FP16 precision. Short sequences ($\le 510$ aa) are batched directly; long sequences ($> 510$ aa, including giant tegument UL36 at 3,164 aa) are processed using a sliding window of 480 residues with 64-residue overlap, mean-pooling chunk embeddings into a 1,024-d representation.
* **[`features/protbert_embeddings.npy`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/protbert_embeddings.npy)**:
  - **Shape:** `(9298, 1024)` | **Dtype:** `float32` | **Size:** 36.32 MB | **NaNs/Infs:** 0.
* **[`features/feature_extraction_report.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/features/feature_extraction_report.md)**:
  - **Purpose:** Feature extraction audit report summarizing dimensions, extraction timings, hardware utilization, and finite value checks.

---

## 5. `benchmarking/` — Machine Learning Benchmarking & Error Analyses

Contains cross-validation drivers, statistical hypothesis tests, and performance tables:

* **[`benchmarking/run_benchmarking.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/run_benchmarking.py)**:
  - **Purpose:** Master benchmarking script evaluating Logistic Regression, Random Forest, and GPU-accelerated PyTorch MLP architectures across all feature modalities over 100 CV folds and leave-species-out splits.
* **[`benchmarking/within_species_results.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/within_species_results.csv)**:
  - **Purpose:** Full evaluation results across all 36 within-species configurations (Macro-F1, Balanced Accuracy, MCC, per-class AUPRC, and Recall).
* **[`benchmarking/cross_species_results.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/cross_species_results.csv)**:
  - **Purpose:** Zero-shot transfer results for HSV-1 $\to$ HSV-2 and HSV-2 $\to$ HSV-1 across all feature sets.
* **[`benchmarking/feature_comparison.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/feature_comparison.md)**:
  - **Purpose:** Markdown table summarizing feature ablation against the Majority Class baseline.
* **[`benchmarking/benchmark_report.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/benchmark_report.md)**:
  - **Purpose:** Comprehensive benchmark report detailing optimal hyperparameters, per-class AUPRC, and confusion patterns.
* **[`benchmarking/run_visualizations_and_analysis.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/run_visualizations_and_analysis.py)**:
  - **Purpose:** Generates UMAP/t-SNE manifold projections, trains global Random Forest for Gini feature importance, executes hard-negative error analysis, and computes paired statistical significance tests across 100 folds.
* **[`benchmarking/error_analysis.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/error_analysis.csv)**:
  - **Purpose:** Top 20 most confident prediction errors from ProtBERT + Logistic Regression, documenting biological edge cases (e.g., UL15 terminase misclassified as replication).
* **[`benchmarking/statistical_tests.csv`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/statistical_tests.csv)**:
  - **Purpose:** Paired $t$-test and Wilcoxon signed-rank test statistics demonstrating extreme significance ($p < 10^{-16}$) of ProtBERT over classical features.
* **[`benchmarking/curation_summary_stats.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/curation_summary_stats.json)**:
  - **Purpose:** High-level record counts per source file, inferred taxonomy, and recommended curation status.
* **[`benchmarking/curation_details_stats.json`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/benchmarking/curation_details_stats.json)**:
  - **Purpose:** Detailed quality audit tallies (ambiguities, stop codons, length anomalies, and vector contaminations).

---

## 6. `results/figures/` — Publication Visualizations

All figures are rendered in 300 DPI high resolution with consistent typography and color palettes:

* **[`results/figures/fig1_confusion_matrix.png`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/results/figures/fig1_confusion_matrix.png)**:
  - Normalized confusion matrix of the champion ProtBERT + Logistic Regression model across within-species grouped cross-validation, demonstrating high sensitivity on glycoproteins and replication machinery.
* **[`results/figures/fig2_roc_pr_curves.png`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/results/figures/fig2_roc_pr_curves.png)**:
  - Multi-class One-vs-Rest Receiver Operating Characteristic (ROC) and Precision-Recall (PR) curves across all 5 viral functional classes.
* **[`results/figures/fig3_feature_comparison.png`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/results/figures/fig3_feature_comparison.png)**:
  - Comparative bar chart demonstrating Macro-F1 performance across feature representations (Physicochemical, K-mer TF-IDF, ProtBERT, Combined) and classifier architectures compared against the Majority Baseline.
* **[`results/figures/fig4_cross_species_transfer.png`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/results/figures/fig4_cross_species_transfer.png)**:
  - Cross-species generalization gap analysis contrasting within-species grouped cross-validation against zero-shot transfer (HSV-1 $\to$ HSV-2 and HSV-2 $\to$ HSV-1).
* **[`results/figures/fig5_umap_tsne_protbert.png`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/results/figures/fig5_umap_tsne_protbert.png)**:
  - Side-by-side two-dimensional UMAP (cosine metric, $n=30$, $\text{min\_dist}=0.1$) and t-SNE (cosine metric, perplexity=30) manifold projections of all 9,298 ProtBERT embeddings colored by functional class.
* **[`results/figures/fig6_physchem_importance.png`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/results/figures/fig6_physchem_importance.png)**:
  - Horizontal bar chart of the top 15 physicochemical descriptors ranked by Random Forest Mean Decrease in Impurity (Gini Importance), highlighting tyrosine molar fraction, aromaticity, and hydropathy (GRAVY).

---

## 7. `root/` — Root Level Files & Global Compilers

* **[`README.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/README.md)**:
  - Publication-ready overview, pipeline flowchart, benchmark summary table, reproduction commands, and data availability statement.
* **[`generate_curation_report.py`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/generate_curation_report.py)**:
  - Standalone report compiler formatting Markdown curation audit summaries.
* **[`docs/DOCUMENTATION.md`](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/docs/DOCUMENTATION.md)**:
  - This comprehensive file-by-file technical manifest.
