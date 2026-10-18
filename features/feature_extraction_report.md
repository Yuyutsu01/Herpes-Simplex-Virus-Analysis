# Feature Extraction Audit & Manifest Report

This report documents the extraction, dimensional scaling, precision validation, and hardware telemetry for all feature representations extracted from [hsv_training_set.tsv](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/hsv_training_set.tsv) ($N = 9,298$ records).

## 1. Hardware & Execution Telemetry

| Hardware / Software Component | Specification |
| :--- | :--- |
| **GPU Model** | **NVIDIA GeForce RTX 4050 Laptop GPU** |
| **Total VRAM** | **6.00 GB GDDR6 (6,438,780,928 bytes)** |
| **NVIDIA Driver** | **592.82** |
| **CUDA Driver Version** | **13.1 (Compute Capability 8.9)** |
| **PyTorch Build** | **torch 2.5.1+cu121 (CUDA 12.1 enabled)** |
| **Inference Precision** | **Mixed Precision (`torch.float16`) via `torch.autocast`** |
| **Total Pipeline Wall Time** | **9.24 minutes (554.61 seconds)** |

## 2. Feature Representation Modalities

| Modality | Output Artifact | Dimensions | Data Type | Disk Size | Extraction Time | Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Physicochemical** | [C:\Users\shiva\OneDrive\Desktop\Herpes-Simplex-Virus-Analysis\features\physchem_features.npy](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/C:\Users\shiva\OneDrive\Desktop\Herpes-Simplex-Virus-Analysis\features\physchem_features.npy) | `(9298, 26)` | `float32` | 0.92 MB | 3.95 s | **PASS** |
| **K-mer TF-IDF + SVD** | [C:\Users\shiva\OneDrive\Desktop\Herpes-Simplex-Virus-Analysis\features\kmer_tfidf_features.npy](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/C:\Users\shiva\OneDrive\Desktop\Herpes-Simplex-Virus-Analysis\features\kmer_tfidf_features.npy) | `(9298, 128)` | `float32` | 4.54 MB | 20.50 s | **PASS** |
| **ProtBERT Embeddings** | [C:\Users\shiva\OneDrive\Desktop\Herpes-Simplex-Virus-Analysis\features\protbert_embeddings.npy](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/C:\Users\shiva\OneDrive\Desktop\Herpes-Simplex-Virus-Analysis\features\protbert_embeddings.npy) | `(9298, 1024)` | `float32` | 36.32 MB | 8.84 min | **PASS** |

## 3. Modality Methodological Summaries

### A. Physicochemical Descriptors (26 Features)
Computed via `Bio.SeqUtils.ProtParam.ProteinAnalysis` across all 9,298 records:
- **Global Descriptors (6):** `sequence_length`, `molecular_weight`, `aromaticity`, `instability_index`, `gravy` (hydropathy), `isoelectric_point` ($pI$).
- **Amino Acid Molar Fractions (20):** Molar fractions ($f_a = \frac{n_a}{L}$) for standard 20 amino acids in alphabetical order (`A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y`).
- **Feature Names Manifest:** [physchem_feature_names.json](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/physchem_feature_names.json).

### B. K-mer Frequency TF-IDF with SVD (128 Features)
Captures local secondary structural motifs and conserved linear sequence elements:
- **Tokenization:** Tripeptides (3-mers) and pentapeptides (5-mers) extracted with sliding stride = 1.
- **TF-IDF Weighting:** Filtered at `min_df=5`; top 2,000 3-mers and top 5,000 5-mers (7,000 total sparse features).
- **Latent Semantic Analysis:** `TruncatedSVD(n_components=128, random_state=42)` capturing **96.21%** of cumulative explained variance.
- **Standardization:** Zero-mean, unit-variance scaled via `StandardScaler` to float32.

### C. ProtBERT Contextual Embeddings (1,024 Features)
High-capacity transformer representations from `Rostlab/prot_bert` (420M parameters):
- **GPU Acceleration:** Executed on the RTX 4050 in FP16 mixed precision, maintaining steady 90–92% GPU compute utilization with VRAM footprint $\le 1.85$ GB.
- **Standard Sequences ($\le 510$ aa):** 4,575 sequences processed directly in batches of 8 with attention-mask mean pooling.
- **Long Sequences ($> 510$ aa):** 4,723 multi-domain sequences (including giant tegument UL36 at 3,164 aa) processed using a sliding window of 480 residues with a 64-residue overlap (stride = 416). Token embeddings mean-pooled per chunk and subsequently averaged across chunks into a holistic 1024-d representation.

## 4. Dataset Alignment Manifest

- **Record IDs Manifest:** [record_ids.json](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/record_ids.json)
- Contains **9,298** record identifiers matching the exact row index of `physchem_features.npy`, `kmer_tfidf_features.npy`, and `protbert_embeddings.npy` (row $i$ corresponds to `record_ids[i]`).
- Perfectly aligned with CV fold manifests [within_species_groupsplit.json](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/within_species_groupsplit.json) and [leave_species_out.json](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/leave_species_out.json).

## 5. Quality Audit Summary

| Quality Check | Target | Observed | Status |
| :--- | :--- | :--- | :--- |
| Physicochemical NaN/Inf | 0 | 0 NaNs, 0 Infs | **PASS** |
| K-mer TF-IDF NaN/Inf | 0 | 0 NaNs, 0 Infs | **PASS** |
| ProtBERT Embeddings NaN/Inf | 0 | 0 NaNs, 0 Infs | **PASS** |
| Row Count Consistency | Exactly 9,298 per matrix | 9,298 rows each | **PASS** |
| Alignment with Training TSV | 100% record ID match | 9,298 / 9,298 matched | **PASS** |

**FEATURE EXTRACTION COMPLETE: YES**
