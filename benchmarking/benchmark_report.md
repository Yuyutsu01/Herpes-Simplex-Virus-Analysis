# Machine Learning Benchmark Report: HSV Proteome Functional Classification

This report documents the rigorous evaluation of Logistic Regression, Random Forest, and Deep Multi-Layer Perceptron (MLP) architectures on the 9,298-record HSV-1/HSV-2 training corpus across four distinct sequence representations.

## 1. Feature Representation Matrices

| Feature Representation | Matrix Artifact | Dimension | Description |
| :--- | :--- | :--- | :--- |
| **Physicochemical** | [physchem_features.npy](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/physchem_features.npy) | `(9298, 26)` | 6 macro-descriptors (Mw, pI, GRAVY, Instability, Aromaticity, Length) + 20 AA molar fractions |
| **K-mer TF-IDF** | [kmer_tfidf_features.npy](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/kmer_tfidf_features.npy) | `(9298, 128)` | Sparse 3-mer and 5-mer TF-IDF vectorization reduced via TruncatedSVD (96.21% explained variance) |
| **ProtBERT Embeddings** | [protbert_embeddings.npy](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/protbert_embeddings.npy) | `(9298, 1024)` | Deep contextual pLM representations from `Rostlab/prot_bert` (FP16 on RTX 4050 GPU) |
| **Combined Concatenation** | Multi-modal hstack | `(9298, 1178)` | Full joint concatenation of Physicochemical + K-mer TF-IDF + ProtBERT vectors |

## 2. Best Performing Model Configuration

- **Top Architecture:** **Logistic Regression**
- **Optimal Feature Space:** **Protbert** (`1024` dimensions)
- **Within-Species Macro-F1 (Overall):** **0.6331 ± 0.1383**
- **Within-Species Balanced Accuracy:** **0.6695 ± 0.1255**
- **Matthews Correlation Coefficient (MCC):** **0.6341 ± 0.1499**

### Per-Class Performance Summary (Best Model - Overall CV)

| Functional Class | AUPRC (Mean ± Std) | Recall / Sensitivity (Mean ± Std) |
| :--- | :--- | :--- |
| **Capsid_Assembly** | 0.7850 ± 0.2255 | 0.8233 ± 0.2215 |
| **Glycoprotein_Envelope** | 0.9344 ± 0.1053 | 0.7775 ± 0.2569 |
| **Regulatory** | 0.5253 ± 0.4125 | 0.4965 ± 0.4672 |
| **Replication** | 0.8641 ± 0.1755 | 0.6905 ± 0.2954 |
| **Tegument** | 0.6742 ± 0.2244 | 0.5598 ± 0.2594 |

## 3. Within-Species Grouped CV Benchmark Results

Evaluation across 10 repeats of 5-fold cross-validation grouped by `ortholog_group` (50 evaluation folds per species, zero homology leakage):

| Features   | Model               |   Macro-F1 (Mean) |   Macro-F1 (Std) |   BalAcc (Mean) |   BalAcc (Std) |   MCC (Mean) |   MCC (Std) |
|:-----------|:--------------------|------------------:|-----------------:|----------------:|---------------:|-------------:|------------:|
| physchem   | Logistic Regression |         0.147852  |        0.113381  |        0.228984 |      0.102164  |  0.0440198   |   0.128756  |
| physchem   | Random Forest       |         0.373118  |        0.144125  |        0.412945 |      0.145561  |  0.276291    |   0.19404   |
| physchem   | MLP Classifier      |         0.067708  |        0.0300649 |        0.199725 |      0.012392  |  0.000153192 |   0.0253666 |
| kmer       | Logistic Regression |         0.224698  |        0.118862  |        0.27264  |      0.13267   |  0.0877249   |   0.16684   |
| kmer       | Random Forest       |         0.313422  |        0.143877  |        0.356162 |      0.143161  |  0.207321    |   0.187735  |
| kmer       | MLP Classifier      |         0.196698  |        0.0920427 |        0.240456 |      0.107429  |  0.0419017   |   0.158283  |
| protbert   | Logistic Regression |         0.633133  |        0.138321  |        0.669515 |      0.125478  |  0.63405     |   0.149922  |
| protbert   | Random Forest       |         0.57161   |        0.130445  |        0.60713  |      0.12609   |  0.575141    |   0.154731  |
| protbert   | MLP Classifier      |         0.593434  |        0.132769  |        0.628049 |      0.122125  |  0.599816    |   0.164921  |
| combined   | Logistic Regression |         0.173692  |        0.117994  |        0.24849  |      0.110033  |  0.0696343   |   0.142418  |
| combined   | Random Forest       |         0.576346  |        0.126355  |        0.61187  |      0.116783  |  0.572555    |   0.135765  |
| combined   | MLP Classifier      |         0.0685456 |        0.0336741 |        0.202279 |      0.0262844 |  0.00889941  |   0.0346885 |

## 4. Cross-Species Zero-Shot Transfer Results

Models trained exclusively on one species and evaluated on the other held-out species:

| Evaluation Split   | Features   | Model               |   Macro-F1 |   Balanced Accuracy |      MCC |
|:-------------------|:-----------|:--------------------|-----------:|--------------------:|---------:|
| HSV-1 -> HSV-2     | physchem   | Random Forest       |   0.84094  |            0.876417 | 0.798662 |
| HSV-2 -> HSV-1     | physchem   | Random Forest       |   0.918337 |            0.908892 | 0.894591 |
| HSV-1 -> HSV-2     | kmer       | Random Forest       |   0.749506 |            0.762024 | 0.750901 |
| HSV-2 -> HSV-1     | kmer       | Random Forest       |   0.735147 |            0.722771 | 0.721362 |
| HSV-1 -> HSV-2     | protbert   | Logistic Regression |   0.890268 |            0.89638  | 0.869801 |
| HSV-2 -> HSV-1     | protbert   | Logistic Regression |   0.940368 |            0.946268 | 0.937194 |
| HSV-1 -> HSV-2     | combined   | Random Forest       |   0.92445  |            0.924141 | 0.907354 |
| HSV-2 -> HSV-1     | combined   | Random Forest       |   0.955866 |            0.944905 | 0.955304 |

## 5. Test Set Class Distributions

Demonstrating model robustness against natural proteome class imbalances:

| species   |   Capsid_Assembly |   Glycoprotein_Envelope |   Regulatory |   Replication |   Tegument |   Total |
|:----------|------------------:|------------------------:|-------------:|--------------:|-----------:|--------:|
| HSV-1     |              1047 |                    1277 |          633 |          1810 |       1518 |    6285 |
| HSV-2     |               416 |                     846 |          228 |           757 |        766 |    3013 |

## 6. Generated Publication Figures

- [fig1_confusion_matrix.png](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/fig1_confusion_matrix.png): Aggregate confusion matrix of the best model under grouped cross-validation.
- [fig2_roc_pr_curves.png](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/fig2_roc_pr_curves.png): Multi-class One-vs-Rest ROC and Precision-Recall curves across all 5 classes.
- [fig3_feature_comparison.png](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/fig3_feature_comparison.png): Comparative bar chart across feature representations and classifiers.
- [fig4_cross_species_transfer.png](file:///c:/Users/shiva/OneDrive/Desktop/Herpes-Simplex-Virus-Analysis/fig4_cross_species_transfer.png): Cross-species transferability and generalization gap.

**BENCHMARKING COMPLETE: YES**
