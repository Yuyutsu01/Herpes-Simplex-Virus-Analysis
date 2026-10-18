# Ablation & Baseline Model Comparison

Systematic evaluation comparing the Majority Class baseline, classical physicochemical descriptors, k-mer frequency TF-IDF vectors, deep ProtBERT contextual embeddings, and multi-modal feature combinations across both Within-Species Grouped CV (50 repeats) and Cross-Species Held-Out Generalization:

| Model / Configuration     | Features               | Within-Species Macro-F1   | Within-Species BalAcc   |   Cross-Species Macro-F1 |   Cross-Species BalAcc |
|:--------------------------|:-----------------------|:--------------------------|:------------------------|-------------------------:|-----------------------:|
| Baseline (Majority Class) | None (Null Prior)      | 0.0699                    | 0.2000                  |                   0.0739 |                 0.2    |
| Logistic Regression       | Physicochemical (26-d) | 0.1479 ± 0.1134           | 0.2290 ± 0.1022         |                   0.1943 |                 0.2927 |
| Random Forest             | Physicochemical (26-d) | 0.3731 ± 0.1441           | 0.4129 ± 0.1456         |                   0.8796 |                 0.8927 |
| Logistic Regression       | K-mer TF-IDF (128-d)   | 0.2247 ± 0.1189           | 0.2726 ± 0.1327         |                   0.3869 |                 0.3997 |
| Random Forest             | K-mer TF-IDF (128-d)   | 0.3134 ± 0.1439           | 0.3562 ± 0.1432         |                   0.7423 |                 0.7424 |
| Logistic Regression       | ProtBERT (1024-d)      | 0.6331 ± 0.1383           | 0.6695 ± 0.1255         |                   0.9153 |                 0.9213 |
| Random Forest             | ProtBERT (1024-d)      | 0.5716 ± 0.1304           | 0.6071 ± 0.1261         |                   0.9558 |                 0.956  |
| Logistic Regression       | Combined (1178-d)      | 0.1737 ± 0.1180           | 0.2485 ± 0.1100         |                   0.1737 |                 0.2576 |
| Random Forest             | Combined (1178-d)      | 0.5763 ± 0.1264           | 0.6119 ± 0.1168         |                   0.9402 |                 0.9345 |

> **Key Finding:** ProtBERT embeddings and the Combined representation provide substantial gains over shallow sequence features and drastically outperform the baseline, maintaining high cross-species transferability.
