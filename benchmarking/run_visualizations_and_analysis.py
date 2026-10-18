#!/usr/bin/env python3
"""
run_visualizations_and_analysis.py

Generates:
  1. fig5_umap_tsne_protbert.png: Side-by-side UMAP and t-SNE plots of ProtBERT embeddings.
  2. fig6_physchem_importance.png: Horizontal bar chart of top 15 Gini importances from Random Forest.
  3. error_analysis.csv: Top 20 most confident misclassifications from ProtBERT + Logistic Regression.
  4. statistical_tests.csv: Paired t-tests and Wilcoxon tests on Macro-F1 across 100 CV folds.
"""

import json
import os
import sys
import time
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import umap
from sklearn.manifold import TSNE
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score


def main():
    print("====================================================================")
    print("STARTING VISUALIZATIONS & ERROR ANALYSIS PIPELINE")
    print("====================================================================")
    t_start = time.time()
    
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
    TRAINING_SET_DIR = os.path.join(REPO_ROOT, "training_set")
    FEATURES_DIR = os.path.join(REPO_ROOT, "features")
    BENCHMARKING_DIR = SCRIPT_DIR
    FIGURES_DIR = os.path.join(REPO_ROOT, "results", "figures")
    os.makedirs(FIGURES_DIR, exist_ok=True)
    
    # 1. Load Data
    print("Loading datasets and feature matrices...")
    train_df = pd.read_csv(os.path.join(TRAINING_SET_DIR, "hsv_training_set.tsv"), sep="\t")
    with open(os.path.join(TRAINING_SET_DIR, "record_ids.json"), "r", encoding="utf-8") as f:
        record_ids = json.load(f)
    rid_to_idx = {rid: i for i, rid in enumerate(record_ids)}
    
    with open(os.path.join(FEATURES_DIR, "physchem_feature_names.json"), "r", encoding="utf-8") as f:
        physchem_names = json.load(f)
        
    protbert_mat = np.load(os.path.join(FEATURES_DIR, "protbert_embeddings.npy"))
    physchem_mat = np.load(os.path.join(FEATURES_DIR, "physchem_features.npy"))
    kmer_mat = np.load(os.path.join(FEATURES_DIR, "kmer_tfidf_features.npy"))
    
    CLASS_NAMES = sorted(train_df["functional_class"].unique().tolist())
    cls_to_idx = {c: i for i, c in enumerate(CLASS_NAMES)}
    y_all = np.array([cls_to_idx[c] for c in train_df["functional_class"]], dtype=np.int64)
    
    # Elegant, distinguishable color palette for the 5 classes
    class_palette = {
        "Capsid_Assembly": "#2b5c8f",        # Slate Deep Blue
        "Glycoprotein_Envelope": "#2ca02c",  # Forest Green
        "Regulatory": "#d62728",             # Crimson Red
        "Replication": "#ff7f0e",            # Amber Orange
        "Tegument": "#9467bd"                # Royal Purple
    }
    
    # ------------------------------------------------------------------
    # DELIVERABLE 1: UMAP & t-SNE OF PROTBERT EMBEDDINGS
    # ------------------------------------------------------------------
    print("\n====================================================================")
    print("DELIVERABLE 1: UMAP & t-SNE EMBEDDING VISUALIZATION")
    print("====================================================================")
    
    # 1. UMAP
    print("Computing UMAP projection (n_neighbors=30, min_dist=0.1, metric='cosine')...")
    t0 = time.time()
    reducer_umap = umap.UMAP(
        n_neighbors=30,
        min_dist=0.1,
        metric="cosine",
        random_state=42
    )
    umap_coords = reducer_umap.fit_transform(protbert_mat)
    print(f"  UMAP completed in {time.time() - t0:.2f} s. Shape: {umap_coords.shape}")
    
    # 2. t-SNE
    print("Computing t-SNE projection (perplexity=30, metric='cosine')...")
    t0 = time.time()
    reducer_tsne = TSNE(
        n_components=2,
        perplexity=30,
        metric="cosine",
        init="random",
        random_state=42,
        n_jobs=-1
    )
    tsne_coords = reducer_tsne.fit_transform(protbert_mat)
    print(f"  t-SNE completed in {time.time() - t0:.2f} s. Shape: {tsne_coords.shape}")
    
    # 3. Create Side-by-Side Figure
    print("Rendering fig5_umap_tsne_protbert.png (300 DPI)...")
    fig, (ax_umap, ax_tsne) = plt.subplots(1, 2, figsize=(16, 7.5), dpi=300)
    
    plot_df = pd.DataFrame({
        "UMAP-1": umap_coords[:, 0],
        "UMAP-2": umap_coords[:, 1],
        "tSNE-1": tsne_coords[:, 0],
        "tSNE-2": tsne_coords[:, 1],
        "functional_class": train_df["functional_class"]
    })
    
    # Subplot A: UMAP
    sns.scatterplot(
        data=plot_df, x="UMAP-1", y="UMAP-2", hue="functional_class",
        palette=class_palette, alpha=0.65, s=22, edgecolor="none", ax=ax_umap
    )
    ax_umap.set_title("A. UMAP Projection (Cosine Metric, N=9,298)", fontsize=13, weight="bold", pad=10)
    ax_umap.set_xlabel("UMAP Dimension 1", fontsize=11, weight="bold")
    ax_umap.set_ylabel("UMAP Dimension 2", fontsize=11, weight="bold")
    ax_umap.grid(True, linestyle="--", alpha=0.5)
    
    # Subplot B: t-SNE
    sns.scatterplot(
        data=plot_df, x="tSNE-1", y="tSNE-2", hue="functional_class",
        palette=class_palette, alpha=0.65, s=22, edgecolor="none", ax=ax_tsne
    )
    ax_tsne.set_title("B. t-SNE Projection (Cosine Metric, N=9,298)", fontsize=13, weight="bold", pad=10)
    ax_tsne.set_xlabel("t-SNE Dimension 1", fontsize=11, weight="bold")
    ax_tsne.set_ylabel("t-SNE Dimension 2", fontsize=11, weight="bold")
    ax_tsne.grid(True, linestyle="--", alpha=0.5)
    
    # Shared Legend
    handles, labels = ax_umap.get_legend_handles_labels()
    ax_umap.legend_.remove()
    ax_tsne.legend_.remove()
    fig.suptitle("ProtBERT Embeddings of HSV-1/HSV-2 Proteome", fontsize=16, weight="bold", y=0.98)
    fig.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.93),
        ncol=5, frameon=True, fontsize=10, title="HSV Functional Class",
        title_fontsize=11
    )
    plt.tight_layout(rect=[0, 0, 1, 0.88])
    fig5_path = os.path.join(FIGURES_DIR, "fig5_umap_tsne_protbert.png")
    plt.savefig(fig5_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {fig5_path}.")
    
    # ------------------------------------------------------------------
    # DELIVERABLE 2: PHYSICOCHEMICAL FEATURE IMPORTANCE
    # ------------------------------------------------------------------
    print("\n====================================================================")
    print("DELIVERABLE 2: PHYSICOCHEMICAL FEATURE IMPORTANCE")
    print("====================================================================")
    
    print("Fitting Random Forest (n_estimators=500, class_weight='balanced') on full dataset...")
    rf_global = RandomForestClassifier(
        n_estimators=500,
        class_weight="balanced",
        n_jobs=-1,
        random_state=42
    )
    rf_global.fit(physchem_mat, y_all)
    
    importances = rf_global.feature_importances_
    feat_imp_df = pd.DataFrame({
        "feature": physchem_names,
        "importance": importances
    }).sort_values(by="importance", ascending=False).reset_index(drop=True)
    
    print("Top 15 Physicochemical Features by Gini Importance:")
    print(feat_imp_df.head(15).to_string(index=False))
    
    # Clean feature names for publication plotting
    clean_labels = {
        "sequence_length": "Sequence Length (aa)",
        "molecular_weight": "Molecular Weight (Da)",
        "isoelectric_point": "Isoelectric Point (pI)",
        "gravy": "GRAVY (Hydropathy)",
        "instability_index": "Instability Index",
        "aromaticity": "Aromaticity Fraction"
    }
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        clean_labels[f"molar_fraction_{aa}"] = f"Molar Fraction: {aa}"
        
    top15_df = feat_imp_df.head(15).copy()
    top15_df["clean_name"] = top15_df["feature"].map(lambda x: clean_labels.get(x, x))
    
    print("Rendering fig6_physchem_importance.png (300 DPI)...")
    fig, ax = plt.subplots(figsize=(9, 6.5), dpi=300)
    sns.barplot(
        data=top15_df, x="importance", y="clean_name",
        palette="viridis", ax=ax, edgecolor="black", linewidth=0.8
    )
    ax.set_title("Top 15 Physicochemical Feature Importances (Random Forest, N=9,298)", fontsize=13, weight="bold", pad=12)
    ax.set_xlabel("Mean Decrease in Impurity (Gini Importance)", fontsize=11, weight="bold")
    ax.set_ylabel("Physicochemical Descriptor", fontsize=11, weight="bold")
    ax.grid(True, linestyle="--", alpha=0.5, axis="x")
    
    # Annotate bar values
    for p in ax.patches:
        val = p.get_width()
        ax.annotate(f"{val:.3f}", (val + 0.003, p.get_y() + p.get_height() / 2.),
                    ha="left", va="center", fontsize=9, color="black", weight="bold")
        
    ax.set_xlim([0, top15_df["importance"].max() * 1.15])
    plt.tight_layout()
    fig6_path = os.path.join(FIGURES_DIR, "fig6_physchem_importance.png")
    plt.savefig(fig6_path, dpi=300)
    plt.close()
    print(f"Saved {fig6_path}.")
    
    # ------------------------------------------------------------------
    # DELIVERABLE 3: ERROR ANALYSIS TABLE
    # ------------------------------------------------------------------
    print("\n====================================================================")
    print("DELIVERABLE 3: ERROR ANALYSIS TABLE")
    print("====================================================================")
    
    print("Training ProtBERT + Logistic Regression on full dataset to evaluate confidence landscape...")
    lr_full = LogisticRegression(
        class_weight="balanced",
        max_iter=2000,
        solver="lbfgs",
        random_state=42
    )
    lr_full.fit(protbert_mat, y_all)
    
    probs_full = lr_full.predict_proba(protbert_mat)
    preds_full = np.argmax(probs_full, axis=1)
    confidences = np.max(probs_full, axis=1)
    
    misclass_mask = (preds_full != y_all)
    n_errors = misclass_mask.sum()
    error_rate = n_errors / len(y_all) * 100
    print(f"Total misclassifications on full training set: {n_errors} / {len(y_all)} ({error_rate:.2f}%)")
    
    error_df = pd.DataFrame({
        "record_id": train_df["record_id"],
        "species": train_df["species"],
        "gene": train_df["gene"],
        "true_class": [CLASS_NAMES[y] for y in y_all],
        "predicted_class": [CLASS_NAMES[p] for p in preds_full],
        "confidence": confidences
    })[misclass_mask].sort_values(by="confidence", ascending=False).reset_index(drop=True)
    
    top20_errors = error_df.head(20).copy()
    top20_errors["confidence"] = top20_errors["confidence"].round(4)
    
    error_csv_path = os.path.join(BENCHMARKING_DIR, "error_analysis.csv")
    top20_errors.to_csv(error_csv_path, index=False)
    print(f"Saved {error_csv_path} with top 20 most confident misclassifications.")
    
    print("\n--- TOP 10 MOST CONFIDENT MISCLASSIFICATIONS ---")
    print(top20_errors.head(10).to_string(index=False))
    
    # Confusion pattern analysis
    print("\n--- SYSTEMATIC GENE-FAMILY CONFUSION SUMMARY ---")
    gene_confusion = error_df.groupby(["gene", "true_class", "predicted_class"]).size().reset_index(name="count")
    gene_confusion = gene_confusion.sort_values(by="count", ascending=False)
    for _, row in gene_confusion.head(10).iterrows():
        print(f"  * Gene {row['gene']} ({row['true_class']}): misclassified as '{row['predicted_class']}' ({row['count']} isolate records)")
        
    # ------------------------------------------------------------------
    # DELIVERABLE 4: STATISTICAL SIGNIFICANCE (PAIRED TESTS ACROSS 100 FOLDS)
    # ------------------------------------------------------------------
    print("\n====================================================================")
    print("DELIVERABLE 4: STATISTICAL SIGNIFICANCE TESTING ACROSS 100 FOLDS")
    print("====================================================================")
    
    with open(os.path.join(TRAINING_SET_DIR, "within_species_groupsplit.json"), "r", encoding="utf-8") as f:
        within_splits = json.load(f)
        
    print("Evaluating paired fold-level Macro-F1 scores across 100 CV folds...")
    scores_protbert_lr = []
    scores_physchem_rf = []
    scores_kmer_rf = []
    
    for species in ["HSV-1", "HSV-2"]:
        for rep in within_splits["splits"][species]["repeats"]:
            for fold in rep["folds"]:
                tr_idx = [rid_to_idx[r] for r in fold["train_record_ids"]]
                te_idx = [rid_to_idx[r] for r in fold["test_record_ids"]]
                y_tr, y_te = y_all[tr_idx], y_all[te_idx]
                
                fold_seed = rep["repeat_seed"] + fold["fold_index"]
                
                # 1. ProtBERT (LR)
                m_pbert = LogisticRegression(class_weight="balanced", max_iter=2000, solver="lbfgs", random_state=fold_seed)
                m_pbert.fit(protbert_mat[tr_idx], y_tr)
                p1 = m_pbert.predict(protbert_mat[te_idx])
                scores_protbert_lr.append(f1_score(y_te, p1, average="macro", zero_division=0))
                
                # 2. Physicochemical (RF)
                m_phys = RandomForestClassifier(n_estimators=200, class_weight="balanced", n_jobs=-1, random_state=fold_seed)
                m_phys.fit(physchem_mat[tr_idx], y_tr)
                p2 = m_phys.predict(physchem_mat[te_idx])
                scores_physchem_rf.append(f1_score(y_te, p2, average="macro", zero_division=0))
                
                # 3. K-mer (RF)
                m_kmer = RandomForestClassifier(n_estimators=200, class_weight="balanced", n_jobs=-1, random_state=fold_seed)
                m_kmer.fit(kmer_mat[tr_idx], y_tr)
                p3 = m_kmer.predict(kmer_mat[te_idx])
                scores_kmer_rf.append(f1_score(y_te, p3, average="macro", zero_division=0))
                
    scores_protbert_lr = np.array(scores_protbert_lr)
    scores_physchem_rf = np.array(scores_physchem_rf)
    scores_kmer_rf = np.array(scores_kmer_rf)
    
    # Paired Test 1: ProtBERT (LR) vs Physicochemical (RF)
    t_stat1, p_val_ttest1 = stats.ttest_rel(scores_protbert_lr, scores_physchem_rf)
    w_stat1, p_val_wilcoxon1 = stats.wilcoxon(scores_protbert_lr, scores_physchem_rf)
    
    # Paired Test 2: ProtBERT (LR) vs K-mer (RF)
    t_stat2, p_val_ttest2 = stats.ttest_rel(scores_protbert_lr, scores_kmer_rf)
    w_stat2, p_val_wilcoxon2 = stats.wilcoxon(scores_protbert_lr, scores_kmer_rf)
    
    stat_rows = [
        {
            "Comparison": "ProtBERT (LR) vs. Physicochemical (RF)",
            "Model 1 Mean Macro-F1": f"{scores_protbert_lr.mean():.4f} ± {scores_protbert_lr.std():.4f}",
            "Model 2 Mean Macro-F1": f"{scores_physchem_rf.mean():.4f} ± {scores_physchem_rf.std():.4f}",
            "Delta (M1 - M2)": f"+{scores_protbert_lr.mean() - scores_physchem_rf.mean():.4f}",
            "Paired t-statistic": f"{t_stat1:.4f}",
            "p-value (t-test)": f"{p_val_ttest1:.4e}",
            "Wilcoxon statistic": f"{w_stat1:.1f}",
            "p-value (Wilcoxon)": f"{p_val_wilcoxon1:.4e}",
            "Significance": "p < 0.001 (Highly Significant)" if p_val_wilcoxon1 < 1e-3 else "Not Significant"
        },
        {
            "Comparison": "ProtBERT (LR) vs. K-mer TF-IDF (RF)",
            "Model 1 Mean Macro-F1": f"{scores_protbert_lr.mean():.4f} ± {scores_protbert_lr.std():.4f}",
            "Model 2 Mean Macro-F1": f"{scores_kmer_rf.mean():.4f} ± {scores_kmer_rf.std():.4f}",
            "Delta (M1 - M2)": f"+{scores_protbert_lr.mean() - scores_kmer_rf.mean():.4f}",
            "Paired t-statistic": f"{t_stat2:.4f}",
            "p-value (t-test)": f"{p_val_ttest2:.4e}",
            "Wilcoxon statistic": f"{w_stat2:.1f}",
            "p-value (Wilcoxon)": f"{p_val_wilcoxon2:.4e}",
            "Significance": "p < 0.001 (Highly Significant)" if p_val_wilcoxon2 < 1e-3 else "Not Significant"
        }
    ]
    
    stat_df = pd.DataFrame(stat_rows)
    stat_csv_path = os.path.join(BENCHMARKING_DIR, "statistical_tests.csv")
    stat_df.to_csv(stat_csv_path, index=False)
    print(f"Saved {stat_csv_path}.")
    print("\n--- STATISTICAL SIGNIFICANCE RESULTS ---")
    print(stat_df.to_string(index=False))
    
    t_end = time.time()
    print(f"\nAll deliverables completed in {t_end - t_start:.2f} seconds.")
    print("VISUALIZATION COMPLETE: YES")


if __name__ == "__main__":
    main()
