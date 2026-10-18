#!/usr/bin/env python3
"""
extract_kmer_tfidf.py

STEP 3: Extract k-mer TF-IDF features for the 9,298 training sequences.

Workflow:
  1. Generate 3-mers and 5-mers as space-separated strings.
  2. Fit TfidfVectorizer:
     - 3-mers: min_df=5, max_features=2000
     - 5-mers: min_df=5, max_features=5000
  3. Concatenate the two sparse matrices (scipy.sparse.hstack).
  4. Dimensionality reduction: TruncatedSVD(n_components=128, random_state=42).
  5. Feature scaling: StandardScaler().
  
Output:
  - kmer_tfidf_features.npy (Shape: 9298 x 128, float32)
"""

import os
import sys
import time
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
TRAINING_SET_DIR = os.path.join(REPO_ROOT, "training_set")
FEATURES_DIR = SCRIPT_DIR


def extract_kmer_tfidf():
    print("Step 3: Extracting k-mer TF-IDF features from hsv_training_set.tsv...")
    t0 = time.time()
    
    # 1. Load training dataset
    train_df = pd.read_csv(os.path.join(TRAINING_SET_DIR, "hsv_training_set.tsv"), sep="\t")
    n_records = len(train_df)
    print(f"Loaded {n_records} sequences.")
    
    # 2. Extract space-separated 3-mers and 5-mers
    print("Generating 3-mer and 5-mer token strings...")
    k3_docs = []
    k5_docs = []
    
    for i, seq in enumerate(train_df["sequence"]):
        seq_str = str(seq).strip().upper()
        L = len(seq_str)
        
        # 3-mers
        k3_tokens = [seq_str[j:j+3] for j in range(L - 2)]
        k3_docs.append(" ".join(k3_tokens))
        
        # 5-mers
        k5_tokens = [seq_str[j:j+5] for j in range(L - 4)]
        k5_docs.append(" ".join(k5_tokens))
        
        if (i + 1) % 3000 == 0 or (i + 1) == n_records:
            print(f"  Tokenized {i + 1}/{n_records} sequences...")
            
    # 3. Fit TfidfVectorizers
    print("Fitting TfidfVectorizer for 3-mers (max_features=2000, min_df=5)...")
    vec3 = TfidfVectorizer(min_df=5, max_features=2000, lowercase=False)
    X3 = vec3.fit_transform(k3_docs)
    print(f"  3-mer TF-IDF matrix shape: {X3.shape}")
    
    print("Fitting TfidfVectorizer for 5-mers (max_features=5000, min_df=5)...")
    vec5 = TfidfVectorizer(min_df=5, max_features=5000, lowercase=False)
    X5 = vec5.fit_transform(k5_docs)
    print(f"  5-mer TF-IDF matrix shape: {X5.shape}")
    
    # 4. Concatenate sparse matrices
    print("Concatenating 3-mer and 5-mer feature spaces...")
    X_concat = sp.hstack([X3, X5], format="csr")
    print(f"  Concatenated sparse shape: {X_concat.shape}")
    
    # 5. Apply TruncatedSVD (n_components=128)
    print("Applying TruncatedSVD (n_components=128, random_state=42)...")
    svd = TruncatedSVD(n_components=128, random_state=42)
    X_svd = svd.fit_transform(X_concat)
    explained_var = svd.explained_variance_ratio_.sum()
    print(f"  TruncatedSVD shape: {X_svd.shape}, Explained variance: {explained_var:.2%}")
    
    # 6. Apply StandardScaler
    print("Applying StandardScaler...")
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_svd).astype(np.float32)
    
    # 7. Validation checks
    assert X_scaled.shape == (9298, 128), f"Unexpected shape: {X_scaled.shape}"
    assert not np.isnan(X_scaled).any(), "NaN values found in k-mer TF-IDF features!"
    assert not np.isinf(X_scaled).any(), "Inf values found in k-mer TF-IDF features!"
    assert X_scaled.dtype == np.float32, f"Expected float32, got {X_scaled.dtype}"
    
    # 8. Save output
    out_path = os.path.join(FEATURES_DIR, "kmer_tfidf_features.npy")
    np.save(out_path, X_scaled)
    print(f"Saved {out_path} with shape {X_scaled.shape} and dtype {X_scaled.dtype}.")
    
    t1 = time.time()
    print(f"Step 3 completed successfully in {t1 - t0:.2f} seconds.")


if __name__ == "__main__":
    extract_kmer_tfidf()
