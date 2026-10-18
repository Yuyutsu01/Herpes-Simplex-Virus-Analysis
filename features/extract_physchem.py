#!/usr/bin/env python3
"""
extract_physchem.py

STEP 2: Extract physicochemical features for the 9,298 training sequences.

Features (26 total, float32):
  1. sequence_length (float)
  2. molecular_weight (float)
  3. aromaticity (float)
  4. instability_index (float)
  5. gravy (float)
  6. isoelectric_point (float)
  7-26. 20 amino acid molar fractions (A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)

Outputs:
  - physchem_features.npy (Shape: 9298 x 26, float32)
  - physchem_feature_names.json
"""

import json
import os
import sys
import time
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Alphabetical list of 20 standard amino acids
AA_LIST = sorted(list("ACDEFGHIKLMNPQRSTVWY"))

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
TRAINING_SET_DIR = os.path.join(REPO_ROOT, "training_set")
FEATURES_DIR = SCRIPT_DIR


def extract_physicochemical_features():
    print("Step 2: Extracting physicochemical features from hsv_training_set.tsv...")
    t0 = time.time()
    
    # 1. Load training dataset
    train_df = pd.read_csv(os.path.join(TRAINING_SET_DIR, "hsv_training_set.tsv"), sep="\t")
    n_records = len(train_df)
    print(f"Loaded {n_records} training records.")
    
    # Define feature names
    feature_names = [
        "sequence_length",
        "molecular_weight",
        "aromaticity",
        "instability_index",
        "gravy",
        "isoelectric_point"
    ] + [f"molar_fraction_{aa}" for aa in AA_LIST]
    
    print(f"Total features per sequence: {len(feature_names)} (expected: 26)")
    assert len(feature_names) == 26, "Feature count must be exactly 26"
    
    features_matrix = np.zeros((n_records, 26), dtype=np.float32)
    
    # 2. Extract features per sequence
    for i, seq in enumerate(train_df["sequence"]):
        seq_str = str(seq).strip().upper()
        L = len(seq_str)
        anal = ProteinAnalysis(seq_str)
        
        mw = anal.molecular_weight()
        aro = anal.aromaticity()
        instab = anal.instability_index()
        gravy = anal.gravy()
        pi = anal.isoelectric_point()
        
        # Amino acid counts to exact molar fractions
        aa_counts = anal.count_amino_acids()
        aa_fractions = [aa_counts[aa] / L for aa in AA_LIST]
        
        row_features = [float(L), mw, aro, instab, gravy, pi] + aa_fractions
        features_matrix[i, :] = np.array(row_features, dtype=np.float32)
        
        if (i + 1) % 2000 == 0 or (i + 1) == n_records:
            print(f"  Processed {i + 1}/{n_records} sequences...")
            
    # 3. Validation checks
    assert features_matrix.shape == (9298, 26), f"Unexpected shape: {features_matrix.shape}"
    assert not np.isnan(features_matrix).any(), "NaN values found in physicochemical features!"
    assert not np.isinf(features_matrix).any(), "Inf values found in physicochemical features!"
    
    # 4. Save outputs
    npy_path = os.path.join(FEATURES_DIR, "physchem_features.npy")
    np.save(npy_path, features_matrix)
    print(f"Saved {npy_path} with shape {features_matrix.shape} and dtype {features_matrix.dtype}.")
    
    names_path = os.path.join(FEATURES_DIR, "physchem_feature_names.json")
    with open(names_path, "w", encoding="utf-8") as f:
        json.dump(feature_names, f, indent=2)
    print(f"Saved {names_path}.")
    
    t1 = time.time()
    print(f"Step 2 completed successfully in {t1 - t0:.2f} seconds.")


if __name__ == "__main__":
    extract_physicochemical_features()
