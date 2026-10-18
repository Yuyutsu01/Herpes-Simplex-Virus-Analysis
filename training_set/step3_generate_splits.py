#!/usr/bin/env python3
"""
step3_generate_splits.py

DELIVERABLE 3: Precompute group-aware train/test split manifests.

Outputs:
  A. within_species_groupsplit.json
     - Within-species 5-fold CV repeated 10 times (50 folds per species).
     - Stratified by functional_class (5 classes).
     - Grouped by ortholog_group (zero group leakage).
     - Contains record_ids for train and test.
     - Seed = 42 stored in file.

  B. leave_species_out.json
     - Two splits: HSV-1 -> HSV-2 and HSV-2 -> HSV-1.
     - record_ids for train and test.
     - Logs class distributions of train and test.
"""

import json
import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TRAINING_SET_DIR = SCRIPT_DIR


def generate_splits():
    print("Generating precomputed split manifests...")
    train_df = pd.read_csv(os.path.join(TRAINING_SET_DIR, "hsv_training_set.tsv"), sep="\t")
    print(f"Loaded training set with {len(train_df)} records.")
    
    # ------------------------------------------------------------------
    # MANIFEST A: within_species_groupsplit.json
    # ------------------------------------------------------------------
    SEED = 42
    N_SPLITS = 5
    N_REPEATS = 10
    
    within_species_manifest = {
        "metadata": {
            "description": "Within-species group-aware 5-fold cross-validation repeated 10 times",
            "seed": SEED,
            "n_splits": N_SPLITS,
            "n_repeats": N_REPEATS,
            "grouping_key": "ortholog_group",
            "stratification_key": "functional_class",
            "species_included": ["HSV-1", "HSV-2"]
        },
        "splits": {}
    }
    
    for species in ["HSV-1", "HSV-2"]:
        sp_df = train_df[train_df["species"] == species].copy().reset_index(drop=True)
        print(f"\nProcessing {species} ({len(sp_df)} records across {sp_df['ortholog_group'].nunique()} groups)...")
        
        # Get unique groups with their functional_class
        group_df = sp_df[["ortholog_group", "functional_class"]].drop_duplicates().reset_index(drop=True)
        print(f"{species} group-level class distribution:")
        print(group_df["functional_class"].value_counts().to_dict())
        
        sp_repeats = []
        
        for rep in range(N_REPEATS):
            rep_seed = SEED + rep
            sgkf = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=rep_seed)
            
            # Perform split at the group level to ensure clean grouping and stratification
            folds_data = []
            for fold_idx, (train_grp_idx, test_grp_idx) in enumerate(
                sgkf.split(group_df, group_df["functional_class"], group_df["ortholog_group"])
            ):
                train_groups = set(group_df.iloc[train_grp_idx]["ortholog_group"])
                test_groups = set(group_df.iloc[test_grp_idx]["ortholog_group"])
                
                # Verify zero group leakage
                assert len(train_groups.intersection(test_groups)) == 0, f"Leakage detected in {species} rep {rep} fold {fold_idx}"
                
                train_records = sp_df[sp_df["ortholog_group"].isin(train_groups)]["record_id"].tolist()
                test_records = sp_df[sp_df["ortholog_group"].isin(test_groups)]["record_id"].tolist()
                
                # Check class representation in test fold
                test_subset = sp_df[sp_df["record_id"].isin(test_records)]
                test_class_dist = test_subset["functional_class"].value_counts().to_dict()
                
                folds_data.append({
                    "fold_index": fold_idx,
                    "train_record_ids": train_records,
                    "test_record_ids": test_records,
                    "train_groups": sorted(list(train_groups)),
                    "test_groups": sorted(list(test_groups)),
                    "train_count": len(train_records),
                    "test_count": len(test_records),
                    "test_class_distribution": test_class_dist
                })
                
            sp_repeats.append({
                "repeat_index": rep,
                "repeat_seed": rep_seed,
                "folds": folds_data
            })
            
        within_species_manifest["splits"][species] = {
            "total_records": len(sp_df),
            "total_groups": int(sp_df["ortholog_group"].nunique()),
            "repeats": sp_repeats
        }
        
    out_within_path = os.path.join(TRAINING_SET_DIR, "within_species_groupsplit.json")
    with open(out_within_path, "w", encoding="utf-8") as f:
        json.dump(within_species_manifest, f, indent=2)
    print(f"\nSaved {out_within_path}.")
    
    # ------------------------------------------------------------------
    # MANIFEST B: leave_species_out.json
    # ------------------------------------------------------------------
    hsv1_records = train_df[train_df["species"] == "HSV-1"]["record_id"].tolist()
    hsv2_records = train_df[train_df["species"] == "HSV-2"]["record_id"].tolist()
    
    hsv1_class_dist = train_df[train_df["species"] == "HSV-1"]["functional_class"].value_counts().to_dict()
    hsv2_class_dist = train_df[train_df["species"] == "HSV-2"]["functional_class"].value_counts().to_dict()
    
    leave_species_manifest = {
        "metadata": {
            "description": "Cross-species held-out generalization splits (HSV-1 vs HSV-2)",
            "note": "Stratification across species is biologically fixed by isolate availability; class distributions logged."
        },
        "splits": {
            "hsv1_to_hsv2": {
                "description": "Train on all HSV-1 records, evaluate zero-shot transfer on all HSV-2 records",
                "train_species": "HSV-1",
                "test_species": "HSV-2",
                "train_record_ids": hsv1_records,
                "test_record_ids": hsv2_records,
                "train_count": len(hsv1_records),
                "test_count": len(hsv2_records),
                "train_class_distribution": hsv1_class_dist,
                "test_class_distribution": hsv2_class_dist
            },
            "hsv2_to_hsv1": {
                "description": "Train on all HSV-2 records, evaluate zero-shot transfer on all HSV-1 records",
                "train_species": "HSV-2",
                "test_species": "HSV-1",
                "train_record_ids": hsv2_records,
                "test_record_ids": hsv1_records,
                "train_count": len(hsv2_records),
                "test_count": len(hsv1_records),
                "train_class_distribution": hsv2_class_dist,
                "test_class_distribution": hsv1_class_dist
            }
        }
    }
    
    out_leave_path = os.path.join(TRAINING_SET_DIR, "leave_species_out.json")
    with open(out_leave_path, "w", encoding="utf-8") as f:
        json.dump(leave_species_manifest, f, indent=2)
    print(f"Saved {out_leave_path}.")
    print("\nDeliverable 3 split manifests generation complete.")


if __name__ == "__main__":
    generate_splits()
