#!/usr/bin/env python3
"""
step2_build_training_set.py

DELIVERABLE 2: Build the full labeled training corpus (hsv_training_set.tsv).

Rows to include:
  A. All 149 KEEP records from hsv_annotation.tsv (provenance = "canonical")
  B. Every MERGE record with assignment_status == "ASSIGNED" and passes_quality_filter == "yes"
     (provenance = "merged_isolate", labels inherited from canonical parent)

Rows to exclude:
  - any record with recommended_status in {EXCLUDE, REVIEW}
  - any MERGE record with assignment_status == "UNRESOLVED"
  - any record with passes_quality_filter == "no" (ambiguous residues > 0% or stop codons)

Columns (exactly 16 columns in specified order):
  record_id, accession, species, gene, protein_name, sequence, sequence_length,
  temporal_class, functional_class, membrane, ortholog_group, provenance,
  canonical_protein_id, identity_to_canonical_pct, source_file, notes

5-class functional scheme:
  Regulatory           : RL2, RS1, UL54, US1, US12, RL1
  Replication          : UL2, UL5, UL8, UL9, UL12, UL23, UL29, UL30, UL39, UL40, UL42, UL50, UL52, UL70
  Capsid_Assembly      : UL6, UL15, UL17, UL18, UL19, UL25, UL26, UL26.5, UL28, UL31, UL32, UL33, UL34, UL35, UL38
  Tegument             : UL3, UL4, UL7, UL11, UL13, UL14, UL16, UL21, UL24, UL36, UL37, UL41, UL46, UL47, UL48, UL49, UL51, UL55, UL56, US2, US3, US10, US11
  Glycoprotein_Envelope: UL1, UL10, UL20, UL22, UL27, UL43, UL44, UL45, UL49A, UL53, US4, US5, US6, US7, US8, US8A, US9
"""

import os
import sys
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data")
CURATION_DIR = os.path.join(REPO_ROOT, "curation")
TRAINING_SET_DIR = SCRIPT_DIR

# Define the standard 20 amino acids
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# 5-class functional mapping
FIVE_CLASS_MAP = {
    "Regulatory": ["RL2", "RS1", "UL54", "US1", "US12", "RL1"],
    "Replication": [
        "UL2", "UL5", "UL8", "UL9", "UL12", "UL23", "UL29",
        "UL30", "UL39", "UL40", "UL42", "UL50", "UL52", "UL70"
    ],
    "Capsid_Assembly": [
        "UL6", "UL15", "UL17", "UL18", "UL19", "UL25", "UL26",
        "UL26.5", "UL28", "UL31", "UL32", "UL33", "UL34", "UL35", "UL38"
    ],
    "Tegument": [
        "UL3", "UL4", "UL7", "UL11", "UL13", "UL14", "UL16", "UL21",
        "UL24", "UL36", "UL37", "UL41", "UL46", "UL47", "UL48", "UL49",
        "UL51", "UL55", "UL56", "US2", "US3", "US10", "US11"
    ],
    "Glycoprotein_Envelope": [
        "UL1", "UL10", "UL20", "UL22", "UL27", "UL43", "UL44", "UL45",
        "UL49A", "UL53", "US4", "US5", "US6", "US7", "US8", "US8A", "US9"
    ]
}

# Invert mapping: gene -> 5-class functional label
GENE_TO_5CLASS = {}
for func_class, genes in FIVE_CLASS_MAP.items():
    for g in genes:
        GENE_TO_5CLASS[g] = func_class


def load_raw_fasta_sequences():
    """Loads all sequences from original FASTA sources in index order."""
    fasta_files = [
        os.path.join(DATA_DIR, "HSV2_non_redundant (2).fasta"),
        os.path.join(DATA_DIR, "non_redundant.fasta")
    ]
    sequences = []
    for fpath in fasta_files:
        with open(fpath, "r", encoding="utf-8", errors="ignore") as f:
            curr_h = None
            curr_s = []
            for line in f:
                line_str = line.strip()
                if not line_str:
                    continue
                if line_str.startswith(">"):
                    if curr_h is not None:
                        sequences.append("".join(curr_s).upper())
                    curr_h = line_str[1:]
                    curr_s = []
                else:
                    curr_s.append(line_str)
            if curr_h is not None:
                sequences.append("".join(curr_s).upper())
    return sequences


def build_training_set():
    print("Building expanded training set (hsv_training_set.tsv)...")
    
    # 1. Load inputs
    curation_df = pd.read_csv(os.path.join(CURATION_DIR, "hsv_curation_records.tsv"), sep="\t")
    annot_df = pd.read_csv(os.path.join(CURATION_DIR, "hsv_annotation.tsv"), sep="\t")
    merge_map_df = pd.read_csv(os.path.join(TRAINING_SET_DIR, "merge_to_canonical.tsv"), sep="\t")
    raw_seqs = load_raw_fasta_sequences()
    
    # Pre-index canonical metadata by protein_id
    annot_by_pid = {}
    for _, r in annot_df.iterrows():
        pid = r["protein_id"]
        annot_by_pid[pid] = {
            "protein_id": pid,
            "species": r["species"],
            "gene": r["gene"],
            "protein_name": r["protein_name"],
            "temporal_class": r["temporal_class"],
            "membrane": r["membrane"],
            "ortholog_group": f"OG_{r['gene']}"
        }
        
    # Pre-index canonical record_index from curation_df where recommended_status == 'KEEP'
    keep_df = curation_df[curation_df["recommended_status"] == "KEEP"].copy()
    keep_df["base_acc"] = keep_df["accession"].apply(lambda x: x.split(".")[0])
    pid_to_curation_row = {}
    for _, r in keep_df.iterrows():
        pid_to_curation_row[r["base_acc"]] = r

    training_rows = []
    
    # ------------------------------------------------------------------
    # PART A: Canonical KEEP Records (149 rows)
    # ------------------------------------------------------------------
    print("Processing 149 canonical KEEP rows...")
    for pid, meta in annot_by_pid.items():
        cur_row = pid_to_curation_row[pid]
        rec_idx = cur_row["record_index"]
        seq = raw_seqs[rec_idx - 1]
        gene = meta["gene"]
        func_class = GENE_TO_5CLASS[gene]
        
        training_rows.append({
            "record_id": f"C_{rec_idx}",
            "accession": cur_row["accession"],
            "species": meta["species"],
            "gene": gene,
            "protein_name": meta["protein_name"],
            "sequence": seq,
            "sequence_length": len(seq),
            "temporal_class": meta["temporal_class"],
            "functional_class": func_class,
            "membrane": meta["membrane"],
            "ortholog_group": meta["ortholog_group"],
            "provenance": "canonical",
            "canonical_protein_id": pid,
            "identity_to_canonical_pct": 100.00,
            "source_file": cur_row["source_file"],
            "notes": "Canonical reference protein"
        })
        
    num_canonical = len(training_rows)
    print(f"Added {num_canonical} canonical records.")
    
    # ------------------------------------------------------------------
    # PART B: Merged Isolate Records (ASSIGNED and passes_quality_filter == 'yes')
    # ------------------------------------------------------------------
    print("Processing ASSIGNED merged isolate rows...")
    assigned_clean_df = merge_map_df[
        (merge_map_df["assignment_status"] == "ASSIGNED") &
        (merge_map_df["passes_quality_filter"] == "yes")
    ]
    print(f"Found {len(assigned_clean_df)} ASSIGNED clean isolate records.")
    
    for _, r in assigned_clean_df.iterrows():
        rec_idx = r["merge_record_index"]
        cur_row = curation_df.iloc[rec_idx - 1]
        seq = raw_seqs[rec_idx - 1]
        canon_pid = r["canonical_protein_id"]
        parent_meta = annot_by_pid[canon_pid]
        gene = parent_meta["gene"]
        func_class = GENE_TO_5CLASS[gene]
        
        training_rows.append({
            "record_id": f"M_{rec_idx}",
            "accession": r["merge_accession"],
            "species": r["merge_species"],
            "gene": gene,
            "protein_name": parent_meta["protein_name"],
            "sequence": seq,
            "sequence_length": len(seq),
            "temporal_class": parent_meta["temporal_class"],
            "functional_class": func_class,
            "membrane": parent_meta["membrane"],
            "ortholog_group": parent_meta["ortholog_group"],
            "provenance": "merged_isolate",
            "canonical_protein_id": canon_pid,
            "identity_to_canonical_pct": r["identity_pct"],
            "source_file": cur_row["source_file"],
            "notes": f"Inherited from canonical {canon_pid} (identity={r['identity_pct']}%)"
        })
        
    # Convert to DataFrame
    cols_order = [
        "record_id", "accession", "species", "gene", "protein_name",
        "sequence", "sequence_length", "temporal_class", "functional_class",
        "membrane", "ortholog_group", "provenance", "canonical_protein_id",
        "identity_to_canonical_pct", "source_file", "notes"
    ]
    train_df = pd.DataFrame(training_rows)[cols_order]
    
    # ------------------------------------------------------------------
    # VALIDATION CHECKS FOR DELIVERABLE 2
    # ------------------------------------------------------------------
    log_lines = []
    log_lines.append("=== HSV TRAINING SET VALIDATION LOG ===")
    
    # Check 1: exactly 149 rows with provenance == 'canonical'
    canon_subset = train_df[train_df["provenance"] == "canonical"]
    check1_pass = len(canon_subset) == 149
    log_lines.append(f"Check 1: Exactly 149 rows with provenance == 'canonical': {'PASS' if check1_pass else 'FAIL'} (count={len(canon_subset)})")
    
    # Check 2: every canonical row's gene matches hsv_annotation.tsv
    canon_gene_match = all(
        annot_by_pid[pid]["gene"] == g
        for pid, g in zip(canon_subset["canonical_protein_id"], canon_subset["gene"])
    )
    log_lines.append(f"Check 2: Every canonical row's gene matches hsv_annotation.tsv: {'PASS' if canon_gene_match else 'FAIL'}")
    
    # Check 3: every merged_isolate row has a valid canonical_protein_id
    merged_subset = train_df[train_df["provenance"] == "merged_isolate"]
    merged_pid_valid = set(merged_subset["canonical_protein_id"]).issubset(set(annot_by_pid.keys()))
    log_lines.append(f"Check 3: Every merged_isolate row has a valid canonical_protein_id: {'PASS' if merged_pid_valid else 'FAIL'}")
    
    # Check 4: no row has non-standard amino acids
    non_std_count = sum(
        1 for seq in train_df["sequence"]
        if any(c not in STANDARD_AA for c in seq)
    )
    check4_pass = non_std_count == 0
    log_lines.append(f"Check 4: No row has non-standard amino acids: {'PASS' if check4_pass else 'FAIL'} (non-standard count={non_std_count})")
    
    # Check 5: no row has internal '*'
    stop_count = sum(1 for seq in train_df["sequence"] if "*" in seq)
    check5_pass = stop_count == 0
    log_lines.append(f"Check 5: No row has internal '*': {'PASS' if check5_pass else 'FAIL'} (stop count={stop_count})")
    
    # Check 6: sum of rows == 149 + (# ASSIGNED MERGE rows passing quality filter)
    expected_sum = 149 + len(assigned_clean_df)
    check6_pass = len(train_df) == expected_sum
    log_lines.append(f"Check 6: Sum of rows == 149 + ASSIGNED clean MERGE rows ({expected_sum}): {'PASS' if check6_pass else 'FAIL'} (total rows={len(train_df)})")
    
    # Check 7: Class distribution breakdown
    log_lines.append("\nCounts by (species, functional_class):")
    class_matrix = pd.crosstab(train_df["species"], train_df["functional_class"], margins=True)
    log_lines.append(class_matrix.to_string())
    
    # Print validation log
    full_log_str = "\n".join(log_lines)
    print("\n" + full_log_str + "\n")
    
    # Write validation log to disk
    log_path = os.path.join(TRAINING_SET_DIR, "hsv_training_set_validation.log")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(full_log_str + "\n")
    print(f"Saved {log_path}.")
    
    # If any validation check fails, stop and raise error
    all_checks = [check1_pass, canon_gene_match, merged_pid_valid, check4_pass, check5_pass, check6_pass]
    if not all(all_checks):
        print("ERROR: One or more validation checks failed! Halting before file output.", file=sys.stderr)
        sys.exit(1)
        
    # Write hsv_training_set.tsv
    out_tsv = os.path.join(TRAINING_SET_DIR, "hsv_training_set.tsv")
    train_df.to_csv(out_tsv, sep="\t", index=False)
    print(f"SUCCESS: Saved {out_tsv} with {len(train_df)} rows and {len(train_df.columns)} columns.")


if __name__ == "__main__":
    build_training_set()
