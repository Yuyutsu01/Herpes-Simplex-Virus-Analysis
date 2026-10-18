#!/usr/bin/env python3
"""
step1_merge_mapping.py

DELIVERABLE 1: Map near-duplicate MERGE records to their canonical KEEP representatives.
Criteria:
  1. Same species (HSV-1 / HSV-2) and same gene_or_protein_guess.
  2. Length ratio within 10%: abs(L_merge - L_canon) / max(L_merge, L_canon) <= 0.10.
  3. Sequence identity >= 90.0% via global pairwise alignment (Bio.Align.PairwiseAligner).
  4. Quality filter check: standard 20-AA alphabet without internal stops.

Outputs:
  - merge_to_canonical.tsv
  - merge_unresolved.md
"""

import os
import sys
import time
from collections import defaultdict
import pandas as pd
from Bio.Align import PairwiseAligner

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
DATA_DIR = os.path.join(REPO_ROOT, "data")
CURATION_DIR = os.path.join(REPO_ROOT, "curation")
TRAINING_SET_DIR = SCRIPT_DIR

# Standard 20 amino acid alphabet
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def load_raw_fasta_sequences():
    """
    Loads sequences from source FASTAs in the exact order parsed during curation.
    Returns: list of sequence strings aligned with 1-based record_index.
    """
    fasta_files = [
        os.path.join(DATA_DIR, "HSV2_non_redundant (2).fasta"),
        os.path.join(DATA_DIR, "non_redundant.fasta")
    ]
    sequences = []
    for fpath in fasta_files:
        print(f"Loading {fpath}...")
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
    print(f"Total raw sequences loaded: {len(sequences)}")
    return sequences


def load_canonical_references():
    """
    Loads the 149 canonical reference sequences and metadata.
    Returns:
      canon_map: dict of (species, gene) -> (protein_id, sequence, length)
      annot_pids: set of valid canonical protein_ids
    """
    annot = pd.read_csv(os.path.join(CURATION_DIR, "hsv_annotation.tsv"), sep="\t")
    annot_pids = set(annot["protein_id"])
    
    ref_seqs = {}
    with open(os.path.join(CURATION_DIR, "hsv1_hsv2_curated_reference.fasta"), "r", encoding="utf-8") as f:
        cid, cseq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cid:
                    ref_seqs[cid] = "".join(cseq).upper()
                cid = line[1:].split("|")[0]
                cseq = []
            else:
                cseq.append(line)
        if cid:
            ref_seqs[cid] = "".join(cseq).upper()
            
    canon_map = {}
    for _, r in annot.iterrows():
        pid = r["protein_id"]
        seq = ref_seqs[pid]
        canon_map[(r["species"], r["gene"])] = (pid, seq, len(seq))
        
    print(f"Loaded {len(canon_map)} canonical reference species-gene pairs.")
    return canon_map, annot_pids


def get_kmers(seq, k=5):
    """Generates set of k-mers for Jaccard calculation."""
    if len(seq) < k:
        return set()
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))


def compute_jaccard(seq1, seq2, k=5):
    """Computes Jaccard similarity between two sequences."""
    k1 = get_kmers(seq1, k)
    k2 = get_kmers(seq2, k)
    if not k1 or not k2:
        return 0.0
    return len(k1 & k2) / len(k1 | k2)


def run_merge_mapping():
    t0 = time.time()
    
    # 1. Load data
    raw_seqs = load_raw_fasta_sequences()
    canon_map, valid_pids = load_canonical_references()
    curation_df = pd.read_csv(os.path.join(CURATION_DIR, "hsv_curation_records.tsv"), sep="\t")
    
    # Filter for MERGE records
    merge_mask = curation_df["recommended_status"] == "MERGE"
    merge_df = curation_df[merge_mask].copy()
    total_merge_records = len(merge_df)
    print(f"Processing {total_merge_records} MERGE records...")
    
    # Configure global Needleman-Wunsch aligner
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    
    mapping_rows = []
    unresolved_records = []
    
    assigned_count = 0
    unresolved_count = 0
    clean_quality_count = 0
    
    for i, (_, row) in enumerate(merge_df.iterrows()):
        rec_idx = row["record_index"]
        acc = row["accession"]
        species = row["organism_inferred"]
        gene = row["gene_or_protein_guess"]
        merge_seq = raw_seqs[rec_idx - 1]
        lm = len(merge_seq)
        
        # Check standard amino acids and stop codons
        has_ambig = any(c not in STANDARD_AA for c in merge_seq)
        has_stop = "*" in merge_seq
        
        # Lookup candidate canonical parent
        if (species, gene) not in canon_map:
            # Species or gene unknown/mismatched
            mapping_rows.append({
                "merge_record_index": rec_idx,
                "merge_accession": acc,
                "merge_species": species,
                "merge_gene": gene,
                "canonical_protein_id": "none",
                "identity_pct": 0.00,
                "coverage_pct": 0.00,
                "assignment_status": "UNRESOLVED",
                "passes_quality_filter": "no",
                "notes": f"UNRESOLVED: No canonical representative for species '{species}' and gene '{gene}'"
            })
            unresolved_records.append({
                "record_index": rec_idx,
                "accession": acc,
                "species": species,
                "gene": gene,
                "length": lm,
                "reason": f"No canonical representative for species '{species}' and gene '{gene}'"
            })
            unresolved_count += 1
            continue
            
        canon_pid, canon_seq, lc = canon_map[(species, gene)]
        
        # Check length ratio rule: abs(lm - lc) / max(lm, lc) <= 0.10
        l_diff_ratio = abs(lm - lc) / max(lm, lc)
        if l_diff_ratio > 0.10:
            coverage = round(min(lm, lc) / lc * 100, 2)
            # Mathematical upper bound on global identity is min(lm, lc) / max(lm, lc) * 100
            max_possible_ident = round((min(lm, lc) / max(lm, lc)) * 100, 2)
            reason = f"Length ratio discrepancy ({l_diff_ratio:.1%} > 10.0%, len={lm} vs canon={lc})"
            mapping_rows.append({
                "merge_record_index": rec_idx,
                "merge_accession": acc,
                "merge_species": species,
                "merge_gene": gene,
                "canonical_protein_id": canon_pid,
                "identity_pct": max_possible_ident,
                "coverage_pct": coverage,
                "assignment_status": "UNRESOLVED",
                "passes_quality_filter": "no",
                "notes": f"UNRESOLVED: {reason}"
            })
            unresolved_records.append({
                "record_index": rec_idx,
                "accession": acc,
                "species": species,
                "gene": gene,
                "length": lm,
                "reason": reason
            })
            unresolved_count += 1
            continue
            
        # Compute k-mer Jaccard
        jaccard_val = compute_jaccard(merge_seq, canon_seq, 5)
        
        # Compute exact sequence alignment identity and coverage
        if lm == lc:
            matches = sum(1 for a, b in zip(merge_seq, canon_seq) if a == b)
            ident_pct = round(matches / lm * 100, 2)
            coverage_pct = 100.00
        else:
            aln = aligner.align(merge_seq, canon_seq)[0]
            matches = sum(1 for a, b in zip(aln[0], aln[1]) if a == b)
            ident_pct = round(matches / max(lm, lc) * 100, 2)
            coverage_pct = round(min(lm, lc) / lc * 100, 2)
            
        # Determine assignment status and quality pass
        if ident_pct >= 90.00:
            status = "ASSIGNED"
            assigned_count += 1
            
            # Quality filter check for training eligibility
            if has_stop or has_ambig:
                passes_quality = "no"
                ambig_chars = sorted(set(c for c in merge_seq if c not in STANDARD_AA))
                note_str = f"review_dropped: ambiguous residues {ambig_chars} detected" if has_ambig else "review_dropped: internal stop codon"
            else:
                passes_quality = "yes"
                clean_quality_count += 1
                note_str = f"Near-duplicate of canonical {canon_pid} (identity={ident_pct}%, jaccard={jaccard_val:.3f})"
                
            mapping_rows.append({
                "merge_record_index": rec_idx,
                "merge_accession": acc,
                "merge_species": species,
                "merge_gene": gene,
                "canonical_protein_id": canon_pid,
                "identity_pct": ident_pct,
                "coverage_pct": coverage_pct,
                "assignment_status": status,
                "passes_quality_filter": passes_quality,
                "notes": note_str
            })
        else:
            status = "UNRESOLVED"
            unresolved_count += 1
            reason = f"Identity {ident_pct}% < 90.0% threshold (jaccard={jaccard_val:.3f})"
            mapping_rows.append({
                "merge_record_index": rec_idx,
                "merge_accession": acc,
                "merge_species": species,
                "merge_gene": gene,
                "canonical_protein_id": canon_pid,
                "identity_pct": ident_pct,
                "coverage_pct": coverage_pct,
                "assignment_status": status,
                "passes_quality_filter": "no",
                "notes": f"UNRESOLVED: {reason}"
            })
            unresolved_records.append({
                "record_index": rec_idx,
                "accession": acc,
                "species": species,
                "gene": gene,
                "length": lm,
                "reason": reason
            })
            
    # Write merge_to_canonical.tsv
    out_df = pd.DataFrame(mapping_rows)
    cols = [
        "merge_record_index", "merge_accession", "merge_species", "merge_gene",
        "canonical_protein_id", "identity_pct", "coverage_pct",
        "assignment_status", "passes_quality_filter", "notes"
    ]
    out_df = out_df[cols]
    tsv_out_path = os.path.join(TRAINING_SET_DIR, "merge_to_canonical.tsv")
    out_df.to_csv(tsv_out_path, sep="\t", index=False)
    print(f"Saved {tsv_out_path} with {len(out_df)} rows.")
    
    # ------------------------------------------------------------------
    # VALIDATION CHECKS FOR DELIVERABLE 1
    # ------------------------------------------------------------------
    v_total_rows = len(out_df) == total_merge_records
    
    assigned_df = out_df[out_df["assignment_status"] == "ASSIGNED"]
    v_assigned_ident = (assigned_df["identity_pct"] >= 90.0).all()
    
    # Check that canonical_protein_id matches species
    annot_df = pd.read_csv(os.path.join(CURATION_DIR, "hsv_annotation.tsv"), sep="\t")
    pid_to_species = dict(zip(annot_df["protein_id"], annot_df["species"]))
    v_assigned_species = all(
        pid_to_species.get(pid) == sp
        for pid, sp in zip(assigned_df["canonical_protein_id"], assigned_df["merge_species"])
    )
    v_pids_valid = set(assigned_df["canonical_protein_id"]).issubset(valid_pids)
    
    print("\n--- DELIVERABLE 1 VALIDATION ---")
    print(f"1. Total rows == MERGE records ({total_merge_records}): {'PASS' if v_total_rows else 'FAIL'}")
    print(f"2. ASSIGNED rows identity >= 90.0%: {'PASS' if v_assigned_ident else 'FAIL'}")
    print(f"3. ASSIGNED rows match species: {'PASS' if v_assigned_species else 'FAIL'}")
    print(f"4. All canonical_protein_ids exist in hsv_annotation.tsv: {'PASS' if v_pids_valid else 'FAIL'}")
    print(f"ASSIGNED records: {len(assigned_df)}")
    print(f"  - Passes quality filter (training eligible): {clean_quality_count}")
    print(f"  - Fails quality filter (ambiguous/stops): {len(assigned_df) - clean_quality_count}")
    print(f"UNRESOLVED records: {unresolved_count}")
    
    if not (v_total_rows and v_assigned_ident and v_assigned_species and v_pids_valid):
        print("ERROR: Deliverable 1 validation failed!", file=sys.stderr)
        sys.exit(1)
        
    # Write merge_unresolved.md
    unres_df = pd.DataFrame(unresolved_records)
    md_out_path = os.path.join(TRAINING_SET_DIR, "merge_unresolved.md")
    with open(md_out_path, "w", encoding="utf-8") as f:
        f.write("# HSV MERGE Unresolved Records Log\n\n")
        f.write(f"Total MERGE records evaluated: **{total_merge_records}**\\\n")
        f.write(f"Successfully ASSIGNED records: **{len(assigned_df)}**\\\n")
        f.write(f"UNRESOLVED records (dropped from assignment): **{unresolved_count}**\n\n")
        
        # Categorize unresolved reasons
        f.write("## Summary of Unresolved Categories\n\n")
        len_reasons = sum(1 for r in unresolved_records if "Length ratio discrepancy" in r["reason"])
        ident_reasons = sum(1 for r in unresolved_records if "Identity" in r["reason"])
        no_canon_reasons = sum(1 for r in unresolved_records if "No canonical representative" in r["reason"])
        
        f.write(f"- **Length Ratio > 10% (Fragments/Extensions):** {len_reasons} records\n")
        f.write(f"- **Sequence Identity < 90.0%:** {ident_reasons} records\n")
        f.write(f"- **No Canonical Representative (Species/Gene):** {no_canon_reasons} records\n\n")
        
        f.write("## Unresolved Records Table (First 100 Sample)\n\n")
        f.write("| Record Index | Accession | Species | Gene | Length | Reason |\n")
        f.write("| --- | --- | --- | --- | --- | --- |\n")
        for r in unresolved_records[:100]:
            f.write(f"| {r['record_index']} | {r['accession']} | {r['species']} | {r['gene']} | {r['length']} | {r['reason']} |\n")
        if len(unresolved_records) > 100:
            f.write(f"\n*(Truncated: {len(unresolved_records) - 100} additional records omitted for brevity)*\n")
            
    print(f"Saved {md_out_path}.")
    t1 = time.time()
    print(f"Step 1 completed successfully in {t1 - t0:.2f} seconds.")


if __name__ == "__main__":
    run_merge_mapping()
