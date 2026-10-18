"""
step3_export_fasta.py
=====================
Deliverable 3: Export curated reference FASTA and integrity manifest.

Concepts:
- Machine-Readable Downstream Headers:
  Bioinformatics and ML tools benefit from standardized FASTA headers embedding key metadata:
  >{protein_id}|species={HSV-1|HSV-2}|gene={gene}|og={ortholog_group}|{original_header}
  Preserving the original header verbatim at the end maintains complete backwards compatibility
  and provenance tracking back to the original GenBank/UniProt/PDB source.

- Data Integrity & Deterministic Cryptographic Hashing:
  Each amino acid sequence is stripped of terminal stop codons, uppercased, formatted to standard
  80-character line wrapping, and fingerprinted using SHA-256. This ensures bit-for-bit reproducibility
  across repeated runs and across pipeline environments.
"""

import os
import csv
import hashlib
from collections import defaultdict

def wrap_fasta(seq, width=80):
    """
    Wraps sequence string to fixed width (80 chars).
    """
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def export_fasta():
    print("=== Step 3: Exporting Curated Reference FASTA ===")
    
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
    DATA_DIR = os.path.join(REPO_ROOT, "data")
    
    annotation_file = os.path.join(SCRIPT_DIR, "hsv_annotation.tsv")
    curation_file = os.path.join(SCRIPT_DIR, "hsv_curation_records.tsv")
    source_hsv1 = os.path.join(DATA_DIR, "non_redundant.fasta")
    source_hsv2 = os.path.join(DATA_DIR, "HSV2_non_redundant (2).fasta")
    
    fasta_out = os.path.join(SCRIPT_DIR, "hsv1_hsv2_curated_reference.fasta")
    manifest_out = os.path.join(SCRIPT_DIR, "hsv1_hsv2_curated_reference.manifest.tsv")
    excluded_dups_out = os.path.join(SCRIPT_DIR, "hsv1_hsv2_curated_reference.excluded_duplicates.tsv")
    
    # 1. Load annotation TSV (order of truth)
    with open(annotation_file, "r", encoding="utf-8") as f:
        annotations = list(csv.DictReader(f, delimiter="\t"))
    print(f"Loaded {len(annotations)} annotations from {annotation_file}")
    
    # 2. Map protein_id to curation record
    with open(curation_file, "r", encoding="utf-8") as f:
        curation_rows = list(csv.DictReader(f, delimiter="\t"))
        
    keep_curation = {
        r["accession"].split(".")[0]: r
        for r in curation_rows
        if r["recommended_status"] == "KEEP"
    }
    
    # Target accessions to search
    target_pids = set(a["protein_id"] for a in annotations)
    
    # 3. Read source FASTAs and track occurrences
    # Map pid -> list of (header, seq, source_file)
    extracted_records = defaultdict(list)
    
    def scan_source(path, label):
        print(f"Scanning {path}...")
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            curr_h = ""
            curr_s = []
            for line in f:
                if line.startswith(">"):
                    if curr_h:
                        header_str = curr_h
                        first_tok = header_str.split()[0].replace(">", "").strip()
                        # Extract potential accession or PDB/UniProt token
                        tokens = [first_tok]
                        if "|" in first_tok:
                            parts = first_tok.split("|")
                            tokens.extend(parts)
                            if len(parts) >= 3 and parts[0] == "pdb":
                                tokens.append(f"{parts[1]}_{parts[2]}")
                        
                        # Match with target_pids
                        matched_pid = None
                        for t in tokens:
                            clean_t = t.split(".")[0]
                            if clean_t in target_pids:
                                matched_pid = clean_t
                                break
                                
                        if matched_pid:
                            seq_str = "".join(curr_s).upper().rstrip("*")
                            extracted_records[matched_pid].append({
                                "header": header_str,
                                "seq": seq_str,
                                "source": label
                            })
                    curr_h = line[1:].strip()
                    curr_s = []
                else:
                    curr_s.append(line.strip())
            # Final record
            if curr_h:
                header_str = curr_h
                first_tok = header_str.split()[0].replace(">", "").strip()
                tokens = [first_tok]
                if "|" in first_tok:
                    parts = first_tok.split("|")
                    tokens.extend(parts)
                    if len(parts) >= 3 and parts[0] == "pdb":
                        tokens.append(f"{parts[1]}_{parts[2]}")
                matched_pid = None
                for t in tokens:
                    clean_t = t.split(".")[0]
                    if clean_t in target_pids:
                        matched_pid = clean_t
                        break
                if matched_pid:
                    seq_str = "".join(curr_s).upper().rstrip("*")
                    extracted_records[matched_pid].append({
                        "header": header_str,
                        "seq": seq_str,
                        "source": label
                    })

    # Scan both files
    scan_source(source_hsv1, "non_redundant.fasta")
    scan_source(source_hsv2, "HSV2_non_redundant (2).fasta")
    
    # 4. Check that all 149 targets were extracted
    missing_targets = target_pids - set(extracted_records.keys())
    if missing_targets:
        raise ValueError(f"CRITICAL ERROR: Failed to find sequence for targets: {missing_targets}")
        
    print(f"All {len(target_pids)} targets successfully found in source FASTA files.")
    
    # 5. Build FASTA and Manifest
    fasta_entries = []
    manifest_rows = []
    excluded_dup_rows = []
    
    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
    non_std_found = []
    
    for ann in annotations:
        pid = ann["protein_id"]
        species = ann["species"]
        gene = ann["gene"]
        og = ann["ortholog_group"]
        expected_len = int(ann["sequence_length"])
        
        matches = extracted_records[pid]
        
        # Guardrail: If multiple occurrences exist in source files, keep first, log others
        primary_match = matches[0]
        if len(matches) > 1:
            for extra in matches[1:]:
                excluded_dup_rows.append({
                    "protein_id": pid,
                    "species": species,
                    "gene": gene,
                    "duplicate_header": extra["header"],
                    "source_file": extra["source"],
                    "action": "EXCLUDED_AS_REDUNDANT_OCCURRENCE"
                })
                
        orig_header = primary_match["header"]
        seq = primary_match["seq"]
        source_label = primary_match["source"]
        
        # Check sequence length consistency
        if len(seq) != expected_len:
            raise ValueError(f"Length mismatch for {pid}: annotation={expected_len}, sequence={len(seq)}")
            
        # Check standard amino acids
        bad_chars = set(seq) - standard_aa
        if bad_chars:
            non_std_found.append((pid, bad_chars))
            
        # Compute SHA-256 hash
        sha256_hash = hashlib.sha256(seq.encode("ascii")).hexdigest()
        
        # Build structured header
        # Format: >{protein_id}|species={HSV-1|HSV-2}|gene={gene}|og={ortholog_group}|{original_header}
        new_header = f">{pid}|species={species}|gene={gene}|og={og}|{orig_header}"
        
        fasta_entries.append(f"{new_header}\n{wrap_fasta(seq)}")
        
        manifest_rows.append({
            "protein_id": pid,
            "species": species,
            "gene": gene,
            "ortholog_group": og,
            "length": len(seq),
            "header_source": source_label,
            "sha256_of_sequence": sha256_hash
        })

    # Validate standard amino acids
    if non_std_found:
        raise ValueError(f"Sequences contain non-standard characters: {non_std_found}")
        
    # Write FASTA
    with open(fasta_out, "w", encoding="utf-8") as f:
        f.write("\n".join(fasta_entries) + "\n")
    print(f"Generated {fasta_out} ({os.path.getsize(fasta_out):,} bytes, {len(fasta_entries)} records).")
    
    # Write Manifest TSV
    manifest_cols = [
        "protein_id",
        "species",
        "gene",
        "ortholog_group",
        "length",
        "header_source",
        "sha256_of_sequence"
    ]
    with open(manifest_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=manifest_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Generated {manifest_out} ({os.path.getsize(manifest_out):,} bytes).")
    
    # Write Excluded Duplicates TSV
    dup_cols = ["protein_id", "species", "gene", "duplicate_header", "source_file", "action"]
    with open(excluded_dups_out, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=dup_cols, delimiter="\t")
        writer.writeheader()
        writer.writerows(excluded_dup_rows)
    print(f"Generated {excluded_dups_out} ({len(excluded_dup_rows)} duplicate occurrences logged).")

    # Final validation checks
    assert len(fasta_entries) == 149
    hsv1_fasta_count = sum(1 for m in manifest_rows if m["species"] == "HSV-1")
    hsv2_fasta_count = sum(1 for m in manifest_rows if m["species"] == "HSV-2")
    assert hsv1_fasta_count == 75
    assert hsv2_fasta_count == 74
    print("All Step 3 validation checks passed successfully!")

if __name__ == "__main__":
    export_fasta()
