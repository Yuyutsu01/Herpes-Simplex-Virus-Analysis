"""
run_full_curation.py
====================
Complete Bioinformatics Curation Pipeline for HSV-1 & HSV-2 Proteome Analysis.

Concepts:
- FASTA Protein Curation: Aggregating diverse clinical isolates and legacy database
  entries requires parsing heterogeneous headers (UniProt, PDB, PRF, PIR, GenBank, Patents),
  stripping annotation clutter, inferring species strictly from textual evidence, identifying
  PCR fragments / truncations, filtering cross-species / vector / host contaminations,
  and resolving redundancy via exact duplicate hashes and near-duplicate clustering.

- Exact Duplicates vs Near-Duplicate Clusters:
  Exact duplicates share 100% identical sequence (SHA-256 hash match).
  Near-duplicates share >=90% sequence identity (due to clinical strain variation or sequencing errors).

- Canonical Proteome Coverage:
  Herpes Simplex Virus genomes contain 71 canonical genes:
  UL1 to UL56 (Unique Long region), US1 to US12 (Unique Short region),
  and RL1, RL2, RS1 (Terminal/Internal Inverted Repeats).
"""

import os
import re
import sys
import json
import hashlib
import difflib
from collections import defaultdict, Counter

# Canonical gene lengths (approximations based on Strain 17 / HG52 reference proteomes)
CANONICAL_LENGTHS = {
    'UL1': 224, 'UL2': 334, 'UL3': 224, 'UL4': 199, 'UL5': 882,
    'UL6': 676, 'UL7': 296, 'UL8': 750, 'UL9': 851, 'UL10': 473,
    'UL11': 99, 'UL12': 512, 'UL13': 518, 'UL14': 232, 'UL15': 735,
    'UL16': 373, 'UL17': 703, 'UL18': 318, 'UL19': 1374, 'UL20': 222,
    'UL21': 535, 'UL22': 838, 'UL23': 376, 'UL24': 269, 'UL25': 580,
    'UL26': 635, 'UL26.5': 329, 'UL27': 904, 'UL28': 785, 'UL29': 1196,
    'UL30': 1235, 'UL31': 306, 'UL32': 596, 'UL33': 130, 'UL34': 275,
    'UL35': 112, 'UL36': 3164, 'UL37': 1123, 'UL38': 465, 'UL39': 1137,
    'UL40': 340, 'UL41': 489, 'UL42': 488, 'UL43': 417, 'UL44': 511,
    'UL45': 172, 'UL46': 718, 'UL47': 693, 'UL48': 490, 'UL49': 301,
    'UL49A': 91, 'UL50': 371, 'UL51': 244, 'UL52': 1058, 'UL53': 338,
    'UL54': 512, 'UL55': 186, 'UL56': 186,
    'US1': 420, 'US2': 291, 'US3': 481, 'US4': 238, 'US5': 92,
    'US6': 394, 'US7': 390, 'US8': 550, 'US8A': 147, 'US9': 90,
    'US10': 312, 'US11': 161, 'US12': 88,
    'RL1': 263, 'RL2': 775, 'RS1': 1298
}

def clean_description(raw_desc):
    """
    Cleans raw FASTA header descriptions by stripping:
    - Bracketed organism tags [Organism name]
    - UniProt prefixes: RecName: Full=, AltName: Full=, Short=..., Flags: ...
    - Redundant whitespace and punctuation.
    """
    d = raw_desc
    # Strip bracketed organism names: [Human alphaherpesvirus 1]
    d = re.sub(r'\[.*?\]', '', d)
    # Strip UniProt annotation prefixes
    d = re.sub(r'(?:RecName|AltName):\s*(?:Full=)?', '', d)
    d = re.sub(r'Flags:\s*[^;]+;?', '', d)
    d = re.sub(r'Short=[^;]+;?', '', d)
    # Strip excess whitespace and semicolons
    d = re.sub(r'\s+', ' ', d).strip(' ;')
    return d if d else "unlabeled description"

def parse_header(header_clean):
    """
    Extracts structured fields from diverse FASTA header formats.
    """
    # 1. Bracketed organism detection
    m_brackets = re.findall(r'\[([^\]]+)\]', header_clean)
    organism_as_written = m_brackets[-1].strip() if m_brackets else "unlabeled"
    
    # Defaults
    accession = "none"
    entry_name = "none"
    raw_desc = header_clean
    
    # 2. UniProt format: >sp|P04413.1|US03_HHV11 ...
    if header_clean.startswith('sp|') or header_clean.startswith('tr|'):
        parts = header_clean.split('|')
        if len(parts) >= 3:
            accession = parts[1]
            entry_parts = parts[2].split(' ', 1)
            entry_name = entry_parts[0]
            raw_desc = entry_parts[1] if len(entry_parts) > 1 else ""
    # 3. PDB format: >pdb|8V1T|A Chain A, ...
    elif header_clean.startswith('pdb|'):
        parts = header_clean.split('|')
        if len(parts) >= 3:
            pdb_id = parts[1]
            chain_parts = parts[2].split(' ', 1)
            chain_id = chain_parts[0]
            accession = f"{pdb_id}_{chain_id}"
            raw_desc = chain_parts[1] if len(chain_parts) > 1 else f"Chain {chain_id}"
    # 4. PRF format: >prf||1813262A ...
    elif header_clean.startswith('prf|'):
        parts = header_clean.split('|')
        if len(parts) >= 3:
            accession = parts[2].split(' ')[0] if parts[2] else (parts[1] if parts[1] else "none")
            raw_desc = header_clean.split(' ', 1)[1] if ' ' in header_clean else ""
    # 5. PIR format: >pir||VGBEDZ ...
    elif header_clean.startswith('pir|'):
        parts = header_clean.split('|')
        if len(parts) >= 3:
            accession = parts[2].split(' ')[0] if parts[2] else (parts[1] if parts[1] else "none")
            raw_desc = header_clean.split(' ', 1)[1] if ' ' in header_clean else ""
    # 6. Standard GenBank / RefSeq / Patent / NCBI format
    else:
        tokens = header_clean.split(' ', 1)
        accession = tokens[0]
        raw_desc = tokens[1] if len(tokens) > 1 else ""

    description_clean = clean_description(raw_desc)
    
    # 3. Organism inference and evidence matching
    hl = header_clean.lower()
    org_l = organism_as_written.lower()
    
    organism_inferred = "unknown"
    evidence_for_inference = "none"
    
    # Check HSV-1 keywords
    if "human alphaherpesvirus 1" in org_l:
        organism_inferred = "HSV-1"
        evidence_for_inference = "Human alphaherpesvirus 1"
    elif "human herpesvirus 1" in org_l or "human herpesvirus 1" in hl:
        organism_inferred = "HSV-1"
        evidence_for_inference = "human herpesvirus 1"
    elif any(k in entry_name for k in ["HHV11", "HHV1F", "HHV1E", "HHV1K", "HHV1M"]):
        organism_inferred = "HSV-1"
        evidence_for_inference = [k for k in ["HHV11", "HHV1F", "HHV1E", "HHV1K", "HHV1M"] if k in entry_name][0]
    elif "hsv-1" in hl or "hsv1" in hl:
        organism_inferred = "HSV-1"
        evidence_for_inference = "HSV-1"
    elif "strain 17" in hl or "strain hzt" in hl:
        organism_inferred = "HSV-1"
        evidence_for_inference = "strain 17" if "strain 17" in hl else "strain Hzt"
        
    # Check HSV-2 keywords
    elif "human alphaherpesvirus 2" in org_l:
        organism_inferred = "HSV-2"
        evidence_for_inference = "Human alphaherpesvirus 2"
    elif "human herpesvirus 2" in org_l or "human herpesvirus 2" in hl:
        organism_inferred = "HSV-2"
        evidence_for_inference = "human herpesvirus 2"
    elif any(k in entry_name for k in ["HHV2", "HHV22", "HHV23", "HHV2H", "HSV2S"]):
        organism_inferred = "HSV-2"
        evidence_for_inference = [k for k in ["HHV2", "HHV22", "HHV23", "HHV2H", "HSV2S"] if k in entry_name][0]
    elif "hsv-2" in hl or "hsv2" in hl or "type ii" in hl:
        organism_inferred = "HSV-2"
        evidence_for_inference = "HSV-2" if "hsv" in hl else "type II"
    elif any(k in org_l for k in ["strain 186", "strain g", "strain ss01", "strain sn03", "strain st04", "strain hg52"]):
        organism_inferred = "HSV-2"
        evidence_for_inference = "HSV-2 strain annotation"
        
    # Check other herpesviruses
    elif any(k in org_l for k in ["betaherpesvirus", "bovine alphaherpesvirus", "gallid alphaherpesvirus", "anatid alphaherpesvirus", "equid alphaherpesvirus", "macacine alphaherpesvirus"]) or any(k in entry_name for k in ["EHV4", "GAHVR"]):
        organism_inferred = "other herpesvirus"
        evidence_for_inference = organism_as_written if organism_as_written != "unlabeled" else entry_name
    elif "caudoviricetes" in org_l:
        organism_inferred = "other herpesvirus"
        evidence_for_inference = "Caudoviricetes sp."
        
    # Check non-viral host / vector / patent / HIV
    elif "homo sapiens" in org_l or "human" in entry_name:
        organism_inferred = "non-viral"
        evidence_for_inference = "Homo sapiens"
    elif "mus musculus" in org_l or "mouse" in entry_name:
        organism_inferred = "non-viral"
        evidence_for_inference = "Mus musculus"
    elif any(k in org_l for k in ["rattus norvegicus", "danio rerio", "pongo abelii", "nomascus leucogenys", "sphingobium"]):
        organism_inferred = "non-viral"
        evidence_for_inference = organism_as_written
    elif "patent" in hl:
        organism_inferred = "non-viral"
        evidence_for_inference = "patent"
    elif any(k in org_l for k in ["expression vector", "cloning vector", "pmttrep"]):
        organism_inferred = "non-viral"
        evidence_for_inference = organism_as_written
    elif "human immunodeficiency virus" in org_l or "gag protein" in hl:
        organism_inferred = "non-viral" # Retrovirus / non-herpesvirus
        evidence_for_inference = "Human immunodeficiency virus 1" if "human immunodeficiency virus" in org_l else "HIV-1 gag"

    return accession, entry_name, description_clean, organism_as_written, organism_inferred, evidence_for_inference

def guess_gene(header, description):
    """
    Assigns canonical HSV gene designation (UL1-UL56, US1-US12, RL1, RL2, RS1, etc.)
    """
    hl = (header + " " + description).lower()
    
    # 1. Direct UL / US / RL / RS regex matching
    m = re.search(r'\bul[-_]?([0-9]{1,2}[a-z]?)\b', hl)
    if m:
        val = m.group(1).upper()
        return f"UL{val}"
    m = re.search(r'\bus[-_]?([0-9]{1,2}[a-z]?)\b', hl)
    if m:
        val = m.group(1).upper()
        return f"US{val}"
    if re.search(r'\brl[-_]?1\b', hl) or 'icp34.5' in hl or 'neurovirulence' in hl:
        return 'RL1'
    if re.search(r'\brl[-_]?2\b', hl) or 'icp0' in hl or 'ie110' in hl or 'icpo' in hl or 'ring finger domain protein' in hl:
        return 'RL2'
    if re.search(r'\brs[-_]?1\b', hl) or 'icp4' in hl or 'ie175' in hl:
        return 'RS1'
        
    # 2. Glycoproteins
    if 'glycoprotein b' in hl or re.search(r'\bgb\b', hl): return 'UL27'
    if 'glycoprotein c' in hl or re.search(r'\bgc\b', hl): return 'UL44'
    if 'glycoprotein d' in hl or re.search(r'\bgd\b', hl): return 'US6'
    if 'glycoprotein e' in hl or re.search(r'\bge\b', hl): return 'US8'
    if 'glycoprotein g' in hl or re.search(r'\bgg\b', hl): return 'US4'
    if 'glycoprotein h' in hl or re.search(r'\bgh\b', hl): return 'UL22'
    if 'glycoprotein i' in hl or re.search(r'\bgi\b', hl): return 'US7'
    if 'glycoprotein k' in hl or re.search(r'\bgk\b', hl): return 'UL53'
    if 'glycoprotein l' in hl or re.search(r'\bgl\b', hl): return 'UL1'
    if 'glycoprotein m' in hl or re.search(r'\bgm\b', hl): return 'UL10'
    if 'glycoprotein n' in hl or re.search(r'\bgn\b', hl): return 'UL49A'
    if 'glycoprotein j' in hl or re.search(r'\bgj\b', hl): return 'US5'
    
    # 3. Core enzymes and structural proteins
    if 'thymidine kinase' in hl or re.search(r'\btk\b', hl): return 'UL23'
    if 'processivity factor' in hl or 'dna polymerase processivity' in hl: return 'UL42'
    if 'dna polymerase' in hl: return 'UL30'
    if 'major capsid' in hl or re.search(r'\bvp5\b', hl): return 'UL19'
    if 'small capsomer' in hl or 'small capsid' in hl or re.search(r'\bvp26\b', hl): return 'UL35'
    if 'capsid triplex subunit 1' in hl or re.search(r'\bvp19c\b', hl): return 'UL38'
    if 'capsid triplex subunit 2' in hl or re.search(r'\bvp23\b', hl): return 'UL18'
    if 'portal' in hl: return 'UL6'
    if 'capsid maturation protease' in hl or re.search(r'\bvp24\b', hl): return 'UL26'
    if 'scaffolding protein' in hl or 'scaffold protein' in hl or re.search(r'\bvp22a\b', hl): return 'UL26.5'
    if 'large tegument' in hl or 'deneddylase' in hl: return 'UL36'
    if 'vp11/12' in hl or 'vp11-12' in hl or 'vp11|12' in hl: return 'UL46'
    if 'vp13/14' in hl or 'vp13-14' in hl or 'vp13|14' in hl: return 'UL47'
    if 'vp16' in hl or 'vmw65' in hl or 'alpha-tif' in hl or 'transactivating tegument' in hl or 'alpha trans-inducing' in hl: return 'UL48'
    if 'vp22' in hl or 'viral protein 22' in hl: return 'UL49'
    if 'tegument host shutoff' in hl or re.search(r'\bvhs\b', hl): return 'UL41'
    if 'origin-binding' in hl: return 'UL9'
    if 'helicase-primase primase' in hl or re.search(r'\bprimase\b', hl): return 'UL52'
    if 'helicase-primase helicase' in hl or re.search(r'\bhelicase\b', hl): return 'UL5'
    if 'helicase-primase' in hl: return 'UL8'
    if 'single-stranded dna-binding' in hl or 'major dna-binding' in hl or 'icp8' in hl: return 'UL29'
    if 'uracil-dna glycosylase' in hl or 'udg' in hl or 'ung' in hl: return 'UL2'
    if 'deoxyuridine' in hl or 'dutpase' in hl: return 'UL50'
    if 'alkaline nuclease' in hl or 'deoxyribonuclease' in hl: return 'UL12'
    if 'terminase subunit 1' in hl: return 'UL28'
    if 'terminase subunit 2' in hl or 'large terminase' in hl: return 'UL15'
    if 'terminase subunit 3' in hl: return 'UL33'
    if 'nuclear egress lamina' in hl: return 'UL31'
    if 'nuclear egress membrane' in hl: return 'UL34'
    if 'inner tegument' in hl or 'capsid assembly' in hl: return 'UL37'
    if 'virion-packaging protein ul17' in hl or 'virion packaging protein ul17' in hl: return 'UL17'
    if 'virion-packaging protein' in hl or 'virion packaging' in hl: return 'UL17'
    if 'virion packaging protein ul25' in hl: return 'UL25'
    if 'ribonucleoside-diphosphate reductase large' in hl or 'ribonucleotide reductase subunit 1' in hl or 'icp6' in hl: return 'UL39'
    if 'ribonucleoside-diphosphate reductase small' in hl or 'ribonucleotide reductase subunit 2' in hl: return 'UL40'
    if 'ribonucleotide reductase' in hl: return 'UL39'
    if 'icp27' in hl or 'ie63' in hl or 'mrna export factor' in hl or 'multifunctional expression regulator' in hl: return 'UL54'
    if 'icp22' in hl or 'ie68' in hl or 'regulatory protein icp 22' in hl: return 'US1'
    if 'icp47' in hl or 'ie12' in hl: return 'US12'
    if 'myristylated tegument' in hl: return 'UL11'
    if 'tegument serine/threonine protein kinase' in hl: return 'US3'
    
    # Non-HSV specific genes
    if 'gag protein' in hl or 'gag' in hl: return 'HIV-1 gag'
    if 'pol protein' in hl: return 'HIV-1 pol'
    if 't cell receptor' in hl or 't-cell receptor' in hl: return 'TCR'
    
    return 'unknown'

def check_contamination(header, description, organism_as_written, organism_inferred):
    """
    Evaluates whether a record constitutes non-HSV contamination.
    """
    hl = (header + " " + description).lower()
    org_l = organism_as_written.lower()
    
    # 1. Patent
    if 'patent' in hl:
        return True, "patent"
    # 2. TCR
    if 't cell receptor' in hl or 't-cell receptor' in hl or re.search(r'\btcr\b', hl):
        return True, "TCR"
    # 3. Antibody / Nanobody
    if any(k in hl for k in ['antibody', 'heavy chain', 'light chain', 'nanobody', 'nb1_gbhsv', 'chain j', 'chain n', 'chain k', 'chain l', 'chain m', 'chain o']):
        if 'envelope' not in hl and 'capsid' not in hl and 'tegument' not in hl:
            return True, "antibody"
    # 4. Expression Vectors & Cloning Vectors
    if any(k in org_l for k in ['expression vector', 'cloning vector', 'pmttrep', 'bac cloning']):
        return True, "vector"
    # 5. Non-viral host proteins / factors
    if any(k in org_l for k in ['homo sapiens', 'mus musculus', 'rattus norvegicus', 'danio rerio', 'pongo', 'nomascus', 'sphingobium']):
        return True, "host"
    if any(k in hl for k in ['complement c3', 'poliovirus receptor', 'herpesvirus entry mediator', 'heparan sulfate', 'ubiquitin carboxyl-terminal hydrolase 7', 'importin subunit alpha-1', 'paired immunoglobulin-like type 2 receptor']):
        return True, "host"
    # 6. Other viruses (CMV, non-HSV alphaherpesviruses, HIV gag/pol)
    if organism_inferred == "other herpesvirus":
        return True, "other-virus"
    if 'human immunodeficiency virus' in org_l or 'gag protein' in hl or 'pol protein' in hl or 'caudoviricetes' in org_l:
        return True, "other-virus"
        
    return False, "none"

def load_fasta_file(path, source_file_name):
    """
    Reads FASTA file into list of dicts.
    """
    records = []
    print(f"Loading {path}...")
    with open(path, 'r', encoding='utf-8', errors='ignore') as f:
        curr_h = None
        curr_s = []
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith('>'):
                if curr_h is not None:
                    records.append({
                        'source_file': source_file_name,
                        'full_header': curr_h,
                        'seq': ''.join(curr_s).upper()
                    })
                curr_h = line_str[1:]
                curr_s = []
            else:
                curr_s.append(line_str)
        if curr_h is not None:
            records.append({
                'source_file': source_file_name,
                'full_header': curr_h,
                'seq': ''.join(curr_s).upper()
            })
    print(f"Loaded {len(records)} records from {source_file_name}")
    return records

def run_curation():
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
    DATA_DIR = os.path.join(REPO_ROOT, "data")
    CURATION_DIR = SCRIPT_DIR
    BENCHMARKING_DIR = os.path.join(REPO_ROOT, "benchmarking")
    
    # 1. Load both FASTA files
    f1_name = 'HSV2_non_redundant (2).fasta'
    f2_name = 'non_redundant.fasta'
    f1 = os.path.join(DATA_DIR, f1_name)
    f2 = os.path.join(DATA_DIR, f2_name)
    
    raw_records = load_fasta_file(f1, f1_name) + load_fasta_file(f2, f2_name)
    total_records = len(raw_records)
    print(f"Total FASTA records to curate: {total_records}")
    
    # 2. Sequential record indexing & metadata extraction
    curated = []
    seq_hash_map = {}
    exact_dup_groups = defaultdict(list)
    
    for idx, r in enumerate(raw_records, 1):
        full_h = r['full_header']
        seq = r['seq']
        src = r['source_file']
        
        # Parse fields
        accession, entry_name, desc_clean, org_as_written, org_inferred, evidence = parse_header(full_h)
        gene_guess = guess_gene(full_h, desc_clean)
        
        # Sequence statistics
        seq_len = len(seq)
        ambig_count = sum(seq.count(c) for c in ['X', 'B', 'Z', 'U', '*'])
        has_stop = "yes" if '*' in seq else "no"
        
        # Exact duplicate hash (ignoring trailing/internal stop if stripped)
        seq_norm = seq.replace('*', '')
        seq_hash = hashlib.sha256(seq_norm.encode('ascii')).hexdigest()
        
        if seq_hash in seq_hash_map:
            first_idx = seq_hash_map[seq_hash]
            exact_duplicate_of = str(first_idx)
            exact_dup_groups[first_idx].append(idx)
        else:
            seq_hash_map[seq_hash] = idx
            exact_duplicate_of = "none"
            
        # Contamination check
        is_contam, contam_reason = check_contamination(full_h, desc_clean, org_as_written, org_inferred)
        contamination_flag = "yes" if is_contam else "no"
        
        # Partial check
        # Header explicit keywords
        hl = full_h.lower()
        is_partial = "no"
        if any(k in hl for k in ['partial', 'fragment', 'fragments', 'truncated', 'short', 'peptide']):
            is_partial = "yes"
        elif gene_guess in CANONICAL_LENGTHS:
            exp_len = CANONICAL_LENGTHS[gene_guess]
            if seq_len < 0.65 * exp_len:
                is_partial = "yes"
        elif seq_len < 50:
            is_partial = "yes"
            
        curated.append({
            'record_index': idx,
            'source_file': src,
            'full_header': full_h,
            'accession': accession,
            'entry_name': entry_name,
            'description_clean': desc_clean,
            'organism_as_written': org_as_written,
            'organism_inferred': org_inferred,
            'evidence_for_inference': evidence,
            'gene_or_protein_guess': gene_guess,
            'sequence_length_aa': seq_len,
            'ambiguous_residues': ambig_count,
            'has_stop_codon': has_stop,
            'exact_duplicate_of': exact_duplicate_of,
            'near_duplicate_cluster': 'singleton', # updated during clustering
            'is_partial': is_partial,
            'contamination_flag': contamination_flag,
            'contamination_reason': contam_reason,
            'recommended_status': 'UNSET', # updated below
            'notes': '',
            'seq_norm': seq_norm
        })

    print("Phase 1 complete: Metadata and exact duplicates extracted.")
    
    # 3. Near-Duplicate Greedy Clustering (>90% identity)
    # Optimized greedy clustering:
    # 1. Partition sequences by gene and length bucket (L // 25).
    # 2. For identical lengths, calculate Hamming difference with early termination.
    # 3. For unequal lengths (difference <= 10%), use sliding window identity test.
    print("Phase 2: Running near-duplicate greedy clustering (>90% identity)...")
    
    def check_near_dup(s1, s2, ls, lr):
        if abs(lr - ls) / max(lr, ls) > 0.10:
            return False
        if lr == ls:
            max_d = int(0.10 * lr)
            diffs = 0
            for a, b in zip(s1, s2):
                if a != b:
                    diffs += 1
                    if diffs > max_d:
                        return False
            return True
        # Unequal lengths within 10%
        short_s, long_s = (s1, s2) if ls < lr else (s2, s1)
        len_short, len_long = len(short_s), len(long_s)
        shift_max = len_long - len_short
        for shift in range(shift_max + 1):
            m = sum(1 for a, b in zip(short_s, long_s[shift:shift+len_short]) if a == b)
            if (m / len_long) >= 0.90:
                return True
        return False

    # Group records by partition
    partitions = defaultdict(list)
    for r in curated:
        g = r['gene_or_protein_guess']
        l_bucket = r['sequence_length_aa'] // 25
        partitions[(g, l_bucket)].append(r)
        
    cluster_counter = 1
    near_dup_clusters_dict = defaultdict(list)
    
    for (gene, l_buck), items in partitions.items():
        part_clusters = [] # list of (rep_s, lr, [members])
        for r in items:
            s = r['seq_norm']
            ls = len(s)
            matched = False
            for rep_tuple in part_clusters:
                rep_s, lr, members = rep_tuple
                if check_near_dup(s, rep_s, ls, lr):
                    members.append(r)
                    matched = True
                    break
            if not matched:
                part_clusters.append((s, ls, [r]))
                
        # Assign cluster labels
        for rep_s, lr, members in part_clusters:
            if len(members) > 1:
                c_id = f"CLUST_{cluster_counter:04d}"
                cluster_counter += 1
                for mem in members:
                    mem['near_duplicate_cluster'] = c_id
                    near_dup_clusters_dict[c_id].append(mem['record_index'])
            else:
                members[0]['near_duplicate_cluster'] = "singleton"

    print(f"Phase 2 complete: Formed {len(near_dup_clusters_dict)} near-duplicate clusters (>90% identity).")

    # 4. Determine Recommended Status (KEEP / MERGE / EXCLUDE / REVIEW) & Notes
    # Rules:
    # - EXCLUDE: Any record with contamination_flag == "yes" (host, vector, antibody, TCR, patent, other-virus).
    # - REVIEW: Any record with organism_inferred == "unknown", or ambiguous unverified headers.
    # - KEEP vs MERGE:
    #   For each viral gene in HSV-1 and HSV-2:
    #   We choose the highest-quality, canonical, full-length representative record to KEEP.
    #   All other records in that gene's duplicate groups / near-duplicate clusters are set to MERGE.
    print("Phase 3: Selecting representative KEEP records and setting MERGE / EXCLUDE / REVIEW status...")
    
    # Pre-select best canonical representative for each (organism, gene)
    # Priority for representative:
    # 1. Non-partial
    # 2. Reviewed UniProt (starts with sp|)
    # 3. RefSeq / Strain 17 / HG52 reference
    # 4. Zero ambiguous residues
    # 5. Length closest to canonical length
    best_candidate = {} # (org, gene) -> record
    
    for r in curated:
        if r['contamination_flag'] == 'yes':
            r['recommended_status'] = 'EXCLUDE'
            r['notes'] = f"Excluded: {r['contamination_reason']} contamination"
            continue
            
        if r['organism_inferred'] == 'unknown':
            r['recommended_status'] = 'REVIEW'
            r['notes'] = "Review: Organism unconfirmed in header; manual verification required"
            continue
            
        if r['organism_inferred'] not in ['HSV-1', 'HSV-2']:
            r['recommended_status'] = 'EXCLUDE'
            r['notes'] = f"Excluded: non-target species ({r['organism_inferred']})"
            continue
            
        # Target viral record
        org = r['organism_inferred']
        gene = r['gene_or_protein_guess']
        
        # If gene is unknown, mark for curator review
        if gene == 'unknown':
            r['recommended_status'] = 'REVIEW'
            r['notes'] = "Review: Uncharacterized/hypothetical protein; requires manual curation"
            continue
        
        # Scoring record quality
        score = 0
        if r['is_partial'] == 'no': score += 1000
        if r['full_header'].startswith('sp|'): score += 500
        if 'strain 17' in r['full_header'].lower() or 'hg52' in r['full_header'].lower(): score += 300
        if r['ambiguous_residues'] == 0: score += 200
        else: score -= r['ambiguous_residues'] * 10
        
        if gene in CANONICAL_LENGTHS:
            diff = abs(r['sequence_length_aa'] - CANONICAL_LENGTHS[gene])
            score -= diff
        else:
            score += min(r['sequence_length_aa'], 500)
            
        key = (org, gene)
        if key not in best_candidate or score > best_candidate[key][0]:
            best_candidate[key] = (score, r)

    # Set KEEP and MERGE
    kept_records = set(rec['record_index'] for score, rec in best_candidate.values())
    
    for r in curated:
        idx = r['record_index']
        if r['recommended_status'] in ['EXCLUDE', 'REVIEW']:
            continue
            
        org = r['organism_inferred']
        gene = r['gene_or_protein_guess']
        
        if idx in kept_records:
            r['recommended_status'] = 'KEEP'
            r['notes'] = f"Selected representative for {org} {gene}"
        else:
            # If it is an exact duplicate
            if r['exact_duplicate_of'] != 'none':
                r['recommended_status'] = 'MERGE'
                r['notes'] = f"Merge: Exact sequence duplicate of record {r['exact_duplicate_of']}"
            else:
                r['recommended_status'] = 'MERGE'
                r['notes'] = f"Merge: Redundant clinical isolate/variant of {org} {gene}"

    print("Phase 3 complete: All 22,710 records assigned statuses.")
    
    # 5. Compute Summaries A through F
    counts = {
        'total_records': total_records,
        'records_per_source_file': dict(Counter(r['source_file'] for r in curated)),
        'records_per_organism_inferred': dict(Counter(r['organism_inferred'] for r in curated)),
        'records_per_recommended_status': dict(Counter(r['recommended_status'] for r in curated)),
        'exact_duplicate_groups_count': len(exact_dup_groups),
        'near_duplicate_clusters_count': len(near_dup_clusters_dict),
        'partial_records_count': sum(1 for r in curated if r['is_partial'] == 'yes'),
        'contamination_flagged_count': sum(1 for r in curated if r['contamination_flag'] == 'yes'),
        'contamination_by_reason': dict(Counter(r['contamination_reason'] for r in curated if r['contamination_flag'] == 'yes'))
    }
    
    # Save summary json
    summary_json_path = os.path.join(BENCHMARKING_DIR, 'curation_summary_stats.json')
    with open(summary_json_path, 'w', encoding='utf-8') as f:
        json.dump(counts, f, indent=2)
    print(f"Saved {summary_json_path}")
    
    # 6. Write Full Master TSV File (All 22,710 records)
    tsv_filename = os.path.join(CURATION_DIR, 'hsv_curation_records.tsv')
    print(f"Writing master TSV table to {tsv_filename}...")
    columns = [
        'record_index', 'source_file', 'full_header', 'accession', 'entry_name',
        'description_clean', 'organism_as_written', 'organism_inferred',
        'evidence_for_inference', 'gene_or_protein_guess', 'sequence_length_aa',
        'ambiguous_residues', 'has_stop_codon', 'exact_duplicate_of',
        'near_duplicate_cluster', 'is_partial', 'contamination_flag',
        'contamination_reason', 'recommended_status', 'notes'
    ]
    
    with open(tsv_filename, 'w', encoding='utf-8') as f:
        f.write('\t'.join(columns) + '\n')
        for r in curated:
            row = [str(r[c]) for c in columns]
            f.write('\t'.join(row) + '\n')
    print(f"Master TSV written ({os.path.getsize(tsv_filename):,} bytes).")
    
    # 7. Write Full Markdown Table File (All 22,710 records)
    md_table_filename = os.path.join(CURATION_DIR, 'hsv_curation_table.md')
    print(f"Writing full markdown table to {md_table_filename}...")
    with open(md_table_filename, 'w', encoding='utf-8') as f:
        f.write("# HSV-1 & HSV-2 Proteome Analysis — Complete Master Curation Table\n\n")
        f.write(f"Total curated records: **{total_records:,}**\n\n")
        f.write("| " + " | ".join(columns) + " |\n")
        f.write("| " + " | ".join(["---"] * len(columns)) + " |\n")
        for r in curated:
            # Escape pipe characters in text fields
            row = [str(r[c]).replace('|', '\\|') for c in columns]
            f.write("| " + " | ".join(row) + " |\n")
    print(f"Master Markdown table written ({os.path.getsize(md_table_filename):,} bytes).")

    # 8. Prepare Structured Lists for Sections B, C, D, E, F
    print("Preparing structured section datasets...")
    
    # Section B: KEEP LIST
    keep_list = [r for r in curated if r['recommended_status'] == 'KEEP']
    keep_list_hsv1 = [r for r in keep_list if r['organism_inferred'] == 'HSV-1']
    keep_list_hsv2 = [r for r in keep_list if r['organism_inferred'] == 'HSV-2']
    
    # Section C: MERGE LIST (Sample groups + all exact duplicates)
    merge_groups = []
    for first_idx, dup_indices in exact_dup_groups.items():
        rep = curated[first_idx - 1]
        dups = [curated[i - 1] for i in dup_indices]
        merge_groups.append({
            'rep_record': rep,
            'duplicate_records': dups,
            'reason': '100% exact identical sequence match'
        })
        
    # Section D: EXCLUDE LIST
    exclude_list = [r for r in curated if r['recommended_status'] == 'EXCLUDE']
    
    # Section E: COVERAGE MAP
    canonical_genes = [f"UL{i}" for i in range(1, 57)] + [f"US{i}" for i in range(1, 13)] + ['RL1', 'RL2', 'RS1']
    
    hsv1_kept_by_gene = {r['gene_or_protein_guess']: r for r in keep_list_hsv1}
    hsv2_kept_by_gene = {r['gene_or_protein_guess']: r for r in keep_list_hsv2}
    
    # Section F: DATA QUALITY WARNINGS
    warnings_ambig = [r for r in curated if r['sequence_length_aa'] > 0 and (r['ambiguous_residues'] / r['sequence_length_aa']) > 0.05]
    warnings_stop = [r for r in curated if r['has_stop_codon'] == 'yes']
    warnings_short = [r for r in curated if r['sequence_length_aa'] < 50]
    warnings_long = [r for r in curated if r['sequence_length_aa'] > 2000]
    warnings_unlabeled = [r for r in curated if r['organism_as_written'] == 'unlabeled']
    warnings_inconsistent = [r for r in curated if r['organism_as_written'] == 'Human alphaherpesvirus 1' and 'gag' in r['full_header'].lower()]

    summary_details = {
        'keep_hsv1_count': len(keep_list_hsv1),
        'keep_hsv2_count': len(keep_list_hsv2),
        'merge_count': sum(1 for r in curated if r['recommended_status'] == 'MERGE'),
        'exclude_count': len(exclude_list),
        'review_count': sum(1 for r in curated if r['recommended_status'] == 'REVIEW'),
        'warnings_ambig_count': len(warnings_ambig),
        'warnings_stop_count': len(warnings_stop),
        'warnings_short_count': len(warnings_short),
        'warnings_long_count': len(warnings_long),
        'warnings_unlabeled_count': len(warnings_unlabeled),
        'warnings_inconsistent_count': len(warnings_inconsistent)
    }
    
    details_json_path = os.path.join(BENCHMARKING_DIR, 'curation_details_stats.json')
    with open(details_json_path, 'w', encoding='utf-8') as f:
        json.dump(summary_details, f, indent=2)
    print(f"Saved {details_json_path}")
    print("All curation steps completed successfully!")

if __name__ == '__main__':
    run_curation()
