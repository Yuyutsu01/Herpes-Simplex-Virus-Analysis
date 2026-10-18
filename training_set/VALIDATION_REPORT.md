# Comprehensive Validation Report: HSV-1 & HSV-2 ML Proteomics Dataset
## Final Multi-Tiered Verification & Readiness Assessment

This document details the final validation results for the four preparatory deliverables of the Herpes Simplex Virus Machine Learning Proteomics study. All evaluations were executed reproducibly against local inputs and cached NCBI RefSeq genomes (`NC_001806.2` and `NC_001798.2`).

---
### 1. Results of Deliverable 1: RefSeq Genomic Length Verification

Both flagged suspicious lengths were globally aligned against canonical CDS translations extracted from GenBank RefSeq reference strain records:

#### Target 1: HSV-2 UL31 (Kept Accession: `WZW62192.1`, Length: 306 aa)
- **Verdict:** **`CONFIRMED_FULL_LENGTH`**
- **Evaluation:** Global pairwise alignment against the canonical HSV-2 HG52 RefSeq CDS (`YP_009137183.1`, 305 aa) confirms authentic full-length architecture with **91.15% sequence identity** (278/305 matching residues) and **100.0% query coverage**. The sequence initiates precisely at Met1 (`MYDTDPHRRG...`), terminates at the canonical C-terminal stop codon (`...YATFDEFPPP*`), and possesses zero internal stop codons. The 1-amino acid difference reflects a natural, strain-specific insertion (`SAAGGTLGVVRR` vs `SAAGNP-GVERR`) in the flexible N-terminal basic domain. Reference genome cross-examination confirms that canonical HSV-1 UL31 is also 306 aa (`YP_009137106.1`). The 231-aa entry in the initial HSV-1 set (`SBO07674.1`) was an N-terminally truncated GenBank submission starting at residue 76, confirming `WZW62192.1` as the authentic full-length protein.

#### Target 2: HSV-2 UL27 (Kept Accession: `QAU10510.1`, Length: 904 aa)
- **Verdict:** **`CONFIRMED_FULL_LENGTH`**
- **Evaluation:** Global pairwise alignment against the canonical HSV-2 HG52 RefSeq CDS (`YP_009137179.1`, 904 aa) demonstrates **100.0% sequence identity** across all 904 amino acids (904/904 matches, 0 mismatches, 0 gaps). The reading frame initiates at Met1 (`MRQGAPAR...`), terminates at the canonical stop (`...VWERR*`), and contains zero internal stops or sequencing ambiguities. Both HSV-1 and HSV-2 canonical gB proteins are exactly 904 aa (`YP_009137102.1` and `YP_009137179.1`). The shorter 742-aa length in the initial HSV-1 record (`SBO07646.1`) was an N-terminally truncated fragment lacking the signal sequence, confirming `QAU10510.1` as the complete, intact 904-aa glycoprotein B precursor.

---
### 2. Row Counts and Class Balance Table (`hsv_annotation.tsv`)

- **Total Rows:** **149** (75 HSV-1 + 74 HSV-2)
- **Gene Redundancy:** 0 duplicates; every gene appears exactly once per species

| Biological Axis | Class Label | HSV-1 Count | HSV-2 Count | Total | Class Proportion |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **Temporal** | `IE` | 5 | 5 | 10 | 6.7% |
| **Temporal** | `E` | 36 | 35 | 71 | 47.7% |
| **Temporal** | `E-L` | 1 | 1 | 2 | 1.3% |
| **Temporal** | `L` | 33 | 33 | 66 | 44.3% |
| **Functional** | `Capsid_Assembly` | 13 | 13 | 26 | 17.4% |
| **Functional** | `Egress` | 2 | 2 | 4 | 2.7% |
| **Functional** | `Glycoprotein_Envelope` | 17 | 17 | 34 | 22.8% |
| **Functional** | `Immune_Evasion` | 2 | 2 | 4 | 2.7% |
| **Functional** | `Replication` | 14 | 13 | 27 | 18.1% |
| **Functional** | `Regulatory` | 4 | 4 | 8 | 5.4% |
| **Functional** | `Tegument` | 23 | 23 | 46 | 30.9% |
| **Membrane** | `yes` | 20 | 20 | 40 | 26.8% |
| **Membrane** | `no` | 55 | 54 | 109 | 73.2% |

---
### 3. FASTA Record Counts and Species Breakdown (`hsv1_hsv2_curated_reference.fasta`)

- **Total FASTA Records Exported:** **149**
- **HSV-1 Records:** **75** records
- **HSV-2 Records:** **74** records
- **Line Wrapping:** Standard 80 characters per line
- **Alphabet Compliance:** 100% standard 20-amino-acid alphabet (zero ambiguous characters: no X, B, Z, U, *, or -)
- **Cryptographic Reproducibility:** Verified; all 149 SHA-256 hashes generated in `hsv1_hsv2_curated_reference.manifest.tsv` match bit-for-bit across runs.

---
### 4. Discrepancy & Concordance Analysis

- **KEEP List Concordance:** 100% concordance between `hsv_curation_report.md` KEEP lists and `hsv_curation_records.tsv`.
- **Source FASTA Presence:** Exactly 149 of 149 target accessions were located in `non_redundant.fasta` and `HSV2_non_redundant (2).fasta`. Zero accessions missing.
- **Duplicate Occurrences in Source Files:** **0 redundant occurrences**; all 149 accessions appear uniquely in their source FASTA files.
- **Species Asymmetry:** Exactly one asymmetric gene detected: **UL70** (Protein ID: `SCL76894`) in HSV-1, representing a legacy GenBank strain 17 annotation for DNA primase (homolog of UL52); absent in HSV-2.

---
### 5. Ready for Machine Learning Checklist

| Evaluation Criterion | Requirement | Observed Status | ML Readiness Verdict |
| :--- | :--- | :--- | :---: |
| `verify_lengths.md` verdicts | No UNVERIFIABLE verdicts | 2 CONFIRMED_FULL_LENGTH | **PASS** |
| `hsv_annotation.tsv` row count | Exactly 149 rows | 149 rows | **PASS** |
| Gene singularity per species | Each gene present exactly once per species | 0 duplicate genes | **PASS** |
| FASTA record count | Exactly 149 records | 149 records | **PASS** |
| Standard 20-aa alphabet | No non-standard residues (X, B, Z, U, *, -) | 100% standard alphabet | **PASS** |
| Curation status resolution | All protein_ids resolve to status == KEEP | 149/149 resolved | **PASS** |
| Species-gene uniqueness | No duplicate (species, gene) pairs | 0 duplicate pairs | **PASS** |

---
### 6. Final Study Determination

**READY FOR FEATURE EXTRACTION: YES**
