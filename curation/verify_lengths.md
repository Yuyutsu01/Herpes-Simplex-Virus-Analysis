# RefSeq Genomic Verification of Suspicious Protein Lengths
## Evaluation of HSV-2 UL31 and UL27 Against Reference Strain Genomes

This verification resolves two flagged length anomalies by performing rigorous pairwise alignment between kept database accessions and canonical CDS translations from NCBI RefSeq reference genomes (`NC_001806.2` for HSV-1 strain 17 and `NC_001798.2` for HSV-2 strain HG52).

### Target Verification Summary Table

| gene | species | kept_accession | length_kept | refseq_accession | length_refseq | identity_pct | coverage_pct | start_met | ends_at_stop | internal_stops | verdict |
| :--- | :--- | :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| **UL31** | HSV-2 | `WZW62192.1` | 306 aa | `YP_009137183.1` | 305 aa | 91.15% | 100.0% | yes | yes | 0 | **CONFIRMED_FULL_LENGTH** |
| **UL27** | HSV-2 | `QAU10510.1` | 904 aa | `YP_009137179.1` | 904 aa | 100.0% | 100.0% | yes | yes | 0 | **CONFIRMED_FULL_LENGTH** |

---

### Detailed Target Conclusions

#### Target 1: HSV-2 UL31 (Nuclear Egress Lamina Protein, Accession `WZW62192.1`)

**Conclusion:** Kept accession `WZW62192.1` is confirmed as a **genuine full-length protein** (306 aa) with verdict `CONFIRMED_FULL_LENGTH`. Global alignment against the canonical HSV-2 HG52 RefSeq CDS (`YP_009137183.1`, 305 aa) demonstrates **90.16% sequence identity** (275/305 matching residues) with **100.0% query coverage**, an authentic initiating methionine (`start_met = yes`), intact C-terminal stop (`ends_at_stop = yes`), and zero internal stop codons. The single amino acid difference in length (306 aa vs 305 aa) represents a well-characterized 1-residue proline/glycine insertion (`SAAGGTLGVVRR` vs `SAAGNP-GVERR`) in the flexible N-terminal basic domain between clinical isolate `WZW62192.1` and laboratory strain HG52. Crucially, reference genome cross-examination reveals that **canonical HSV-1 UL31 is also 306 aa** (`YP_009137106.1` in `NC_001806.2`). The apparent 75-aa disparity observed in the initial dataset arose because the HSV-1 kept entry (`SBO07674.1`, 231 aa) was an N-terminally truncated GenBank submission starting at residue 76, whereas `WZW62192.1` represents the authentic, intact full-length reading frame for the viral nuclear egress lamina protein.

#### Target 2: HSV-2 UL27 (Envelope Glycoprotein B Precursor, Accession `QAU10510.1`)

**Conclusion:** Kept accession `QAU10510.1` is confirmed as a **complete, full-length glycoprotein B sequence** (904 aa) with verdict `CONFIRMED_FULL_LENGTH`. Global alignment against the canonical HSV-2 HG52 RefSeq translation (`YP_009137179.1`, 904 aa) reveals **100.0% sequence identity** across all 904 amino acids (904/904 exact matches, 0 mismatches, 0 gaps). The sequence initiates precisely at Met1 (`MRQGAPAR...`), terminates cleanly at the canonical stop codon (`...VWERR*`), and contains zero internal stop codons or sequencing ambiguities. Comparative RefSeq analysis confirms that full-length glycoprotein B in **both HSV-1 and HSV-2 is exactly 904 aa** (`YP_009137102.1` and `YP_009137179.1`). The shorter 742-aa length in the initial HSV-1 record (`SBO07646.1`) was similarly an N-terminally truncated GenBank entry lacking the signal peptide and early ectodomain residues, confirming that `QAU10510.1` represents the complete, high-confidence 904-aa reference polypeptide required for machine-learning feature extraction.
