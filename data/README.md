# HK-SARS-CoV-2-genomic-epidemiology-*data*

---
This directory hosts the data for different analysis of the study. Specifically, The phylodynamics analysis, BDSS estimation and transmission bottleneck estimation. 

## Description:

### seq_local_analysis/hk_case_primer_masked_2021-05-08.afa: 
	The HK sequences for local analysis. 
	Local masking strategy: The head and tail 100 nt are masked. Ten other 10 error-prone primer binding sites are also masked. 

### seq_global_analysis/hk_case_primer_masked_GISAID_masked_2021-05-08.afa: 
	The HK sequences for global analysis (harmonized masking strategy with GISIAD data).
	Global masking strategy: In addition to the local masking strategy, another 73 masking sites applied in GISIAD data were also masked.

### seq_global_analysis/seq_global_all.fasta:
	The global reference sequences used for global analysis (n = 2301).
	It includes two parts:
		1. The most-similar-to-HK global sequences from GISAID (p-distance <= 3 to at least one HK sequence, n = 1022). (The first 1022 sequences)
		2. The representative sequences corresponding to the 1279 PANGO lineages (n = 1279). (the last 1279 sequences)
	The sequences are already aligned to the reference genome, and the masking strategy is consistent with the HK sequences for global analysis, no further alignment is needed.

### seq_global_analysis/seq_global_all_rmdup.fasta:
	Removed 6 duplicated sequences in the "seq_global_analysis/seq_global_all.fasta". (n = 2295)

### seq_global_analysis/metadata_global_top3hit.csv:
	The corresponding metadata for the 385 GISAID most-similar-to-HK sequences (top 3 hits of each HK sequence).

### seq_global_analysis/metadata_pango_ref_1279.csv:
	The corresponding metadata for the above GISAID most-similar-to-HK sequences.

### metadata_2021-04-30.xlsx: 
	Metadata for all HK cases up to 2021-01-26, the new column "reclassification" replaced the previous "Classification" for more accurate importation information.

### PANGO_lineage_2021-05-08.csv: 
	The PANGO lineage assignment for all the HK sequences using PANGOLIN (v 2.3.9 with pangoLEARN 2021-04-23). Please use this information when referring to PANGO lineages.

### nextclade_lineage_2021-05-08.csv: 
	Nextclade lineage assignment for all the HK sequences.

---