# Ssd1_CRAC_demult_dedup_20190114_all

This directory contains selected output of the CRAC pipeline.
We only included output that was directly needed for downstream analysis, e.g. read counts, pileups.
We did *not* include large intermediate files, such as trimmed fastq or aligned bam files.

## bedgraph_genomecov

Bedgraph files showing counts of (deduplicated) aligned reads on the whole genome, separately for minus and plus strand reads.

* 20190114_Ssd1_CRAC_trimmed_NNNAGAGC_SSD1_3_42_minus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNAGAGC_SSD1_3_42_plus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNCTAGC_SSD1_4_42_minus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNCTAGC_SSD1_4_42_plus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNGACTTAGC_BY4741_minus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNGACTTAGC_BY4741_plus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30_minus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30_plus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNTGGAGC_SSD1_4_30_minus.bedgraph
* 20190114_Ssd1_CRAC_trimmed_NNNTGGAGC_SSD1_4_30_plus.bedgraph

## multicov_analyses

Tab-separated gff-like text file with stranded counts of (deduplicated) aligned reads with on each full-length transcript, and transcript co-ordinates.

* allsample_transcriptcounts.txt

## pyCalculateFDRs_analyses

Tab-separated gff-like text files showing peak locations, heights, and false discovery rate.

* 20190114_Ssd1_CRAC_trimmed_NNNAGAGC_SSD1_3_42_output_FDRs.gtf
* 20190114_Ssd1_CRAC_trimmed_NNNCTAGC_SSD1_4_42_output_FDRs.gtf
* 20190114_Ssd1_CRAC_trimmed_NNNGACTTAGC_BY4741_output_FDRs.gtf
* 20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30_output_FDRs.gtf
* 20190114_Ssd1_CRAC_trimmed_NNNTGGAGC_SSD1_4_30_output_FDRs.gtf

## pyPileup_analyses

Tab-separated text files showing "pileups" of read counts with deletions and mutations, and sequence, for specific Ssd1-targeted genes only. 

* 20190114_Ssd1_CRAC_trimmed_NNNAGAGC_SSD1_3_42_pileups.txt
* 20190114_Ssd1_CRAC_trimmed_NNNCTAGC_SSD1_4_42_pileups.txt
* 20190114_Ssd1_CRAC_trimmed_NNNGACTTAGC_BY4741_pileups.txt
* 20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30_pileups.txt
* 20190114_Ssd1_CRAC_trimmed_NNNTGGAGC_SSD1_4_30_pileups.txt
