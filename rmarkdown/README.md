# rmarkdown

Analysis scripts in rmarkdown (Rmd) format, that produce many figures in the manuscript.
These mostly run *after* the CRAC pipeline has completed.


## process_RNAseq.Rmd

Preprocesses RNA-seq data measuring mRNA levels in conditions matched to CRAC experiments. These data are from Bresson et al, 2020, and archived in GEO: GSE148166.

The script takes fastq files as input, and ends with a text-format coverage file on RNA transcripts: `results/RNAseq_stressmatched_transcriptcounts.txt`.

It relies on genome sequences and transcript annotations from `input_annotation`, and uses the same novoalign index as `src/CRAC_pipeline_SE_demult_dedup.py`.

This script runs completely independently of the CRAC pipeline.
It needs to be run before `normalise_Ssd1_CRAC_counts_vs_RNAseq.Rmd`, which relies on the transcript counts that this script outputs.


## normalise_Ssd1_CRAC_counts_vs_RNAseq.Rmd

Analyse relative *enrichment* of Ssd1 on mRNA, by processing time & condition matched RNA seq and Ssd1 CRAC data. Then we crudely normalize Ssd1 *counts* against RNA-seq on full-length transcripts.

This script relies on CRAC transcript counts produced by `multiBamCov` running in the CRAC pipeline, and found in `Ssd1_CRAC_demult_dedup_20190114_all/multicov_analyses/allsample_transcriptcounts.txt`. So it must be run after  `src/CRAC_pipeline_SE_demult_dedup.py`, see the repository `README.md` for details. 

This script also relies on RNA-seq transcript counts produced by `multiBamCov` running `process_RNAseq.Rmd`, found in `results/RNAseq_stressmatched_transcriptcounts.txt`. It must be run after `process_RNAseq.Rmd`.


## compare_target_transcripts.Rmd

Takes the list of transcripts enriched in Ssd1 binding in the CRAC dataset, and compares them with Ssd1-targeted transcripts from previous publications that do RIP-ChIP or RIP-seq:

* Hogan 2008 
* Jansen 2009
* Hose 2020

These publicly available data need to be downloaded to the `input_targets` directory, see `input_targets/README.md` for details.

This script also relies on lists of Ssd1-enriched transcripts produced by `normalise_Ssd1_CRAC_counts_vs_RNAseq.Rmd`, found in `results/Ssd1_targets_TPMratio4x_30C.txt`
and `results/Ssd1_targets_TPMratio4x_42C.txt`.


## bedgraph_genome_plots.Rmd - zoomed-out read counts on transcripts

This script produces genome browser-style figures of SSD1 CRAC profiles from bedgraph data. It focuses on 30°C rep A (`Ssd1_3_30`) and B (`Ssd1_4_30`), showing profiles on extended full-length transcripts.

This script relies on bedgraph files for both plus and minus strands, produced by `genomeCoverageBed` while running the pipeline, put in the directory `/Ssd1_CRAC_demult_dedup_20190114_all/bedgraph_genomecov/`. So it must be run after  `src/CRAC_pipeline_SE_demult_dedup.py`, see the repository `README.md` for details.


## pileup_Ssd1_closeup.Rmd - close-up view of read counts and deletions/mutations on transcripts

Pileup plots for Ssd1 CRAC data, only on specific regions of specific Ssd1 target genes that are specified in `input_annotation/Ssd1TargetGeneNamesOnly.txt`. These show detailed profiles of read counts, including nucleotide-specific mutations and deletions, along selected transcripts. Search for instances of the `CNYUCNYU` Ssd1-associated motif in those transcripts, and the `CCAACU` upstream motif.

This script relies on "pileup" files in tab-separated text format produced by `pyPileup` script while running the pipeline, put in the directory `/Ssd1_CRAC_demult_dedup_20190114_all/pyPileup_analyses/`. So it must be run after  `src/CRAC_pipeline_SE_demult_dedup.py`, see the repository `README.md` for details.


## peakanalysis_Ssd1.Rmd

This analysis searches for motifs associated with Ssd1-bound RNA from the CRAC datasets.
This script relies on Ssd1-bound peak data in .gtf (gff-like) format produced by `pyCalculateFDRs.py` script while running the pipeline, put in the directory `/Ssd1_CRAC_demult_dedup_20190114_all/pyCalculateFDRs_analyses/`. 

The objectives are to:
* filter Ssd1 CRAC hits by read count and width in addition to FDR
* generate fasta files of peak sequences
* discover motifs enriched in filtered peaks

It takes as input: 
* locations of Ssd1-bound peaks, with peak height and false discovery rate (FDR) as output by `pyCalculateFDRs.py` in gff-like format
* genome sequence in fasta format (EF4.74 annotation, R64-1-1 genome build)

It outputs, to directory `results`:
* sequences associated with filtered peak list, in fasta format
* MEME motif analysis of sequences of top 100 peaks.
