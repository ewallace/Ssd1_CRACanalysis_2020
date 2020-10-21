# input_annotation

This directory contains genome sequences, annotations, and gene lists used to process CRAC and RNA-seq data.

## abundant_verified_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff

Messenger RNA transcript locations from end-to-end in .gff format. This is a modified version of the transcript maps from Saccharomyces Genome Database (Ng et al. 2020), using the "abundant transcript" data derived from (Pelechano, Wei, and Steinmetz 2013). We added default-length 25nt 5′UTRs, and 125nt 3′UTRs, as recorded in [yeastutrgff repository](https://github.com/ewallace/yeastutrgff).

## gff_ncRNAs_abundantverifiedmRNAparts.gff

Transcript annotations with non-coding RNAs, and mRNAs broken down into parts (5′UTR,CDS, intron, 3′UTR). For details, see [yeastutrgff repository](https://github.com/ewallace/yeastutrgff).

## Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf

Another gff of transcript features, from David Tollervey's lab. This includes mRNA annotations from EF4.74, and many classes of ncRNA including unstable transcripts (CUTs, SUTs).

## Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa

*Saccharomyces cerevisiae* genome sequence in .fasta format from ensemble release EF4.74, which is identical to SGD release R64-1-1. We edited the chromosome names to be `chrI, chrII, ...`.

Some of the scripts run in this repository create an index file for this fasta,  `Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.fai`. We do not include this in the repository as it can be made fresh locally.

## Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab

This is the data from `Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa`, but in a tab-separated variable format required by some of the tools in pyCRAC suite.

## Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths

Table of lengths of the chromosomes in `Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa`.

## scer-mrna-protein-absolute-estimate.txt

File with gene names, ORF names, and estimates of absolute molecule counts from:

Csárdi, Gábor et al. (2016), Data from: Accounting for experimental noise reveals that mRNA levels, amplified by post-transcriptional processes, largely determine steady-state protein levels in yeast, Dryad, Dataset, https://doi.org/10.5061/dryad.d644f


## Ssd1TargetGeneNamesOnly.txt

The names of select Ssd1 target genes, that are used as input for detailed pileup maps in `pyPileup.py` as run from the pipeline.