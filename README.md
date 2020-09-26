# Ssd1_CRACanalysis_2020

Analysis of CRAC data on Ssd1, August 2020


## The pipeline for read processing

This is "as self--contained as possible". Most files are kept locally, with the exception of the 
* .fastq data (very large)
* .novoindex aligner index (large)

To run the pipeline, use this code:

```
python src/CRAC_pipeline_SE_demult_dedup.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC.fastq  \
  -c input_annotation/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf input_annotation/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 5 \
  --name Ssd1_CRAC_demult_dedup_20190114_all \
  -a input_barcodes/3primeadapter.fasta \
  -b input_barcodes/Ssd1CRACSampleBarcodes.txt \
  --genelist input_annotation/Ssd1TargetGeneNamesOnly.txt \
  --genometab input_annotation/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab \
  --transcriptgff input_annotation/abundant_verified_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff
```

Unfortunately, due to a bug we have not yet been able to fix, this needs to be run *twice* for all tasks to complete.