# Notes on running the CRAC pipeline

This file contains somewhat disorganised notes on multiple pilot attempts to run the CRAC pipeline. 

Development was carried out at https://git.ecdf.ed.ac.uk/ewallac2/crac_pipelines.

```
CRAC_pipeline_SE.py  \
  -f name_of_fastq_file  \
  -c path_to_chromosome_length_file  \
  --gtf name_of_ gtf_file  \
  --novoindex path_to_novoindexfile  \
  -p numberofprocessorstouse \
  --name name_for_run  \
  -a path_to_adapter_file  \
  -b path_to _list_of_barcodes
```

Example from Rosey Bayne, 2019.

```
CRAC_pipeline_SE.py \
  -f 20190114_Ssd1_CRAC.fastq \
  -c Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
  --novoindex Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_Jan19 \
  -a 3primeadapter.fasta \
  -b Ssd1_Barcodes.txt
```

Test run on initial 100000 records:

```
CRAC_pipeline_SE.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC_init100000.fastq \
  -c Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
  --novoindex Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_pilot_init100000_test22Aug2020 \
  -a input/3primeadapter.fasta \
  -b input/Ssd1_Barcodes.txt
```

TO DO: Find these files.

- `Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths`
- `Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf`
- `Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex`

## Small test

Test run on initial 100000 records:

```
src/test_pipeline_flexbar_SE.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC_init100000.fastq \
  -p 8 \
  --name test_pipeline_flexbar_SE_init100000 \
  -a input/3primeadapter.fasta 
```


```
python test_pipeline_flexbar_SE.py \
  -f data/20181101_Ssd1_CRAC_init10000.fastq \
  -p 4 \
  -a data/3primeadapter.fasta 
```

## 10000 records on the CRAC SE pipeline

```
python CRAC_pipeline_SE.py \
  -f data/20181101_Ssd1_CRAC_init10000.fastq \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_pilot_init10000_collapse_24Aug2020 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt
```

```
python CRAC_pipeline_SE.py \
  -f data/20181101_Ssd1_CRAC_init10000.fastq \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --nocollapse \
  --name Ssd1_CRAC_pilot_init10000_NOcollapse_24Aug2020 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt
```


## 10000 records on the CRAC SE demult pipeline

```
python CRAC_pipeline_SE_demultonly.py \
  -f data/20181101_Ssd1_CRAC_init10000.fastq \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_pilot_demultonly_init10000_4Sep2020 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt
```

## 10000 records on the CRAC SE demult pipeline, adding pypileup

```
python CRAC_pipeline_SE_demultonly.py \
  -f data/20181101_Ssd1_CRAC_init10000.fastq \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_pilot_demultonly_init10000_8Sep2020 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt \
  --genelist data/Ssd1TargetGeneNamesOnly.txt \
  --genometab /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab
```

## All Techrep 1 (poor quality) records on the CRAC SE demult pipeline

```
python CRAC_pipeline_SE_demultonly.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2018/Ssd1CRAC-2018-11-02/20181101_Ssd1_CRAC.fastq \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_demultonly_20181101_all \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt \
  --genelist data/Ssd1TargetGeneNamesOnly.txt \
  --genometab /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab
```

## Init 100000 Techrep 2 (good quality) records on the CRAC SE demult pipeline

```
python CRAC_pipeline_SE_demultonly.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC_init100000.fastq  \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_demultonly_20190114_init100000 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt \
  --genelist data/Ssd1TargetGeneNamesOnly.txt \
  --genometab /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab
```

```
python CRAC_pipeline_SE_demultonly.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC_init100000.fastq  \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_demultonly_20190114_init100000_multicov \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt \
  --genelist data/Ssd1TargetGeneNamesOnly.txt \
  --genometab /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab \
  --transcriptgff /homes/ewallac2/yeastutrgff/out/abundant_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff
```

## All Techrep 2 (good quality) records on the CRAC SE demult pipeline

```
python CRAC_pipeline_SE_demultonly.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC.fastq  \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_demultonly_20190114_all_run8Sep2020 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt \
  --genelist data/Ssd1TargetGeneNamesOnly.txt \
  --genometab /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab
```

## Init 100000 Techrep 2 (good quality) records on the CRAC SE demult dedup pipeline

```
python CRAC_pipeline_SE_demult_dedup.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC_init100000.fastq \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_demult_dedup_20190114_init100000 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt \
  --genelist data/Ssd1TargetGeneNamesOnly.txt \
  --genometab /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab \
  --transcriptgff /homes/ewallac2/yeastutrgff/out/abundant_verified_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff
```

## All Techrep 2 (good quality) records on the CRAC SE demult dedup pipeline

```
python CRAC_pipeline_SE_demult_dedup.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC.fastq  \
  -c /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.lengths \
  --gtf data/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.quotefix.gtf \
  --novoindex /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.novoindex \
  -p 8 \
  --name Ssd1_CRAC_demult_dedup_20190114_all_run25Sep2020 \
  -a data/3primeadapter.fasta \
  -b data/Ssd1CRACSampleBarcodes.txt \
  --genelist data/Ssd1TargetGeneNamesOnly.txt \
  --genometab /homes/rbayne2/CRACAnalysisJan2019_1/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.shortChrNames.fa.tab \
  --transcriptgff /homes/ewallac2/yeastutrgff/out/abundant_verified_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff
```

## Create STAR index

```
STAR --runThreadN 12 --runMode genomeGenerate \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--genomeFastaFiles /homes/genomes/s.cerevisiae/sacCer3/sacCer3.fa \
--genomeSAindexNbases 10 \
--sjdbGTFfile /homes/genomes/s.cerevisiae/sacCer3/annotation/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
--sjdbOverhang 49 \
2> /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a/Log.out
```

```
STAR-2.7.5a --runThreadN 12 --runMode genomeGenerate \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.5a \
--genomeFastaFiles /homes/genomes/s.cerevisiae/sacCer3/sacCer3.fa \
--genomeSAindexNbases 10 \
--sjdbGTFfile /homes/genomes/s.cerevisiae/sacCer3/annotation/Saccharomyces_cerevisiae.EF4.74_SGDv64_CUTandSUT_withUTRs_noEstimates_antisense_intergenic_4xlncRNAs_final.pyCheckGTFfile.output.gtf \
--sjdbOverhang 49 \
2> /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.5a/Log.out
```

## Test alignment to STAR 2.7.3a

```
STAR --runThreadN 12 \
--alignIntronMax 1500 \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--readFilesIn Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/demultiplexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.fastq \
--outFileNamePrefix Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_default
```

```
STAR --runThreadN 12 \
--alignIntronMax 1500 \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--readFilesIn Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/demultiplexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.fastq \
--chimSegmentMin 15 --chimOutType WithinBAM --outSAMtype BAM Unsorted \
--outFileNamePrefix Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_chim15_
```

```
STAR --runThreadN 12 \
--alignIntronMax 1500 \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--readFilesIn Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/demultiplexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.fastq \
--outSAMmultNmax 1 \
--outFileNamePrefix Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_multNmax1_
```


```
STAR --runThreadN 12 \
--alignIntronMax 1500 \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--readFilesIn Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/demultiplexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.fastq \
--scoreDelOpen 0 \
--outFileNamePrefix Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_
```

```
STAR --runThreadN 12 \
--alignIntronMax 1500 \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--readFilesIn Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/demultiplexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.fastq \
--scoreDelOpen 0 \
--outSAMmultNmax 1 \
--outFileNamePrefix Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_
```


```
STAR --runThreadN 12 \
--alignIntronMax 1500 \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--readFilesIn Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/demultiplexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.fastq \
--scoreDelOpen 0 \
--outSAMmultNmax 1 \
--scoreGap -10 --outSJfilterCountUniqueMin -1 -1 -1 -1 \
--seedSearchStartLmax 30 \
--outFileNamePrefix Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap10_nonewSJ_seedL30_
```

```
STAR --runThreadN 12 \
--alignIntronMax 1500 \
--genomeDir /homes/ewallac2/EW_genomes/Scer_R64_STAR_2.7.3a \
--readFilesIn Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/demultiplexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.fastq \
--scoreDelOpen 0 \
--outSAMmultNmax 1 \
--scoreGap -20 --outSJfilterCountUniqueMin -1 -1 -1 -1 \
--seedSearchStartLmax 20 \
--outFileNamePrefix Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap20_nonewSJ_seedL20_
```

```
samtools view -b Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_defaultAligned.out.sam | \
samtools sort -@ 3 -O bam -o Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_defaultAligned.out.sorted.bam
samtools index Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_defaultAligned.out.sorted.bam

samtools view -b Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_Aligned.out.sam | \
samtools sort -@ 3 -O bam -o Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_Aligned.out.sorted.bam
samtools index Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_Aligned.out.sorted.bam

samtools view -b Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap10_nonewSJ_seedL30_Aligned.out.sam | \
samtools sort -@ 3 -O bam -o Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap10_nonewSJ_seedL30_Aligned.out.sorted.bam
samtools index Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap10_nonewSJ_seedL30_Aligned.out.sorted.bam

samtools view -b Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap20_nonewSJ_seedL20_Aligned.out.sam | \
samtools sort -@ 3 -O bam -o Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap20_nonewSJ_seedL20_Aligned.out.sorted.bam
samtools index Ssd1_CRAC_demultonly_20190114_all_run8Sep2020_qtrim/aligner_tests/SSD1_3_30_STAR_nodelpen_multNmax1_scoreGap20_nonewSJ_seedL20_Aligned.out.sorted.bam


```

## multiBamCoverage

```
mkdir Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/multicov_analyses/
echo "Ssd1_3_30\tSsd1_4_40\tSsd1_3_42\tSsd1_4_42\tBY4741" > Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/multicov_analyses/samplenames.txt
multiBamCov -s -bed /homes/ewallac2/yeastutrgff/out/abundant_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff \
 -bams \
Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/aligned_bamindexed/20190114_Ssd1_CRAC_trimmed_NNNGTGAGC_SSD1_3_30.bam \
Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/aligned_bamindexed/20190114_Ssd1_CRAC_trimmed_NNNTGGAGC_SSD1_4_30.bam \
Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/aligned_bamindexed/20190114_Ssd1_CRAC_trimmed_NNNAGAGC_SSD1_3_42.bam \
Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/aligned_bamindexed/20190114_Ssd1_CRAC_trimmed_NNNCTAGC_SSD1_4_42.bam \
Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/aligned_bamindexed/20190114_Ssd1_CRAC_trimmed_NNNGACTTAGC_BY4741.bam \
 > Ssd1_CRAC_demultonly_20190114_all_dedup_run13Sep2020/multicov_analyses/allsample_transcriptcounts_nopipeline.txt
```