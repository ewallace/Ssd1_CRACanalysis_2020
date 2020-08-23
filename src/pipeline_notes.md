# Notes on running the CRAC pipeline



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
test_pipeline_flexbar_SE.py \
  -f /homes/ewallac2/mount/datastore/wallace_rna/bigdata/2019/Ssd1CRAC-2019-01-15/20190114_Ssd1_CRAC_init100000.fastq \
  -p 8 \
  --name test_pipeline_flexbar_SE_init100000 \
  -a input/3primeadapter.fasta 
```