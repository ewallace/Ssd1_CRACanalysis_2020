# GOterms

Gene Ontology terms for lists of genes producing Ssd1-bound mRNAs, from resources at Saccharomyces Genome Database https://www.yeastgenome.org/

Run by Edward Wallace on 6th February 2021, from gene lists output by the script `deseq2_Ssd1_CRAC_vs_RNAseq.Rmd`.
* `Ssd1_enrichment_DeSeq2_30C_2x_padj0p05_genesonly.txt`
* `Ssd1_enrichment_DeSeq2_42C_2x_padj0p05_genesonly.txt`

From`Ssd1_targets_DeSeq2_42C_4x_padj0p05_genesonly.txt`, removed TAR1/YLR154W-C removed because the name is ambiguous. This is anyway antisense to rRNA.

## GO terms enriched/overrepresented in Ssd1 targets

The GO Term Finder (Version 0.86) searches for significant shared GO terms, or parents of those GO terms, used to describe the genes in your list to help you discover what the genes may have in common. 
https://www.yeastgenome.org/goTermFinder

* `Ssd1_targets_DeSeq2_30C_component_GOenriched.txt`
* `Ssd1_targets_DeSeq2_30C_function_GOenriched.txt`
* `Ssd1_targets_DeSeq2_30C_process_GOenriched.txt`
* `Ssd1_targets_DeSeq2_42C_component_GOenriched.txt`
* `Ssd1_targets_DeSeq2_42C_function_GOenriched.txt`
* `Ssd1_targets_DeSeq2_42C_process_GOenriched.txt`

## GO slim terms

The GO Slim Mapper maps annotations of a group of genes to more general terms and/or bins them into broad categories, i.e. GO Slim terms. It is not enrichment.

* `Ssd1_targets_DeSeq2_30C_component_slimTab.txt`
* `Ssd1_targets_DeSeq2_30C_function_slimTab.txt`
* `Ssd1_targets_DeSeq2_30C_process_slimTab.txt`
