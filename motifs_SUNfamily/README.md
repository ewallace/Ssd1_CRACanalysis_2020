# motifs_SUNfamily

Edward Wallace, July 2020, revised June 2021

# Summary

Analysis of S. cerevisiae SIM1/UTH1/NCA3/SUN4 homologs in some ascomycetes, to look for conservation of Ssd1 regulation. 
I searched for the Ssd1 recognition element CNYTCNYT in their transcript sequences, 5'UTRs, Coding sequence, and 3'UTRs
We also looked for upstream-element-like motifs (CCAACT) immediately preceeding those.

We looked in species:
- Saccharomyces cerevisiae (Sc)
- Candida albicans (Ca)
- Aspergillus fumigatus (Af)
- Neurospora crassa (Nc)
- Schizosaccharomyces pombe (Sp)

These are the genes in PANTHERDB family [PTHR31316:SF0](http://www.pantherdb.org/panther/family.do?clsAccession=PTHR31316:SF0).

# Code used for multiple sequence alignments and tree

These lines in bash were used to generate the multiple sequence alignments (with MAFFT) and approximate phylogenetic trees (with fasttree).

```bash
mafft --maxiterate 1000 --genafpair SUN4_homologs_proteinsequences.fasta > SUN4_homologs_proteinsequences_mafft_geinsi.fasta

fasttree SUN4_homologs_proteinsequences_mafft_geinsi.fasta > SUN4_homologs_proteinsequences_mafft_geinsi_fasttree.nwk
```

# Contents of folder

- SUN4_homologs_proteinsequences.fasta - just the protein sequences for the relevant homologs
- SUN4_homologs_features.txt - features associated with Ssd1 regulation, compiled from locus_maps below
- SUN4_homologs_proteinsequences_mafft_geinsi.fasta - multiple sequence alignment made by MAFFT (see code above)
- SUN4_homologs_proteinsequences_mafft_geinsi_fasttree.nwk - phylogenetic tree made by fasttree (see code above)
- SUN4_homologs_tree_Ssd1targets.Rmd - makes figures from those data
- SUN4_homologs_tree_Ssd1targets.html - html output with results and figures

## locus_maps -  Snapgene .dna and Genbank .gb files with motifs illustrated

Files like `xxxx_genomic1000.dna` is a snapgene-format file of genomic sequence for gene xxxx including introns and other features 1000nt upstream of start codon and 1000nt downstream of stop.
Files like `xxxx_genomic1000.gb` are genbank-format (plain-text) files exported from the original snapgene file.

For Candida albicans SUN41 we had to go further upstream because its 5'UTR including intron is longer than 1000nt.
These files are based on sequences retrieved from FungiDB.

- AfSUN1_SUN4homolog_AFUA_7G05450_genomic1500up_1000down.dna
- CaSIM1_UTH1LDO_C1_13940W_A_C_albicans_SC5314.dna
- CaSUN41_SUN4LDO_C6_00820W_A_flanking_intergenic.dna
- Nc-ghx-3_SUN4homolog_NCU02668_genomic_1500up_1000down.dna
- ScNCA3_YJL116C_genomic1000.dna
- ScSIM1_YIL123W_genomic1000.dna
- ScSUN4_YNL066W_genomic1000.dna
- ScUTH1_YKR042W_genomic1000.dna
- SpPSU1_SPAC1002.13c_genomic1000.dna
- SpPSU2_SPBC2G2.17c_genomic1000.dna
