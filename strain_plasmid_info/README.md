# strain_plasmid_info

This folder contains information on yeast strain construction for Ssd1 analysis.

Maps are all in two files, one in snapgene format (`.dna`) and also in plain-text genbank format (`.gb`).
We just give the stem of the filename in the tables.

## guideRNAPlasmids - Plasmid maps for guide RNA plasmids for CRISPR

These are the plasmid maps for guide RNA plasmids based on pML104 from [Laughery & Wyrick](https://doi.org/10.1002/cpmb.110).
Each filename starts with the plasmid id, e.g. pRB0001.


| File | Description |
|-|-|
| pRB0001_CtermSsd1gRNA1 | gRNA plasmid used to add the  Cterminal HFtag to the  SSD1ORF  |
| pRB0002_NtermSsd1gRNA1 | gRNA plasmid used to add the Nterminal FHtag +/- the first 338 amino acids of Ssd1 |
| pRB0003_Cl1_RNB_gRNA1 | gRNA plasmid used to make point mutations in the RNB cluster |
| pRB0004_Cl2_CSD-side-gRNA1 | gRNA plasmid used to make the K622E and K623E point mutations in the CSD-side cluster |
| pRB0005_Cl3_CSD-top_gRNA1 | gRNA plasmid used to make the point mutations in the CSD-top cluster |
| pRB0006_Cl4_CSD1-insert_gRNA1 | gRNA plasmid used to make the point mutations in the CSD1-insert cluster |
| pRB0007_Cl2_CSD-sideR549EgRNA1 | gRNA plasmid used to add the R549E point mutation to the CSD-side cluster |


## Ssd1LocusMaps - Ssd1 locus maps detailing sequence and edits

These are maps of the Ssd1 locus and edited locus in _S. cerevisiae_ yeast.
They are based on the SSD1/YDR293C ORF with 1000nt of genomic sequence on each side, wild-type sequence for S288C yeast from the file `orf_genomic_1000.fasta` as downloaded from the [saccharomyces genome database](https://www.yeastgenome.org/), genome release R64-2-1.

The scarless tag templates (in snapgene format) show the guide RNAs as well as the repair templates used to edit the locus.

| File | Description
|-|-|
| SSD1_HTP_URA3selPlus | Tandem tagged HTP with K. lactis SSD1 UTR and K. lactis URA3 selection cassette. |
| Ssd1-CtermHFconstruct | Ssd1 scarless tag with HF tag at C-terminus |
| NTermFH-Ssd1_construct | Ssd1 scarless tag with FH at N-terminus |
| NTermFH-Ssd1_delN338_construct | Ssd1 with amino acids 1-338 deleted and tagged with FH at N-terminus |
| Ssd1Cl1_RNB_Mutations | Ssd1 triple point mutants on RNB domain: R745E, R746E, R748E |
| PartialSsd1Cl2_CSD-side_Mutations | Ssd1 double point mutants on CSD-side: K622E, R623E  |
| Ssd1Cl2_CSD-side_Mutations_plusR549E | Ssd1 triple point mutants on CSD-side: R549E, K622E, R623E |
| Ssd1Cl3_CSD-top_Mutations | Ssd1 double point mutants on CSD-top: W583A, K585A |
| Ssd1Cl4_CSD1-insert_Mutations | Ssd1 point mutants on CSD1 insert: K505E, Q506E, Q510E |
