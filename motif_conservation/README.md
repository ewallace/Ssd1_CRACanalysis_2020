
## motif_conservation

Motif conservation analysis in alternative fungi.
This takes the `CNYUCNYU` RNA motif found via CRAC analysis (DNA: `CNYTCNYT`), as calculated by code in `rmarkdown/peakanalysis_Ssd1.html`.
It searches for this `CNYTCNYT` in 5'UTR regulatory sequences in _Candida albicans_, _Schizosaccharomyces pombe_, and _Aspergillus fumigatus_, downloaded from [FungiDB](https://fungidb.org/)
For each of these organisms, the script reports gene ontology (GO) analysis of the list of genes that have at least 2 copies of `CNYTCNYT` in the 100nt of genomic sequence upstream of the annotated start codon.

We took data from fungidb for:

* C. albicans SC5314 (filenames start Ca)
* A.fumigatus Af293 (filenames start Af)
* S. pombe 972h- (filenames start Sp)

## src

Contains R markdown (.Rmd) files implementing the motif counting and analysis.
Also, the html files output by this analysis.
Full explanations are given in the files themselves.

* Afumigatus-Ssd1-motifs-fixedupstream.Rmd/.html
* Calbicans-Ssd1-motifs-fixedupstream.Rmd/.html
* Spombe-Ssd1-motifs-fixedupstream.Rmd/.html

## data

All genes in an organism [found by fungidb's "my strategies"](https://fungidb.org/fungidb/app/workspace/strategies/),
then selected genomic sequence from 1000nt to 1nt upstream of ATG start codon. Downloaded January/February 2021.

* Afumigatus_Af293_ATG_upstream_1000nt.fasta
* Calbicans_SC5314_ATG_upstream_1000nt.fasta
* Spombe_972h-_ATG_upstream_1000nt.fasta


## results

Lists of genes with at least 1 or at least 2 CNYTCNYT motifs in the 100nt upstream of ATG,
output by the .Rmd files in src.

* Af_id_list_CNYTCNYT1_up100.txt
* Af_id_list_CNYTCNYT2_up100.txt
* Ca_id_list_CNYTCNYT1_up100.txt
* Ca_id_list_CNYTCNYT2_up100.txt
* Sp_id_list_CNYTCNYT1_up100.txt
* Sp_id_list_CNYTCNYT2_up100.txt

We also do GO searches for the lists with at least 2 CNYTCNYT motifs in the 100nt upstream.
GO term finder results from Aspergillus genome database, Candida Genome database, and Princeton generic GO term finder, described in the .Rmd files in src.
Note that these were found using webforms so are not "self-contained" in this repository.

* Af_goFinderResult_11373_component_AspGD_id_list_CNYTCNYT2_up100.txt
* Ca_goFinderResult_28625_component_CGD_id_list_CNYTCNYT2_up100.txt
* Sp_goFinderResult_component_pombase_id_list_CNYTCNYT2_up100.txt
