Scripts and output for the global multi-mutant analysis (GMMA) by Johansson, Lindorff-Larsen and Winther.

Runs with Rscript version 4.0.2 or newer. 
Dependencies: igraph, minpack.lm

Get the original data file "amino_acid_genotypes_to_brightness_parsed.tsv" from
> Sarkisyan, K., Bolotin, D., Meer, M. et al. Local fitness landscape of the green fluorescent protein.
> Nature 533, 397–401 (2016). https://doi.org/10.1038/nature17995

```bash
Rscript gmma00_fetch_Sarkisyan.r
Rscript gmma01_structure.r amino_acid_genotypes_to_brightness_parsed.tsv
Rscript gmma_assign.r gmma_structured.rda
Rscript gmma02_fit_individual.r
Rscript gmma03_graph.r gmma_structured.rda
Rscript gmma04_global_estimation.r
Rscript gmma05_analysis.r gmma_fit_global.rda


The results are included in the output directory. Note that we use the residue numbering of the
original data which is shifted 2 positions compared to conventional numbering (uniprot, fpbase.org, etc).

Plot from the paper may be generated based on the content of the output directory using R-Studio
notebook, plots.Rmd

Enjoy!
