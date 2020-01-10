
### Brief Description 

This repository contains files associated with our manuscript: "Phylogenetic comparative methods are problematic when applied to gene trees with speciation and duplication nodes: correcting for biases in testing the ortholog conjecture".

The files include:

Manuscript.Rmd - This file contains detailed text of our manuscript and the codes to reproduce our results. This analysis is based on reanalyses of a published paper by Dunn et al. (doi:10.1073/pnas.1707515115). We used the same compara gene tree file (ftp://ftp.ensembl.org/pub/release-75/emf/ensembl-compara/homologies/) as Dunn et al. used for their study. Before running our analysis, one need to run "manuscript_kernel.R" script of Dunn et al. using the required files available in the "kmrr" folder of their github repository. For our reanalysis, we ran the script and stored the output in the “manuscript_dunn.RData” file (https://doi.org/10.5281/zenodo.3604104). 

Premanuscript_run_TM.R - To save time during knitting of the "Manuscript.Rmd" file, one needs to run the script "Premanuscript_run_TM.R" first. The outputs are generated and saved as "Data_TMRR_latest.rda" available in https://doi.org/10.5281/zenodo.3604104.

TM.bib - The Bibliography file of our Rmarkdown file. 

functions_Dunn.R - The script adopted from Dunn et al. and renamed. The script contains functions, which were used to reproduce the results of Dunn et al. Some of the functions were further used in our study.      

funtion_TM_new.R - Additional script of functions written for our analysis.

plos.csl - CSL style used in our Manuscript.Rmd file.
