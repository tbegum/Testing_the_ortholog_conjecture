
### Brief Description 

This repository contains required scripts, and files associated with our manuscript: "Special care is needed in applying phylogenetic comparative methods to gene trees with speciation and duplication nodes".

This includes:

README.md - This contains a brief description about all the scripts, and files provided here.

Model_fitting.R - This script is written to fix time calibration bias of old duplicates in a phylogeny with speciation and duplication nodes. It also performs phylogenetic data modeling in a maximum-likelihood, and in a Bayesian frameworks. The outputs of this script are saved as "Model_fitting_TMRR.Rdata", and is archived at https://doi.org/10.5281/zenodo.4003391. 

Premanuscript_run_TMRR.R - This script is required to reanalyze the results of Dunn et al. (2018) using a Phylogenetic Independent Contrasts (PIC) method. It also provides a guideline to implement PIC on gene trees to study the effect of gene duplication, and to assess the results of PIC method. The outputs are saved as "Analyses_TMRR.Rdata", and is available at https://doi.org/10.5281/zenodo.4003391. 

functions_Dunn.R - To reproduce the results of Dunn et al. (2018), we used their functions ("functions.R") from https://github.com/caseywdunn/comparative_expression_2017. The same is provided here by renaming it as "functions_Dunn.R". 

funtion_TM_new.R - All the functions required to reproduce the results of this study.

Manuscript.Rmd - This R markdown file contains our  manuscript text, and source codes. To knit the markdown file, one needs to have outputs available at zenodo (https://doi.org/10.5281/zenodo.4003391). This includes the output of the study of Dunn et al. (2018), which was reproduced by using the script "manuscript_kernel.R" of them, and by using the necessary files provided on their github (https://github.com/caseywdunn/comparative_expression_2017). The reproduced output of Dunn et al. (2018) is stored at the zenodo link mentioned above as "manuscript_dunn.RData".  

TM.bib - The required bibliography to knit the Rmarkdown file. 

molecular-biology-and-evolution.csl - CSL style used in our Rmarkdown.

