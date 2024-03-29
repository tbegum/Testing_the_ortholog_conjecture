
### Brief Description 

This repository contains all the required scripts, and files associated with ourpublished manuscript: 

Begum T, Robinson-Rechavi M (2021). Special care is needed in applying phylogenetic comparative methods to gene trees with speciation and duplication nodes. Molecular Biology and Evolution. https://doi.org/10.1093/molbev/msaa288.

This git repository includes:

README.md - It contains a brief description about all the required scripts, and files.

Model_fitting.R - This script is written to fix time calibration bias of old duplicates in a phylogeny with atleast one speciation and one duplication nodes. It also performs phylogenetic data modeling in a maximum-likelihood, and in a Bayesian frameworks. The outputs of this script are saved as "Model_fitting_TMRR.Rdata", and is archived at https://doi.org/10.5281/zenodo.4003391. 

Premanuscript_run_TMRR.R - This script is required to reanalyze the results of Dunn et al. (2018) using a Phylogenetic Independent Contrasts (PICs) method. It also provides a guideline to implement PIC on gene trees to study the effect of gene duplication, and to assess the results of PIC method. The outputs are saved as "Analyses_TMRR.Rdata", and is available at https://doi.org/10.5281/zenodo.4003391.The reproduced output of Dunn et al. (2018) is stored at the zenodo link mentioned above as "manuscript_dunn.RData". 

functions_Dunn.R - To reproduce the results of Dunn et al. (2018), we used their functions ("functions.R") from https://github.com/caseywdunn/comparative_expression_2017. The same is provided here by renaming it as "functions_Dunn.R". 

funtion_TM_new.R - All the functions required to reproduce the results of this study.

