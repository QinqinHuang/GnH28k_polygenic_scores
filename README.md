# Transferability of genetic loci and polygenic scores for cardiometabolic traits in British Pakistanis and Bangladeshis

Code used to calculate polygenic scores (PGS) and assess its performance. 

Preprint: https://www.medrxiv.org/content/10.1101/2021.06.22.21259323v1. 
Published in Nature Communications: https://www.nature.com/articles/s41467-022-32095-5


## 1. scores from the PGS Catalog
1 and 2: calculating previously developed scores downloaded from the PGS Catalog

3: assessing the predictive accuracy of the scores

4: quintile analysis


## 2. PGSs constructed using the Clumping and P-value Thresholding method
1: preparing LD reference using the 1000 Genomes EUR superpopulation

2: a function used to match GWAS summary statistics with SNPs in the validation cohort

3: calculating PGS using the Clumping and P-value Thresholding (C+T) method

4: assessing the predictive accuracy of the C+T scores and reporting the best one

## 3. Marquez-Luna metaPGS method
This method (https://onlinelibrary.wiley.com/doi/10.1002/gepi.22083) combines a PGS derived from large EUR GWAS and a PGS derived from smaller GWAS in matched ancestry. We downloaded SAS GWAS from PanUKBB (https://pan.ukbb.broadinstitute.org) to construct a C+T SAS PGS.

1: preparing LD reference using the 1000 Genomes SAS superpopulation

2: calculating SAS PGS using the C+T method

3: selecting the best C+T SAS PGS

4: calcualting the metaPGS and evaluating its performance

## 4. Calculating QRISK3 scores and clinical utility of PGS

1: preparing QRISK3 variables

2: performing multiple imputation to impute missing data in QRISK3 continuous variables, and calculating QRISK3 scores

3: calculating the integrated score combining QRISK3 and PGS, running the Cox regression model to estimate C-index, and calculating the net reclasssification indices (NRI; both categorical and continuous NRI)


#### Code for calculating power for replicating GWAS loci is available at https://github.com/Nsallah1/GH_Manuscript
