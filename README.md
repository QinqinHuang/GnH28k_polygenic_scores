# Transferability of genetic loci and polygenic scores for cardiometabolic traits in British Pakistanis and Bangladeshis

Code used to calculate polygenic scores (PGS) and assess its performance. 

Preprint: https://www.medrxiv.org/content/10.1101/2021.06.22.21259323v1. Under review.


## 1. scores from the PGS Catalog
1 and 2: calculating previously developed scores downloaded from the PGS Catalog

3: assessing the predictive accuracy of the scores

4: quintile analysis


## 2. PGSs constructed using the Clumping and P-value Thresholding method
1: preparing reference EUR data using 1000 Genomes

2: a function used to match GWAS summary statistics with the validation cohort

3: calculating PGS using the Clumping and P-value Thresholding (C+T) method

4: assessing the predictive accuracy of the C+T scores and reporting the best one

## 3. Marquez-Luna metaPGS method
This method (https://onlinelibrary.wiley.com/doi/10.1002/gepi.22083) combines a PGS derived from large EUR GWAS and a PGS derived from smaller GWAS in matched ancestry. We downloaded SAS GWAS from PanUKBB (https://pan.ukbb.broadinstitute.org) to construct a C+T SAS PGS.

1: preparing reference SAS data using 1000 Genomes

2: calculating SAS PGS using the Clumping and P-value Thresholding (C+T) method

3: selecting the best C+T SAS PGS

4: calcualting the metaPGS and evaluating the performance

## 4. Calculating QRISK3 scores and clinical utility of PGS
