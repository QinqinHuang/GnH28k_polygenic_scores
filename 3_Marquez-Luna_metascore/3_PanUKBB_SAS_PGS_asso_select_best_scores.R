#-------------------------------------------------------
# 2021-05-01
# 
# Aims:
#  (1) to calculate performance for each C+T PGS
#  (2) to select the best PGS
#  (3) to calculate boostrapping 95% CIs for the best PGS
#-------------------------------------------------------
### Functions are defined in: https://github.com/QinqinHuang/GnH28k_polygenic_scores/blob/main/My_functions_GnH_PGS.R
source("My_functions_GnH_PGS.R")


#----- set up -----
### update with your path file

# "traitGWAS": a file with all PGS-trait pairs that need to be tested

# "phenotable_2ndexcl": this table contains phenotype and covariate data for unrelated individuals. 
# code for removing related samples can be found in "remove_rel" and "remove_rel_pheno" functions in the "My_functions_GnH_PGS.R" file.


#----- C+T PGS accuracy -----
# Using the phenotype table for unrelated individuals (2nd and closer relatives excluded)
for(ii in 1:nrow(traitGWAS)) {
  
  write(paste0(" ** Phenotype: ", traitGWAS$phenotype[ii]), "../log.txt", append = T)
  
  validate_PGSs(phenotable = phenotable_2ndexcl, cohort = "ELGH", 
                trait = traitGWAS$phenotype[ii], 
                GWAS = traitGWAS$GWAS[ii],
                binary = ifelse(traitGWAS$traittype[ii]=="binary", yes=T, no=F),
                quandata = ifelse(traitGWAS$traittype[ii]!="binary", yes="IRN", no = "raw"),
                whichage = ifelse(traitGWAS$traittype[ii]!="binary", yes="at_event", no = "at_recrt"),
                PGSfolder = "~/P_drive/Polygenic_scores_Feb28k/Marquez-Luna_method/PRS_SAS_PanUKBB/all28k")
}


# summary stats of the full model
filelist = grep("Full", dir("Individual_files"), value = T)

# Load all data
allstats = foreach(ff = filelist, .combine = rbind) %do% {
  d = readRDS(paste0("Individual_files/", ff))
  return(d)
}

# reorder columns
allstats = allstats[, .(phenotype, trait, GWASsumm, ancestry,
                        N, prevalence, agedata, transformation,
                        LDr2, pval_threshold, bestscore,
                        full_R2, null_R2, incR2, incR2_CI_L, incR2_CI_U, incR2_bootstrapest,
                        AUC_full, AUC_null, incAUC, incAUC_CI_L, incAUC_CI_U, incAUC_bootstrapest,
                        PRS_beta, PRS_se, PRS_p )]


# best PGS per trait-GWAS pair
bestPGS = allstats[bestscore == "BestScore" & !is.na(incR2_CI_L)]













