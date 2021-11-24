#-------------------------------------------------------
# 2021-06-01
# 
# Aims:
#  to calculate the metaPGS (Marquez-Luna et al)
#
# Model:
#  pheno ~ EURPRS + SASPRS + cov
#-------------------------------------------------------
### Functions are defined in: https://github.com/QinqinHuang/GnH28k_polygenic_scores/blob/main/My_functions_GnH_PGS.R
source("My_functions_GnH_PGS.R")


MetaModel = formula("pheno ~ PRS_EUR + PRS_SAS +gender+age+I(age^2)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")
MetaModelNULL = formula("pheno ~ gender+age+I(age^2)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")


#----- set up -----
### update with your path file

# "phenotable_2ndexcl": this table contains phenotype and covariate data for unrelated individuals. 
# code for removing related samples can be found in "remove_rel" and "remove_rel_pheno" functions in the "My_functions_GnH_PGS.R" file.


## Load parameters for the best PGS derived from European GWASs
EURPRSbest = readRDS("Best_CTPGS_ELGH.RDS")

## Load parameters for the best PGS derived from PanUKBB SAS GWASs
SASPRSbest = readRDS("Best_PGS_binary_quantitative_pairs_ELGH.RDS")

phenolist = merge(EURPRSbest[,.(trait, phenotype, GWAS_EUR = GWASsumm, 
                                LDr2_EUR = LDr2, pval_threshold_EUR = pval_threshold)], 
                  SASPRSbest[,.(trait, phenotype, GWAS_SAS = GWASsumm,
                                LDr2_SAS = LDr2, pval_threshold_SAS = pval_threshold,
                                isbinary = ifelse(is.na(AUC_full), yes = F, no = T))],
                  by = c("trait","phenotype"))


#----- metaPGS -----
metasumm = foreach(ii = 1:nrow(phenolist), .combine = rbind) %do% {
  
  #-- prepare phenotype (unrelated) and cov table --
  mydd = phenotable_2ndexcl
  mydd$pheno = mydd[, get(phenolist$phenotype[ii])]
  mydd = mydd[!is.na(pheno)]
  
  # IRN for binary traits
  if(!phenolist$isbinary[ii]) {
    mydd$pheno = rankNorm(mydd$pheno)
  }
  
  # Age: age at recruitment or age at measurement
  if(phenolist$isbinary[ii]) {
    mydd[, age := age_at_recrt]
  } else {
    cc = paste0("age_", phenolist$phenotype[ii])
    mydd[, age := get(cc)]
  }
  
  #-- best PGSs for this phenotype -- 
  myEURPRS = load_bestPGS(whichGWASsumm = phenolist$GWAS_EUR[ii], ifonUKSerp = T, 
                          LDr2 = phenolist$LDr2_EUR[ii], pval = phenolist$pval_threshold_EUR[ii],
                          scalePGS = T)
  mySASPRS = load_bestPGS(PGSfilename = paste0("~/P_drive/Polygenic_scores_Feb28k/Marquez-Luna_method/PRS_SAS_PanUKBB/all28k/",phenolist$GWAS_SAS[ii],"_GAsP_bgen_dosage/PRSice_output_",phenolist$LDr2_SAS[ii],".all.score"), ifonUKSerp = T, 
                          LDr2 = phenolist$LDr2_SAS[ii], pval = phenolist$pval_threshold_SAS[ii],
                          scalePGS = T)
  
  myPRS = merge(myEURPRS[,.(ID, PRS_EUR = PRS_scaled)],
                mySASPRS[,.(ID, PRS_SAS = PRS_scaled)], by = "ID")
  
  # merge pheno and PRS 
  combined = merge(mydd, myPRS, by = "ID")
  
  
  # Binary traits
  if(phenolist$isbinary[ii]) {
    # Logistic regression
    myfit_full = glm(MetaModel, data = combined, family = binomial())
    myfit_null = glm(MetaModelNULL, data = combined, family = binomial())
    
    returndd = data.table(N = nrow(combined),
                          prevalence = mean(combined$pheno),
                          Ncases = sum(combined$pheno),
                          AUC_full = roc(combined$pheno, predict(myfit_full, type = c("response")))$auc[1],
                          AUC_null = roc(combined$pheno, predict(myfit_null, type = c("response")))$auc[1],
                          PRS_EUR_beta = summary(myfit_full)$coefficients["PRS_EUR",1],
                          PRS_EUR_se = summary(myfit_full)$coefficients["PRS_EUR",2],
                          PRS_EUR_p = summary(myfit_full)$coefficients["PRS_EUR",4],
                          PRS_SAS_beta = summary(myfit_full)$coefficients["PRS_SAS",1],
                          PRS_SAS_se = summary(myfit_full)$coefficients["PRS_SAS",2],
                          PRS_SAS_p = summary(myfit_full)$coefficients["PRS_SAS",4],
                          full_R2 = NA,
                          null_R2 = NA )
    
    returndd[, incR2 := NA]
    returndd[, incAUC := AUC_full - AUC_null]
  
    # run 1000 bootstraps to estimate 95% CI: incAUC
    incAUC_resamples = vector()
    while(length(incAUC_resamples) <1000) {
      resamples = combined[sample(1:nrow(combined), size = nrow(combined), replace = T)]
      myfit_full = glm(MetaModel, data = resamples, family = binomial())
      myfit_null = glm(MetaModelNULL, data = resamples, family = binomial())
      myAUCfull = roc(resamples$pheno, predict(myfit_full, type = c("response")))$auc[1]
      myAUCnull = roc(resamples$pheno, predict(myfit_null, type = c("response")))$auc[1]
      newest = myAUCfull - myAUCnull
      
      incAUC_resamples = append(incAUC_resamples, newest)
    }
    
    incAUC_resamples = sort(incAUC_resamples[1:1000])
    returndd$incAUC_CI_L = incAUC_resamples[25]
    returndd$incAUC_CI_U = incAUC_resamples[975]
    returndd$incAUC_bootstrapest = list(incAUC_resamples)
    
    returndd$incR2_CI_L = NA
    returndd$incR2_CI_U = NA
    returndd$incR2_bootstrapest = NA
    
  }  else {  # Quantitative traits
    
    # Linear regression
    myfit_full = summary(lm(MetaModel, data = combined))
    myfit_null = summary(lm(MetaModelNULL, data = combined))
    
    returndd = data.table(N = nrow(combined),
                          prevalence = NA,
                          Ncases = NA,
                          AUC_full = NA, 
                          AUC_null = NA,
                          PRS_EUR_beta = myfit_full$coefficients["PRS_EUR",1],
                          PRS_EUR_se = myfit_full$coefficients["PRS_EUR",2],
                          PRS_EUR_p = myfit_full$coefficients["PRS_EUR",4],
                          PRS_SAS_beta = myfit_full$coefficients["PRS_SAS",1],
                          PRS_SAS_se = myfit_full$coefficients["PRS_SAS",2],
                          PRS_SAS_p = myfit_full$coefficients["PRS_SAS",4],
                          full_R2 = myfit_full$r.squared,
                          null_R2 = myfit_null$r.squared)
    
    returndd[, incR2 := full_R2 - null_R2]
    returndd[, incAUC := NA]
    
    
    # run 1000 bootstraps to estimate 95% CI
    incR2_resamples = vector()
    while(length(incR2_resamples) <1000) {
      resamples = combined[sample(1:nrow(combined), size = nrow(combined), replace = T)]
      newest = summary(lm(MetaModel, data = resamples))$r.squared - summary(lm(MetaModelNULL, data = resamples))$r.squared
      
      incR2_resamples = append(incR2_resamples, newest)
    }
    
    incR2_resamples = sort(incR2_resamples[1:1000])
    returndd$incR2_CI_L = incR2_resamples[25]
    returndd$incR2_CI_U = incR2_resamples[975]
    returndd$incR2_bootstrapest = list(incR2_resamples)
    
    returndd$incAUC_CI_L = NA
    returndd$incAUC_CI_U = NA
    returndd$incAUC_bootstrapest = NA
  }
  
  
  returndd = cbind(phenolist[ii,], returndd)
  
  return(returndd)
}










