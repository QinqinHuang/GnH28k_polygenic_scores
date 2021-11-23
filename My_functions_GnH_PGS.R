#--------------------------------------------------
# 2021-11-23
#
# This R script defines functions that I used in
# the analysis.
#
#
# 1. get a list of unrelated samples while keep
# the most samples
# //remove_rel()
#  "related_pairs": ID1, ID2
#     1st and 2nd degree related pairs
#     will merge with the sample list
#  "samplelist": ID
# 
# 2. get a list of unrelated samples; keeping
# the most cases
# //remove_rel_pheno()
#  "related_pairs": ID1, ID2
#     1st and 2nd degree related pairs
#     will merge with pheno data
#  "phenolist": ID, pheno (0:control, 1:case)
#
#
# 3. Load PGS file (output of PRSice2) given the
# file name
# //loadPGSfile()
# "PGSfilename"
# 
# 4. Load the best PGS given C+T parameters
# //load_bestPGS()
# "PGSfilename": default ""; if given, load data from this file
# "ifonUKSerp": TURE - load G&H, FALSE - load eMerge; 1kG EUR as LD
# "whichGWASsumm"
# "LDr2": 0.1, 0.2, 0.5, or 0.8
# "pval"
# "scalePGS": TRUE or FALSE
#
#
# 6. A function to estimate association statistics for binary 
# (incAUC) or quantitative traits (incR2)
# //getstatistics() 
# "mydd": data with pheno, all covariates, and PRS 
# or PRS_scaled
# "binary": TRUE or FALSE
# "quandata": default "raw", or "IRN"
# "scale_PRS": default TRUE 
# "bootstrapping": default "yes_for_siginicant", or "no", "yes"
# "fullmodel": default defined 
# "nullmodel": default defined
#
#
# 7. A function to  
#  (1) calculate incR2 for each C+T PGS
#  (2) select the best C+T PGS - incR2 for cont and incAUC for binary
#  (3) calculate boostrapping 95% CI for the best C+T PGS
#  Doesn't return any value; save two output files in 
# the current working dir.
# //validate_PGSs
# "phenotable": must have ID, cov, and trait as one colname
# "cohort": "ELGH" or "eMergeEUR"
# "trait"
# "GWAS"
# "stratifiedBYsex": default "no", or "male", "female" 
# "binary": default TRUE
# "quandata": default "raw", or "IRN"
# "whichage": default "at_recrt", or "at_event"
# "PGSfolder": default NA
#    if not provided, 
#    eMergeEUR: "~/storage/PRS_ELGH/Polygenic_scores/eMerge_1kGEUR_LD"
#    ELGH: "~/P_drive/Polygenic_scores_Feb28k/all28k"
# "LDthreshold_list": default c(0.1, 0.2, 0.5, 0.8)
#--------------------------------------------------
library(data.table)
library(foreach)
library(pROC)
library(RNOmni)


#------ Association models ------
nPCs = 10
cat("** Correcting for", nPCs, "PCs **\n")

fullmodel = formula(paste0("pheno ~ PRS_scaled + gender + age+I(age^2) +", 
                           paste0("PC",1:nPCs,collapse = "+")))
nullmodel = formula(paste0("pheno ~ gender + age+I(age^2) +", 
                           paste0("PC",1:nPCs,collapse = "+")))

# sex-specific analysis
fullmodel_wogender = formula(paste0("pheno ~ PRS_scaled + age+I(age^2) +", 
                                    paste0("PC",1:nPCs,collapse = "+")))
nullmodel_wogender = formula(paste0("pheno ~ age+I(age^2) +", 
                                    paste0("PC",1:nPCs,collapse = "+")))



#----- Load phenotype, cov data and pairs of related individuals -----
### "phenotable_2ndexcl": phenotype and covariate data for unrelateds 
# samples in rows and variables in columns

### "related_pairs": rows are pairs of related samples, first two columns are sample IDs
#-----------------------------------


# 1. get a list of unrelated samples
remove_rel = function(related_pairs, samplelist) {
  # relevant pairs
  related_pairs = related_pairs[(ID1 %in% samplelist$ID) & (ID2 %in% samplelist$ID)]
  
  cat("  Number of related pairs:", nrow(related_pairs), "\n")
  cat("   Number of unique samples:", length(c(unique(related_pairs$ID2), unique(related_pairs$ID1))), "; unique samples in col1 and col2:", length(unique(related_pairs$ID1)), length(unique(related_pairs$ID2)), "\n")
  
  # For paris with the same phenotype, remove the one with the most relatives each time
  # an empty vector to keep the list of individuals to be removed
  tbrm = c()
  # aim is to make this table empty
  leftpihat = related_pairs
  
  while(nrow(leftpihat) > 0) {
    # get the most frequent one
    rmthis = names(sort(table(c(leftpihat$ID1, leftpihat$ID2)), decreasing = T))[1]
    tbrm = c(tbrm, rmthis)
    leftpihat = leftpihat[which((!ID1 %in% rmthis) & (!ID2 %in% rmthis))]
  }
  
  cat("  ** Number of indivudals to be removed:", length(unique(tbrm)), "\n")
  samplelist = samplelist[!ID %in% tbrm]
  
  return(samplelist)
}


# 2. get a list of unrelated samples, exclude controls preferentially
remove_rel_pheno = function(related_pairs, phenolist) {
  # relevant pairs
  related_pairs = related_pairs[(ID1 %in% phenolist$ID) & (ID2 %in% phenolist$ID)]
  
  cat("  Number of related pairs:", nrow(related_pairs), "\n")
  cat("   Number of unique samples:", length(c(unique(related_pairs$ID2), unique(related_pairs$ID1))), "; unique samples in col1 and col2:", length(unique(related_pairs$ID1)), length(unique(related_pairs$ID2)), "\n")
  
  # For paris with the same phenotype, remove the one with the most relatives each time; if multiple samples have the same number, remove Control if any
  tbrm = c()
  # aim is to make this table empty
  leftpairs = related_pairs
  while(nrow(leftpairs) > 0) {
    temptablesumm = sort(table(c(leftpairs$ID1, leftpairs$ID2)), decreasing = T)
    temptablesumm = temptablesumm[which(temptablesumm == max(temptablesumm))]
    if(length(temptablesumm) == 1) {
      rmthis = names(temptablesumm)
    } else {
      temptablesumm = data.table(ID = names(temptablesumm))
      temptablesumm$pheno = phenolist[match(temptablesumm$ID, phenolist$ID)]$pheno
      rmthis = ifelse(0 %in% temptablesumm$pheno, 
                      yes = temptablesumm[pheno==0]$ID[1],
                      no = temptablesumm$ID[1])
    }
    tbrm = c(tbrm, rmthis)
    leftpairs = leftpairs[which((!ID1 %in% rmthis) & (!ID2 %in% rmthis))]
  }
  
  cat("  ** Number of indivudals to be removed:", length(unique(tbrm)))
  phenolist = phenolist[!ID %in% tbrm]
  print(table(phenolist$pheno))
  cat("   Prevalence:", sum(phenolist$pheno == 1)/nrow(phenolist), "\n")
  return(phenolist)
}


# 3. Load PGS file for a given file name
loadPGSfile = function(PGSfilename) {
  if(!file.exists(PGSfilename)) {
    cat(" PGS file doesn't exist!\n")
    return(NULL)
  }
  
  # skip the header line
  allscores = fread(PGSfilename, header = F, skip = 1)
  
  # column name
  header_allscores = fread(PGSfilename, header = F, nrows = 1)
  
  if(ncol(header_allscores) == ncol(allscores)+1) {
    names(allscores) = c("ID", unlist(header_allscores[1,-1:-2]))
  } else if(ncol(header_allscores) == ncol(allscores) & unlist(header_allscores)[1] == "FID" & unlist(header_allscores)[2] == "IID") {
    names(allscores) = unlist(header_allscores)
    allscores = allscores[, -2]
    names(allscores)[1] = "ID"
  } else {
    stop("  ** check number of scores in the score file!")
  }
  
  return(allscores)
}


# 4. Load the best PGS for given C+T parameters
load_bestPGS = function(PGSfilename = "", whichGWASsumm, ifonUKSerp,
                        LDr2, pval, scalePGS = F) {
  
  # Read all PGSs
  if(PGSfilename == "") {
    if(ifonUKSerp) {
      
      ### you will need to update the file path
      PGSfilename = paste0("~/P_drive/Polygenic_scores_Feb28k/all28k/",whichGWASsumm,"_GAsP_bgen_dosage/PRSice_output_",LDr2,".all.score")
    } else {
      PGSfilename = paste0("~/storage/PRS_ELGH/Polygenic_scores/eMerge_1kGEUR_LD/",whichGWASsumm,"_bgen_dosage/PRSice_output_",LDr2,".all.score")
    }
  }
  cat("  PGS file:", PGSfilename, "\n")
  
  # read the PGS file
  allscores = loadPGSfile(PGSfilename)
  
  
  pval_num = names(allscores)[which(as.numeric(names(allscores)[-1]) == pval) + 1]
  cat("   LD clumping r2", LDr2, "and P-value threshold is", pval_num, "\n")
  bestPGS = data.table(ID = allscores$ID, PRS = allscores[, get(pval_num)])
  
  if(scalePGS) {
    bestPGS$PRS_scaled = scale(bestPGS$PRS)
  }
  
  return(bestPGS)
}



# 6. calcualte association statistics for binary or quantitative traits
getstatistics = function(mydd, binary, quandata = "raw", 
                         scale_PRS = T, bootstrapping = "yes_for_siginicant",
                         fullmodel = fullmodel, nullmodel = nullmodel) {
  
  # Remove missing phenotype data
  mydd = mydd[which(!is.na(mydd$pheno))]
  
  
  # quantitative data might need to be inverse-normal transformed
  if(!binary & quandata == "IRN") {
    mydd$pheno = rankNorm(mydd$pheno)
    cat(" Using inverse transformed data!\n")
  }
  
  # scale PRS
  if(scale_PRS) {
    mydd$PRS_scaled = scale(mydd$PRS)
  }
  
  # Binary traits
  if(binary) {
    # Logistic regression
    myfit_full = glm(fullmodel, data = mydd, family = binomial())
    myfit_null = glm(nullmodel, data = mydd, family = binomial())
    
    returndd = data.table(N = nrow(mydd),
                          prevalence = mean(mydd$pheno),
                          Ncases = sum(mydd$pheno),
                          AUC_full = roc(mydd$pheno, predict(myfit_full, type = c("response")))$auc[1],
                          AUC_null = roc(mydd$pheno, predict(myfit_null, type = c("response")))$auc[1],
                          PRS_beta = summary(myfit_full)$coefficients["PRS_scaled",1],
                          PRS_se = summary(myfit_full)$coefficients["PRS_scaled",2],
                          PRS_p = summary(myfit_full)$coefficients["PRS_scaled",4],
                          full_R2 = NA, null_R2 = NA )
    
    returndd[, incR2 := NA]
    returndd[, incAUC := AUC_full - AUC_null]
  }
  
  # Quantitative triats
  if(!binary) {
    # Linear regression
    myfit_full = summary(lm(fullmodel, data = mydd))
    myfit_null = summary(lm(nullmodel, data = mydd))
    
    returndd = data.table(N = nrow(mydd),
                          prevalence = NA,
                          Ncases = NA,
                          AUC_full = NA, 
                          AUC_null = NA,
                          PRS_beta = myfit_full$coefficients["PRS_scaled",1],
                          PRS_se = myfit_full$coefficients["PRS_scaled",2],
                          PRS_p = myfit_full$coefficients["PRS_scaled",4],
                          full_R2 = myfit_full$r.squared,
                          null_R2 = myfit_null$r.squared)
    
    returndd[, incR2 := full_R2 - null_R2]
    returndd[, incAUC := NA]
  }
  
  # run 1000 bootstraps to estimate 95% CI for incAUC for binary traits and incR2 for quantitative traits
  runb = F
  if(bootstrapping == "yes") {runb = T}
  if(bootstrapping == "yes_for_siginicant" & returndd$PRS_p < 0.05) {runb = T}
  
  if(runb) {
    # (1) quantitative trait - incR2
    if(!binary) {
      incR2_resamples = vector()
      while(length(incR2_resamples) <1000) {
        resamples = mydd[sample(1:nrow(mydd), size = nrow(mydd), replace = T)]
        newest = summary(lm(fullmodel, data = resamples))$r.squared - summary(lm(nullmodel, data = resamples))$r.squared
        incR2_resamples = append(incR2_resamples, newest)
      }
      
      incR2_resamples = sort(incR2_resamples[1:1000])
      returndd$incR2_CI_L = incR2_resamples[25]
      returndd$incR2_CI_U = incR2_resamples[975]
      returndd$incR2_bootstrapest = list(incR2_resamples)
      
      returndd$incAUC_CI_L = NA
      returndd$incAUC_CI_U = NA
      returndd$incAUC_bootstrapest = NA
    } # end of quantitative trait
    
    
    # (2) binary trait - incAUC
    if(binary) {
      incAUC_resamples = vector()
      while(length(incAUC_resamples) <1000) {
        resamples = mydd[sample(1:nrow(mydd), size = nrow(mydd), replace = T)]
        myfit_full = glm(fullmodel, data = resamples, family = binomial())
        myfit_null = glm(nullmodel, data = resamples, family = binomial())
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
    } # end of binary trait
    
  } else {
    returndd$incR2_CI_L = NA
    returndd$incR2_CI_U = NA
    returndd$incR2_bootstrapest = NA
    returndd$incAUC_CI_L = NA
    returndd$incAUC_CI_U = NA
    returndd$incAUC_bootstrapest = NA
  } # end of bootstrap
  
  
  return(returndd)
}



# 7. A function to  (1) calculate incR2 for each C+T PGS, (2) select the best C+T PGS, and (3) calculate boostrapping 95% CI for the best C+T PGS
validate_PGSs = function(phenotable, cohort, trait, GWAS, 
                         stratifiedBYsex = "no",
                         binary = T, quandata = "raw", whichage = "at_recrt",
                         LDthreshold_list = c(0.1, 0.2, 0.5, 0.8),
                         PGSfolder = NA) {
  
  #---- Name of input/output files ----
  # eMerge data are on farm5
  if(cohort == "eMergeEUR") {
    
    ### you will need to update the file path
    if(is.na(PGSfolder)) { PGSfolder = "~/storage/PRS_ELGH/Polygenic_scores/eMerge_1kGEUR_LD" }
    PGSpath = paste0(PGSfolder, "/", GWAS, "_bgen_dosage")
    
    # Null model output
    filename_null = paste0(trait,"_",GWAS,"_Null_model_",ifelse(stratifiedBYsex=="no",yes="",no=paste0(stratifiedBYsex,"_")),quandata,"_age_",whichage,".txt")
    
    # Full model output
    filename_full = paste0(trait,"_",GWAS,"_Full_model_",ifelse(stratifiedBYsex=="no",yes="",no=paste0(stratifiedBYsex,"_")),quandata,"_age_",whichage,".RDS")
  }
  
  # ELGH data are on UKSerp
  if(cohort == "ELGH") {
    
    ### you will need to update the file path
    if(is.na(PGSfolder)) { PGSfolder = "~/P_drive/Polygenic_scores_Feb28k/all28k" }
    PGSpath = paste0(PGSfolder, "/", GWAS, "_GAsP_bgen_dosage")
    
    # Null model output
    filename_null = paste0(trait,"_",cohort,"_",GWAS,"_Null_model_",ifelse(stratifiedBYsex=="no",yes="",no=paste0(stratifiedBYsex,"_")),quandata,"_age_",whichage,".txt")
    
    # Full model output
    filename_full = paste0(trait,"_",cohort,"_",GWAS,"_Full_model_",ifelse(stratifiedBYsex=="no",yes="",no=paste0(stratifiedBYsex,"_")),quandata,"_age_",whichage,".RDS")
  }
  
  
  #----- check if the polygenic scores have been calculated -----
  if(!file.exists(PGSpath)) {
    cat(" ****PGS haven't been calculated.\n")
    return(NULL)
  }
  
  
  #----- prepare phenotype table -----
  # remove samples with missing phenotype data
  mypheno = phenotable
  mypheno$pheno = mypheno[, get(trait)]
  mypheno = mypheno[!is.na(pheno)]
  
  # If focus on males or females only
  if(stratifiedBYsex %in% c("female", "male")) {
    
    # gender is coded as 1 and 2
    if(stratifiedBYsex == "male") {
      mypheno = mypheno[gender == 1]
    } else if(stratifiedBYsex == "female") {
      mypheno = mypheno[gender == 2]
    }
    
    # update models
    fullmodel = fullmodel_wogender
    nullmodel = nullmodel_wogender
    
    cat("  ** Focusing on", stratifiedBYsex, "\n")
  } #end of "stratifiedBYsex"
  
  
  # quantitative data might need to be inverse-normal transformed
  if(!binary) {
    if(quandata == "IRN") {
      mypheno$pheno = rankNorm(mypheno[, get(trait)])
      cat(" Using inverse transformed data!\n")
    } else if(quandata == "raw") {
      cat(" Using raw data!\n")
    }
  }
  
  # Age: age at recruitment or age at measurement
  if(whichage == "at_recrt") {
    mypheno[, age := age_at_recrt]
  } else if(whichage == "at_event") {
    cc = paste0("age_", trait)
    mypheno[, age := get(cc)]
  } else if(whichage == "age2020") {
    # in eMERGE age was calculated as 2020-yob
    if(!"age" %in% names(mypheno)) { stop("  eMerge data do not have age!") }
  }
  
  
  
  # the final table
  if(nPCs == 20) {
    mydata = mypheno[,.(ID, pheno, age, gender, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                       PC11, PC12, PC13, PC14, PC15, PC16, PC17, PC18, PC19, PC20)]
  } else {
    mydata = mypheno[,.(ID, pheno, age, gender, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]
  }
  
  
  
  #----- Null model without PGS -----
  # Binary traits: logistic regression
  if(binary) {
    myfit_null = glm(nullmodel, data = mydata, family = binomial())
    nullstats = summary(myfit_null)$coefficients[-1,]
    nullstats_dt = as.data.table(nullstats)
    names(nullstats_dt) = c("beta","se","zvalue","pval")
    nullstats_dt$predictor = rownames(nullstats)
    nullstats_dt = nullstats_dt[,.(predictor, beta, se, zvalue, pval)]
    # AUC of the null model
    nullstats_dt$AUC = roc(mydata$pheno, predict(myfit_null, type = c("response")))$auc[1]
  } else if(!binary) {  # Quantitative traits: linear regression
    myfit_null = summary(lm(nullmodel, data = mydata))
    nullstats = myfit_null$coefficients[-1,]
    nullstats_dt = as.data.table(nullstats)
    names(nullstats_dt) = c("beta","se","tvalue","pval")
    nullstats_dt$predictor = rownames(nullstats)
    nullstats_dt = nullstats_dt[,.(predictor, beta, se, tvalue, pval)]
  }
  nullstats_dt[, `:=` (ancestry = cohort, phenotype = trait, transformation = quandata,  
                       GWASsumm = GWAS, agedata = whichage)]
  
  # Save the summary statistics
  fwrite(nullstats_dt, file = filename_null, sep = "\t")
  
  
  #----- Full model with C+T PGS -----
  # Four LD clumping thresholds
  asso4LDr2 = foreach(LDthreshold = LDthreshold_list, .combine = rbind) %do% {
    
    # Load all scores at various p-value thresholds
    PGSfilename = paste0(PGSpath,"/PRSice_output_",LDthreshold,".all.score")
    if(!file.exists(PGSfilename)) {return(NULL)}
    
    allscores = loadPGSfile(PGSfilename = PGSfilename)
    
    # All p-value thresholds
    allpvalues = names(allscores)[-1]
    
    
    # Merge phenotype data and PGSs
    combined = merge(mydata, allscores, by = "ID")
    
    
    #--- Validate all scores ---
    allstats = foreach(ii = allpvalues, .combine = rbind) %do% {
      
      # Extract the phenotype and the score
      cc = c(names(mydata), ii)
      mydd = combined[, ..cc]
      
      # scale PRS
      mydd$PRS_scaled = scale(mydd[, get(ii)])
      
      # (1) Binary traits: logistic regression
      if(binary) {
        # logistic regression
        myfit_full = glm(fullmodel, data = mydd, family = binomial())
        
        # summary statistics of the full model, AUC and R2
        returndd = data.table(N = nrow(mydd),
                              prevalence = sum(mydd$pheno==1)/nrow(mydd),
                              LDr2 = LDthreshold,
                              pval_threshold = ii,
                              AUC_full = roc(mydd$pheno, predict(myfit_full, type = c("response")))$auc[1],
                              AUC_null = roc(mydd$pheno, predict(myfit_null, type = c("response")))$auc[1],
                              PRS_beta = summary(myfit_full)$coefficients["PRS_scaled",1],
                              PRS_se = summary(myfit_full)$coefficients["PRS_scaled",2],
                              PRS_p = summary(myfit_full)$coefficients["PRS_scaled",4],
                              full_R2 = NA,
                              null_R2 = NA )
        
        returndd[, incR2 := NA]
        returndd[, incAUC := AUC_full - AUC_null]
      }
      
      # (2) Quantitative traits: linear regression
      if(!binary) {
        # Linear regression
        myfit_full = summary(lm(fullmodel, data = mydd))
        
        # summary statistics of the full model
        returndd = data.table(N = nrow(mydd), prevalence = NA, 
                              LDr2 = LDthreshold,
                              pval_threshold = ii, AUC_full = NA, AUC_null = NA,
                              PRS_beta = myfit_full$coefficients["PRS_scaled",1],
                              PRS_se = myfit_full$coefficients["PRS_scaled",2],
                              PRS_p = myfit_full$coefficients["PRS_scaled",4],
                              full_R2 = myfit_full$r.squared,
                              null_R2 = myfit_null$r.squared )
        
        returndd[, incR2 := full_R2 - null_R2]
        returndd[, incAUC := NA]
      }
      
      return(returndd)
    }
    
    return(allstats)
  }
  
  # Label the best score - incR2 for quantitative and incAUC for binary
  asso4LDr2$bestscore = ""
  if(!binary) {
    asso4LDr2$bestscore[asso4LDr2$incR2 == max(asso4LDr2$incR2)] = "BestScore"
  } else {
    asso4LDr2$bestscore[asso4LDr2$incAUC == max(asso4LDr2$incAUC)] = "BestScore"
  }
  
  
  
  #----- Bootstrap 95% CI for the best PGS -----
  # P-value and LD r2 threshold of the best score
  bestone = asso4LDr2[bestscore == "BestScore"][order(PRS_p)][1,]
  best_pval = bestone$pval_threshold
  best_LDthreshold = bestone$LDr2
  
  # read the best PRS 
  bestPGS = load_bestPGS(PGSfilename = paste0(PGSpath,"/PRSice_output_",best_LDthreshold,".all.score"),
                         LDr2 = best_LDthreshold, pval = best_pval, scalePGS = F)
  
  
  # Merge phenotype data with all scores
  bestPGS = merge(mydata, bestPGS, by = "ID")
  
  # scale PRS
  bestPGS$PRS_scaled = scale(bestPGS$PRS)
  
  # 1000 bootstrap resamples
  # (1) quantitative traits
  if(!binary) {
    incR2_resamples = vector()
    while(length(incR2_resamples) <1000) {
      resamples = bestPGS[sample(1:nrow(bestPGS), size = nrow(bestPGS), replace = T)]
      newest = summary(lm(fullmodel, data = resamples))$r.squared - summary(lm(nullmodel, data = resamples))$r.squared
      incR2_resamples = append(incR2_resamples, newest)
    }
    
    incR2_resamples = sort(incR2_resamples[1:1000])
    asso4LDr2$incR2_CI_L = NA
    asso4LDr2$incR2_CI_L[which(asso4LDr2$pval_threshold == best_pval & 
                                 asso4LDr2$LDr2 == best_LDthreshold)] = incR2_resamples[25]
    asso4LDr2$incR2_CI_U = NA
    asso4LDr2$incR2_CI_U[which(asso4LDr2$pval_threshold == best_pval & 
                                 asso4LDr2$LDr2 == best_LDthreshold)] = incR2_resamples[975]
    
    # return 1000 bootstrap estimates too
    asso4LDr2$incR2_bootstrapest = NA
    asso4LDr2$incR2_bootstrapest[which(asso4LDr2$pval_threshold == best_pval & 
                                   asso4LDr2$LDr2 == best_LDthreshold)] = list(incR2_resamples)
    
    asso4LDr2$incAUC_CI_L = NA
    asso4LDr2$incAUC_CI_U = NA
    asso4LDr2$incAUC_bootstrapest = NA
  }# end of bootstrap for quantitative 
  
  
  # (2) bianry traits - incAUC
  if(binary) {
    incAUC_resamples = vector()
    while(length(incAUC_resamples) <1000) {
      resamples = bestPGS[sample(1:nrow(bestPGS), size = nrow(bestPGS), replace = T)]
      myfit_full = glm(fullmodel, data = resamples, family = binomial())
      myfit_null = glm(nullmodel, data = resamples, family = binomial())
      myAUCfull = roc(resamples$pheno, predict(myfit_full, type = c("response")))$auc[1]
      myAUCnull = roc(resamples$pheno, predict(myfit_null, type = c("response")))$auc[1]
      newest = myAUCfull - myAUCnull
      incAUC_resamples = append(incAUC_resamples, newest)
    }
    
    incAUC_resamples = sort(incAUC_resamples[1:1000])
    asso4LDr2$incAUC_CI_L = NA
    asso4LDr2$incAUC_CI_L[which(asso4LDr2$pval_threshold == best_pval & 
                                 asso4LDr2$LDr2 == best_LDthreshold)] = incAUC_resamples[25]
    asso4LDr2$incAUC_CI_U = NA
    asso4LDr2$incAUC_CI_U[which(asso4LDr2$pval_threshold == best_pval & 
                                 asso4LDr2$LDr2 == best_LDthreshold)] = incAUC_resamples[975]
    
    # return 1000 bootstrap estimates too
    asso4LDr2$incAUC_bootstrapest = NA
    asso4LDr2$incAUC_bootstrapest[which(asso4LDr2$pval_threshold == best_pval & 
                                   asso4LDr2$LDr2 == best_LDthreshold)] = list(incAUC_resamples)
    
    asso4LDr2$incR2_CI_L = NA
    asso4LDr2$incR2_CI_U = NA
    asso4LDr2$incR2_bootstrapest = NA
  }# end of bootstrap for binary 
  
  
  # save the summary statistics
  asso4LDr2[, `:=` (ancestry = cohort, phenotype = trait, transformation = quandata, GWASsumm = GWAS, agedata = whichage)]
  
  saveRDS(asso4LDr2, file = filename_full)
}



