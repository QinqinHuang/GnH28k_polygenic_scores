#------------------------------------------------------
# 2021-06-01
# 
# Aim: 
#   to calculate the integrated score
#   to run Cox regression model
#   to estimate categorical NRI
#   to estimate continuous NRI
#   
# QRISK3 assessment date:
AssessmentDate = as.POSIXct("2010-01-01", format = "%Y-%m-%d")
#------------------------------------------------------
### Functions are defined in: https://github.com/QinqinHuang/GnH28k_polygenic_scores/blob/main/My_functions_GnH_PGS.R
source("My_functions_GnH_PGS.R")

library(survival)


#----- the best PGS for CAD -----
# Load the best score from PGS Catalog (Amit's LDpred score: GPS_CAD_SA, PGS000296)
myPGSID = "PGS000296"
# read the PGS
myPGS = fread(paste0("~/P_drive/Polygenic_scores_Feb28k/PGSCatalog/ELGH_PGS/plink2_score_sum_",myPGSID,".txt"))
# scale the PGS
myPGS$PRS_scaled = scale(myPGS$score_sum)

# regress out 10 PCs
covtable = fread("~/P_drive/ELGH_phenotype_2020Feb/clean_pheno_cov_table/Covariate_table.txt")
CADPRS = merge(myPGS[,.(ID, PRS_scaled)], 
               covtable[, .(ID = paste0("1_",ID), PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)], by = "ID")

fitPCs = summary(lm(PRS_scaled ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = CADPRS))
# scaling the residuals
CADPRS[, CADPRScorr := scale(fitPCs$residuals)]

names(CADPRS)[which(names(CADPRS) == "PRS_scaled")] = "CADPRS"


#----- read QIRSK3 scores -----
QRISK3input = readRDS("../QRISK3_inputdata_2010-01-01__nonmissing3_9477_impute_pmm.RDS")
QRISK3score = fread("../QRISK3_calculated_2010-01-01__nonmissing3_9477_impute_pmm.txt")

identical(QRISK3input$pseudoNHSnumber, QRISK3score$pseudoNHSnumber)
CADpheno = QRISK3input[, .(ID, age_at_recrt, gender, CAD, age_at_QRISK3 = age, CAD_Dx_DateE,
                           QRISK3 = QRISK3score$QRISK3_2017)]

# merge the data with PGS
CADpheno = merge(CADpheno, CADPRS, by = "ID")


#----- keep only unrelated individuals -----
dd = remove_rel_pheno(related_pairs, CADpheno[, .(ID, pheno = CAD)])
CADpheno = CADpheno[ID %in% dd$ID]
# before removing relatives: Ncontrol=9039, Ncases=438
# after: Ncases=420, Ncontrols=7702

# correlation between PRS and QRISK3
cor.test(CADpheno$CADPRScorr, CADpheno$QRISK3)


#----- Calculate the integrated score -----
# Odds of QIRSK3
CADpheno[, OddsQRISK3 := QRISK3/100 / (1-QRISK3/100)]
# logit-transform of QRISK3
CADpheno[, logitQRISK3 := log(OddsQRISK3)]

# sex-specific OR for the CAD PGS
myfit = summary(glm(CAD ~ CADPRScorr + logitQRISK3 + CADPRScorr*logitQRISK3, CADpheno[gender == 0], family = binomial()))
PRScorrbeta_male = myfit$coefficient["CADPRScorr", 1]
intbeta_male = myfit$coefficient["CADPRScorr:logitQRISK3", 1]

myfit = summary(glm(CAD ~ CADPRScorr + logitQRISK3 + CADPRScorr*logitQRISK3, CADpheno[gender == 1], family = binomial()))
PRScorrbeta_female = myfit$coefficient["CADPRScorr", 1]
intbeta_female = myfit$coefficient["CADPRScorr:logitQRISK3", 1]

# calculate IRT 
CADpheno[gender == 0, OddsIRT_sexOR := OddsQRISK3 * exp(PRScorrbeta_male*CADPRScorr + intbeta_male*CADPRScorr*logitQRISK3)]
CADpheno[gender == 1, OddsIRT_sexOR := OddsQRISK3 * exp(PRScorrbeta_female*CADPRScorr + intbeta_female*CADPRScorr*logitQRISK3)]
CADpheno[, IRT_sexOR := OddsIRT_sexOR/(1+OddsIRT_sexOR)*100]



#----- Cox regression and categorcial NRI -----
# time: from assessment to diagnosis or censored
#   case: CAD_Dx_DateE - AssessmentDate
#   control: 2020-02 - AssessmentDate
# status: 1=event, 0=censored
CADpheno[, timeobserved := ifelse(CAD == 1, yes = CAD_Dx_DateE - AssessmentDate,
                                  no = as.POSIXct("2020-02-01", format = "%Y-%m-%d") - AssessmentDate)]


#--- A function to calculate cindex, its se and NRI ---
calculate_cindex_NRI = function(mydd, highriskcutoff = 10, IRTcoln, label = NA, nbootstrap = 1000) {
  mydd[, newscore := mydd[, get(IRTcoln)]]
  
  # Cox regression model to estimate C-index
  # baseline - age and gender
  cox_age_gender = summary(coxph(Surv(timeobserved, CAD) ~ age_at_QRISK3 + gender, 
                                 data = mydd))
  
  # PRS, age, and gender
  cox_PRS_age_gender = summary(coxph(Surv(timeobserved, CAD) ~ CADPRScorr + age_at_QRISK3 + gender, 
                                 data = mydd))
  
  # QRISK3
  cox_QRISK3 = summary(coxph(Surv(timeobserved, CAD) ~ QRISK3, 
                             data = mydd))
  
  # IRT: QRISK3 enhanced with PRS
  cox_IRT = summary(coxph(Surv(timeobserved, CAD) ~ newscore, data = mydd))

  # NRI for cases
  NRI_case = (length(which(mydd$QRISK3 < highriskcutoff & mydd$newscore >= highriskcutoff & mydd$CAD == 1)) - length(which(mydd$QRISK3 >= 10 & mydd$newscore < 10 & mydd$CAD == 1))) / length(which(mydd$CAD == 1))
  
  # NRI for controls
  NRI_ctrl = (length(which(mydd$QRISK3 >= highriskcutoff & mydd$newscore < highriskcutoff & mydd$CAD == 0)) - length(which(mydd$QRISK3 < highriskcutoff & mydd$newscore >= highriskcutoff & mydd$CAD == 0))) / length(which(mydd$CAD == 0))
  
  # bootstrap to estimate CI for NRI
  bootstrapest = foreach(ii = 1:1000, .combine = rbind) %do% {
    newdd = mydd[sample(1:nrow(mydd), nrow(mydd), replace = T)]
    NRI_case_new = (length(which(newdd$QRISK3 < highriskcutoff & newdd$newscore >= highriskcutoff & newdd$CAD == 1)) - length(which(newdd$QRISK3 >= 10 & newdd$newscore < 10 & newdd$CAD == 1))) / length(which(newdd$CAD == 1))
    NRI_ctrl_new = (length(which(newdd$QRISK3 >= highriskcutoff & newdd$newscore < highriskcutoff & newdd$CAD == 0)) - length(which(newdd$QRISK3 < highriskcutoff & newdd$newscore >= highriskcutoff & newdd$CAD == 0))) / length(which(newdd$CAD == 0))
    return(data.table(NRI = NRI_case_new+NRI_ctrl_new, NRI_case = NRI_case_new, 
                      NRI_ctrl = NRI_ctrl_new))
  }
  
  returndd = data.table(dataset = label, N = nrow(mydd), Ncases = nrow(mydd[CAD==1]),
                        cutoff = highriskcutoff, 
                        C_agesex = cox_age_gender$concordance[1],
                        Cse_agesex = cox_age_gender$concordance[2],
                        C_PRSagesex = cox_PRS_age_gender$concordance[1],
                        Cse_PRSagesex = cox_PRS_age_gender$concordance[2],
                        C_QRISK3 = cox_QRISK3$concordance[1],
                        Cse_QRISK3 = cox_QRISK3$concordance[2],
                        C_IRT = cox_IRT$concordance[1],
                        Cse_IRT = cox_IRT$concordance[2],
                        NRI = NRI_case+NRI_ctrl, NRI_CIL = sort(bootstrapest$NRI)[25],
                        NRI_CIU = sort(bootstrapest$NRI)[975],
                        NRI_case = NRI_case, NRI_case_CIL = sort(bootstrapest$NRI_case)[25],
                        NRI_case_CIU = sort(bootstrapest$NRI_case)[975],
                        NRI_control = NRI_ctrl, NRI_control_CIL = sort(bootstrapest$NRI_ctrl)[25],
                        NRI_control_CIU = sort(bootstrapest$NRI_ctrl)[975])
  
  return(returndd)
}


# using sex-specific OR, age cutoff 55
agecutoff = 55
NRIstats = rbind(
  # full
  calculate_cindex_NRI(mydd = CADpheno, IRTcoln = "IRT_sexOR", label = "full"),
  # young
  calculate_cindex_NRI(mydd = CADpheno[age_at_QRISK3 < agecutoff],  IRTcoln = "IRT_sexOR", label = "young"),
  # young male
  calculate_cindex_NRI(mydd = CADpheno[age_at_QRISK3 < agecutoff & gender == 0], IRTcoln = "IRT_sexOR", label = "young_male"),
  # young female
  calculate_cindex_NRI(mydd = CADpheno[age_at_QRISK3 < agecutoff & gender == 1], IRTcoln = "IRT_sexOR", label = "young_female"),
  # old
  calculate_cindex_NRI(mydd = CADpheno[age_at_QRISK3 >= agecutoff], IRTcoln = "IRT_sexOR", label = "old"),
  # old male
  calculate_cindex_NRI(mydd = CADpheno[age_at_QRISK3 >= agecutoff & gender == 0], IRTcoln = "IRT_sexOR", label = "old_male"),
  # old female
  calculate_cindex_NRI(mydd = CADpheno[age_at_QRISK3 >= agecutoff & gender == 1], IRTcoln = "IRT_sexOR", label = "old_female")
)

write.csv(NRIstats, "Cindex_NRI_full_subgroups_impute_ppm_sexOR.csv", row.names = F)



#----- get reclassification numbers for Table S12 -----
getReNum = function(mydata, highriskcutoff = 10) {
  return(
    data.table(N = nrow(mydata), Ncase = nrow(mydata[CAD==1]), Ncontrol = nrow(mydata[CAD==0]),
               N_QRISK3_atrisk = nrow(mydata[QRISK3 >= 10]),
               N_IRT_atrisk = nrow(mydata[IRT_sexOR >= 10]),
               N_IRT_up = nrow(mydata[QRISK3 < 10 & IRT_sexOR >= 10]),
               N_IRT_down = nrow(mydata[QRISK3 > 10 & IRT_sexOR <= 10]),
               N_IRT_up_case = nrow(mydata[QRISK3 < 10 & IRT_sexOR >= 10 & CAD == 1]),
               N_IRT_down_control = nrow(mydata[QRISK3 > 10 & IRT_sexOR <= 10 & CAD == 0]))
  )
}

d = rbind(
  getReNum(CADpheno),
  getReNum(CADpheno[age_at_QRISK3 < 55]),
  getReNum(CADpheno[age_at_QRISK3 < 55 & gender == 0]),
  getReNum(CADpheno[age_at_QRISK3 < 55 & gender == 1]),
  getReNum(CADpheno[age_at_QRISK3 >= 55]),
  getReNum(CADpheno[age_at_QRISK3 >= 55 & gender == 0]),
  getReNum(CADpheno[age_at_QRISK3 >= 55 & gender == 1]) 
)
fwrite(d, "reclassification_number_impute_ppm_sexOR.csv", row.names = F)








#----- Continuous NRI -----
# the sum of proportions of cases and noncases with improved combined score 
# (ie, higher combined score for cases denoted by P[up|case] {where P indicates probability}
# and lower for noncases denoted by P[down|noncase]) minus the sum of proportions with deteriorated combined score (ie, P[up|noncase]) and P[down|case])

#  P[up|case] - P[down|case] + P[down|noncase] - P[up|noncase] 

#----- A function to calculate category-free NRI -----
calculate_contNRI = function(mydd, highriskcutoff = 10, IRTcoln, label = NA, nbootstrap = 1000) {
  mydd[, newscore := mydd[, get(IRTcoln)]]
  
  # NRI for cases
  NRI_case = (length(which(mydd$QRISK3 < mydd$newscore & mydd$CAD == 1)) - length(which(mydd$QRISK3 > mydd$newscore & mydd$CAD == 1))) / length(which(mydd$CAD == 1))
  
  # NRI for controls 
  NRI_ctrl = (length(which(mydd$QRISK3 > mydd$newscore & mydd$CAD == 0)) - length(which(mydd$QRISK3 < mydd$newscore & mydd$CAD == 0))) / length(which(mydd$CAD == 0))
  
  # bootstrap to estimate CI for NRI
  bootstrapest = foreach(ii = 1:1000, .combine = rbind) %do% {
    newdd = mydd[sample(1:nrow(mydd), nrow(mydd), replace = T)]
    NRI_case_new = (length(which(newdd$QRISK3 < newdd$newscore & newdd$CAD == 1)) - length(which(newdd$QRISK3 > newdd$newscore & newdd$CAD == 1))) / length(which(newdd$CAD == 1))
    NRI_ctrl_new = (length(which(newdd$QRISK3 > newdd$newscore & newdd$CAD == 0)) - length(which(newdd$QRISK3 < newdd$newscore & newdd$CAD == 0))) / length(which(newdd$CAD == 0))
    return(data.table(NRI = NRI_case_new+NRI_ctrl_new, NRI_case = NRI_case_new, 
                      NRI_ctrl = NRI_ctrl_new))
  }
  
  returndd = data.table(dataset = label, N = nrow(mydd), Ncases = nrow(mydd[CAD==1]),
                        NRI = NRI_case+NRI_ctrl, NRI_CIL = sort(bootstrapest$NRI)[25],
                        NRI_CIU = sort(bootstrapest$NRI)[975],
                        NRI_case = NRI_case, NRI_case_CIL = sort(bootstrapest$NRI_case)[25],
                        NRI_case_CIU = sort(bootstrapest$NRI_case)[975],
                        NRI_control = NRI_ctrl, NRI_control_CIL = sort(bootstrapest$NRI_ctrl)[25],
                        NRI_control_CIU = sort(bootstrapest$NRI_ctrl)[975])
  
  return(returndd)
}


# using sex-specific OR, age cutoff 55
agecutoff = 55
NRIstats = rbind(
  # full
  calculate_contNRI(mydd = CADpheno, IRTcoln = "IRT_sexOR", label = "full"),
  # young
  calculate_contNRI(mydd = CADpheno[age_at_QRISK3 < agecutoff],  IRTcoln = "IRT_sexOR", label = "young"),
  # young male
  calculate_contNRI(mydd = CADpheno[age_at_QRISK3 < agecutoff & gender == 0], IRTcoln = "IRT_sexOR", label = "young_male"),
  # young female
  calculate_contNRI(mydd = CADpheno[age_at_QRISK3 < agecutoff & gender == 1], IRTcoln = "IRT_sexOR", label = "young_female"),
  # old
  calculate_contNRI(mydd = CADpheno[age_at_QRISK3 >= agecutoff], IRTcoln = "IRT_sexOR", label = "old"),
  # old male
  calculate_contNRI(mydd = CADpheno[age_at_QRISK3 >= agecutoff & gender == 0], IRTcoln = "IRT_sexOR", label = "old_male"),
  # old female
  calculate_contNRI(mydd = CADpheno[age_at_QRISK3 >= agecutoff & gender == 1], IRTcoln = "IRT_sexOR", label = "old_female")
)
NRIstats[,.(dataset, NRI, NRI_case, NRI_control)]

write.csv(NRIstats, "Continuous_NRI_full_subgroups_impute_ppm_sexOR.csv", row.names = F)












