#------------------------------------------------------------
# 2021-04-30
#
# Aims: 
# To do quintile analysis of PGS
#    phenotype ~ PRSgroup + gender + age + age^2 + 10PCs
#
# for binary traits, the cut-offs are defined in the controls
#------------------------------------------------------------
### Functions are defined in: https://github.com/QinqinHuang/GnH28k_polygenic_scores/blob/main/My_functions_GnH_PGS.R
source("My_functions_GnH_PGS.R")


Modelbin = formula("pheno ~ PRS_bin_label+gender+age+I(age^2)+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10")


#----- Load data -----
### "bestPGS": this object contains the best PGS per trait
# scores are from the PGS Catalog 


#----- (1) Quantile analysis for binary traits -----
assostats = foreach(ii = 1:nrow(bestPGS), .combine = rbind) %do% {
  
  myisbinary = !is.na(bestPGS$AUC_full[ii])
  if(!myisbinary) {return(NULL)}
  
  myPGSID = bestPGS$PGS_ID[ii]
  mytrait = bestPGS$phenotype[ii]
  
  
  # read the PGS given the PGS ID
  myPGS = fread(paste0("yourpath/plink2_score_sum_",myPGSID,".txt"))
  # scale 
  myPGS$PRS_scaled = scale(myPGS$score_sum)
  
  # merge PGS to the phenotype table
  mydd = merge(phenotable_2ndexcl, myPGS[,.(ID, PRS_scaled)], by = "ID")
  
  # use age_at_recrt for diseases
  mydd[, age := age_at_recrt]
  
  # remove missing data
  mydd$pheno = mydd[, get(mytrait)]
  mydd = mydd[!is.na(pheno)]
  
  
  
  # thresholds are defined in controls
  combined = mydd
  combined[, PRS_bin_label := cut(PRS_scaled, quantile(combined[pheno==0]$PRS_scaled, probs = 0:5/5), labels = FALSE, include.lowest = TRUE)]
  combined$PRS_bin_label = factor(combined$PRS_bin_label, levels = c(3,1,2,4,5))
  
  # logistic regression
  myfit1 = glm(Modelbin, data = combined, family = binomial(link = "logit"))
  quincoef = summary(myfit1)$coefficients
  
  # OR, CI, and Pvalue
  metrics = as.data.table(quincoef[2:5, -3])
  names(metrics) = c("beta","se","p_value")
  returndd = rbind(metrics[1:2], data.table(beta = 0, se = 0, p_value = 1),
                   metrics[3:4])
  returndd[, `:=` (phenotype = mytrait, PGS_ID = myPGSID, 
                   ancestry = "ELGH", Model = "quintile")]
  returndd$Group = c("bottom20%","20-40%","40-60%","60-80%","top20%")
  returndd[, OR := exp(beta)]
  returndd[, OR_CI_L := exp(beta - 1.96*se)]
  returndd[, OR_CI_U := exp(beta + 1.96*se)]
  returndd[, Pvalue := p_value]
  returndd[, Ncontrols := c(nrow(combined[pheno == 0 & PRS_bin_label == 1]),
                            nrow(combined[pheno == 0 & PRS_bin_label == 2]),
                            nrow(combined[pheno == 0 & PRS_bin_label == 3]),
                            nrow(combined[pheno == 0 & PRS_bin_label == 4]),
                            nrow(combined[pheno == 0 & PRS_bin_label == 5])
  )]
  returndd[, Ncases := c(nrow(combined[pheno == 1 & PRS_bin_label == 1]),
                         nrow(combined[pheno == 1 & PRS_bin_label == 2]),
                         nrow(combined[pheno == 1 & PRS_bin_label == 3]),
                         nrow(combined[pheno == 1 & PRS_bin_label == 4]),
                         nrow(combined[pheno == 1 & PRS_bin_label == 5]) )]
  
  return(returndd[,-1:-3])
}

write.csv(assostats, "Quantile_analysis_binary_traits.csv", row.names = F)




#----- (2) Quantile analysis for quantitative traits -----
assostats = foreach(ii = 1:nrow(bestPGS), .combine = rbind) %do% {
  
  myisbinary = !is.na(bestPGS$AUC_full[ii])
  if(myisbinary) {return(NULL)}
  
  myPGSID = bestPGS$PGS_ID[ii]
  mytrait = bestPGS$phenotype[ii]
  
  
  # read the PGS
  myPGS = fread(paste0("yourpath/plink2_score_sum_",myPGSID,".txt"))
  # scale 
  myPGS$PRS_scaled = scale(myPGS$score_sum)
  
  # merge PGS to the phenotype table
  mydd = merge(phenotable_2ndexcl, myPGS[,.(ID, PRS_scaled)], by = "ID")
  
  # use age_at_event for quantitative trait
  cc = paste0("age_", mytrait)
  mydd[, age := get(cc)]
  
  # remove missing data
  mydd$pheno = mydd[, get(mytrait)]
  mydd = mydd[!is.na(pheno)]
  
  # IRN
  mydd$pheno = rankNorm(mydd$pheno)
  
  
  # Quintile analysis
  combined = mydd
  combined[, PRS_bin_label := cut(PRS_scaled, quantile(combined$PRS_scaled, probs = 0:5/5), labels = FALSE, include.lowest = TRUE)]
  combined$PRS_bin_label = factor(combined$PRS_bin_label, levels = c(3,1,2,4,5))
  
  # linear regression
  myfit1 = lm(Modelbin, data = combined)
  quincoef = summary(myfit1)$coefficients
  
  # OR, CI, and Pvalue
  metrics = as.data.table(quincoef[2:5, -3])
  names(metrics) = c("beta","se","Pvalue")
  returndd = rbind(metrics[1:2], data.table(beta = 0, se = 0, Pvalue = 1),
                   metrics[3:4])
  returndd[, `:=` (phenotype = mytrait, PGS_ID = myPGSID, 
                   ancestry = "ELGH", Model = "quintile")]
  returndd$Group = c("bottom20%","20-40%","40-60%","60-80%","top20%")
  returndd[, N := c(nrow(combined[PRS_bin_label == 1]),
                    nrow(combined[PRS_bin_label == 2]),
                    nrow(combined[PRS_bin_label == 3]),
                    nrow(combined[PRS_bin_label == 4]),
                    nrow(combined[PRS_bin_label == 5]) )]
  
  returndd = returndd[,.(phenotype, PGS_ID, ancestry, Model, Group, N, 
                         beta, se, Pvalue)]
  
  return(returndd)
}

write.csv(assostats, "Quantile_analysis_quantitative_traits.csv", row.names = F)



