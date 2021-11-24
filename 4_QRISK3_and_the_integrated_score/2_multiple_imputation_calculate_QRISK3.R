#------------------------------------------------------
# 2021-02-17
#
# Aims:
#   to impute missing data in QRSIK3 variables using 
#   the "mice" package (multiple imputation)
#
# Assessment date:
AssessmentDate = as.POSIXct("2010-01-01", format = "%Y-%m-%d")
#------------------------------------------------------
### Functions are defined in: https://github.com/QinqinHuang/GnH28k_polygenic_scores/blob/main/My_functions_GnH_PGS.R
source("My_functions_GnH_PGS.R")

library(QRISK3)
library(mice)
library(pheatmap)


#----- Load data -----
## QRISK3 variables with missing data
QRISK3 = readRDS("QRISK3_variables_prior_2010-01-01.RDS")

# Number of non-missing data per individual
QRISK3$Nnonmissing = sapply(1:nrow(QRISK3), function(x) {
  return(length(which(!is.na(c(QRISK3$Height[x],
                               QRISK3$Weight[x],
                               QRISK3$SBP[x],
                               QRISK3$TC[x],
                               QRISK3$HDL[x])))))
})

# keep samples with ≥3 non-missing: N=9477; 438 cases & 9039 controls
QRISK3_3nonmissing = QRISK3[Nnonmissing >= 3]


#----- Heatmap to show missingness (Figure S3) -----
# QRISK3 variables with missing data
varlist = c("Weight", "Height", "TC", "HDL", "SBP", "SBP2ysd", "smoking")

# NA: 1, non-NA: 0
QRISK3varNA = is.na(ddplot[, ..varlist])
QRISK3varNA = 1*QRISK3varNA
QRISK3varNA = as.data.frame(QRISK3varNA)
rownames(QRISK3varNA) = ddplot$pseudoNHSnumber

# heatmap for all individuals 
annotation_col_df = data.frame(Excluded = ddplot$flag)
rownames(annotation_col_df) = ddplot$pseudoNHSnumber
# clustering samples/columns
pdf("Heatmap_missingness_QRISK3_2021Jan__all_clusterCols.pdf", 
    width = 6, height = 4)
pheatmap(t(QRISK3varNA), 
         cluster_cols = T, cluster_rows = F,
         annotation_col = annotation_col_df, labels_col = rep("",nrow(annotation_col_df)),)
dev.off()



#----- Impute continuous values -----
QRISK3toimp = QRISK3_3nonmissing[,.(pseudoNHSnumber, ID, gender, UMAPancestry, age=age_QRISK3, Weight, Height, TC, HDL, SBP, SBP2ysd, smoking)]

# mice - using pmm (Predictive mean matching)
md.pattern(QRISK3toimp) # missing patterns
# m: number of imputed datasets
# maxit: no. of iterations
# pmm: predictive mean match
QRISK3imp = mice(QRISK3toimp, m=5, maxit = 50, method = 'pmm', seed = 500)
# this object has only imputed values for NAs

# get complete data, use the last one (5th)
QRISK3compl = as.data.table(complete(QRISK3imp, 5))

cc = setdiff(names(QRISK3_3nonmissing), names(QRISK3compl))
QRISK3compl = cbind(QRISK3compl, QRISK3_3nonmissing[, ..cc])


saveRDS(QRISK3compl, "QRISK3_inputdata_2010-01-01__nonmissing3_9477_impute_pmm.RDS")


#----- calcualte QRISK3 -----
QRISK3ELGH = QRISK3_2017(data = as.data.frame(QRISK3compl), patid = "pseudoNHSnumber", 
                         gender = "gender", age = "age",
                         atrial_fibrillation = "AF", 
                         atypical_antipsy = "antipsy",
                         regular_steroid_tablets = "steorid", 
                         erectile_disfunction = "ED",
                         migraine = "Migraine", rheumatoid_arthritis = "RA",
                         chronic_kidney_disease = "CKD345", severe_mental_illness = "SMI",
                         systemic_lupus_erythematosis = "SLE",
                         blood_pressure_treatment = "BPmed", diabetes1 = "T1D",
                         diabetes2 = "T2D", weight = "Weight", height = "Height",
                         ethiniciy = "ethnicity", heart_attack_relative = "FH_CHD",
                         cholesterol_HDL_ratio = "TC_HDL_ratio", 
                         systolic_blood_pressure = "SBP",
                         std_systolic_blood_pressure = "SBP2ysd", 
                         smoke = "smoking", townsend = "townsend")
fwrite(QRISK3ELGH, "QRISK3_calculated_2010-01-01__nonmissing3_9477_impute_pmm.txt", sep = "\t")






