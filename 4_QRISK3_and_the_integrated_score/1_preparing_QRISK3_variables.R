#------------------------------------------------------
# 2021-02-17
# 
# Aim: 
# to prepare data for QRISK3 calculation
# 
# Assessment date:
AssessmentDate = as.POSIXct("2010-01-01", format = "%Y-%m-%d")
# 
# *Medication: ≥2 prescriptions, and the latest within 1 month
# *conditions: 1 if the earliest date is before assessment date
# *measurements: the most recent one before assessment date
#------------------------------------------------------
### Functions are defined in: https://github.com/QinqinHuang/GnH28k_polygenic_scores/blob/main/My_functions_GnH_PGS.R
source("My_functions_GnH_PGS.R")
library(QRISK3)


#----- Load preparation data -----
## Weight, Height, TC, HDL, SBP, smoke data; Medications; Conditions
weight = readRDS("Weight_multiple.RDS")
height = readRDS("Height_multiple.RDS")
TC = readRDS("TC_multiple.RDS")
HDL = readRDS("HDL_multiple.RDS")
SBP = readRDS("SBP_clean_multiple_measurements_18yo_range_60_254.RDS")
smoking = readRDS("Smoking_status_2019primary.RDS")

# medications
steorid = readRDS("Steroid_medications.RDS")
antipsy = readRDS("Antipsychotics_medications.RDS")
BPmed = readRDS("BP_medications.RDS")

# conditions considered in the QRISK3 algorithm
conditions = readRDS("Relevant_conditions_QIRSK3_earliest_Dx.RDS")


#----- CAD cases and controls that are considered in QRISK3 -----
CADpheno = phenotable[,.(pseudoNHSnumber, ID, age_at_recrt, gender, CAD = CAD_Khera)]
CADpheno = CADpheno[!is.na(CAD)]

# earliest code
CADdateE = readRDS(".CAD_cases_Khera_earliest_code_date.RDS")

# prevalent cases and incident cases
CADdateE[, casetype := ifelse(CAD_Dx_DateE <= AssessmentDate, yes = "prevalent", no = "incident")]
CADpheno = merge(CADpheno, CADdateE, by = "pseudoNHSnumber", all.x = T)
CADpheno[is.na(casetype), casetype := "control"]

# year of birth and ancestry
cov = fread("SampleIDs_AgeAtRecrt_yob_corrected_gender_ancestry.txt", colClasses = "character")
CADpheno = merge(CADpheno, cov[,.(ID, yob, UMAPancestry)], by = "ID")
CADpheno[, yob := as.numeric(yob)]

# Age at CAD diagnosis
CADpheno[, age_CAD_Dx := year(CAD_Dx_DateE) - yob]

# Age at QRISK3 assessment
CADpheno[, age_QRISK3 := 2009 - yob]

# QRISK3 age range: 25-84
CADpheno[, flag := "QRISK3"]
CADpheno[age_QRISK3 > 84 | age_QRISK3 < 25, flag := "excluded_outside_age_range"]
CADpheno[casetype == "prevalent", flag := "excluded_prevalent_cases"]

CADpheno[, QRISK3_Date := AssessmentDate]

saveRDS(CADpheno, "CAD_case_control_list_QRISK3.RDS")



#----- QRISK3 variables -----
CADpheno = CADpheno[flag == "QRISK3"]

# prepare the matrix for the QRISK3_2017 function 
forQRISK3 = CADpheno

#-----
# functions
# (1) getcondition - diagnosis data
# (2) getmedication - medication data
# (3) latestmeasurement - measurement data

getcondition = function(mydata, whichcondition) {
  # relevant conditions
  myconditions = conditions[Condition == whichcondition]
  tempdd = merge(mydata[,.(pseudoNHSnumber, QRISK3_Date)],
                 myconditions[,.(pseudoNHSnumber, Date)], by = "pseudoNHSnumber")
  
  # if the diagnosis date is earlier than QRISK3 assessment date
  tempdd[, flag := Date <= QRISK3_Date]
  #print(table(tempdd$flag))
  
  mydata$newcol = 0
  mydata[pseudoNHSnumber %in% tempdd[flag==T]$pseudoNHSnumber, newcol := 1]
  #print(table(mydata$newcol))
  names(mydata)[ncol(mydata)] = whichcondition
  
  return(mydata)
}

getmedication = function(mydata, whichmed, mymed) {
  # relevant medication
  tempdd = merge(mydata[,.(pseudoNHSnumber, QRISK3_Date)],
                 mymed[,.(pseudoNHSnumber, dateE, dateL)], by = "pseudoNHSnumber")
  
  # remove data that are after QRISK3 assessment date
  cat("  Remove", nrow(tempdd[QRISK3_Date < dateE]), "rows - prescription after QRISK3 assessment date\n")
  tempdd = tempdd[QRISK3_Date >= dateE]
  
  # 1. remove individuals with only one prescription
  singlepre = tempdd[!pseudoNHSnumber %in% tempdd[duplicated(pseudoNHSnumber)]$pseudoNHSnumber]
  singlepre = singlepre[dateL==dateE]
  cat("  Remove", nrow(singlepre), "individuals with only one prescription before QRISK3 assessment date\n")
  tempdd = tempdd[!pseudoNHSnumber %in% singlepre$pseudoNHSnumber]
  
  # 2. if the QRSIK3 date is within 28 days after the latest prescription date 
  tempdd[, flag := QRISK3_Date <= dateL + 28]
  #print(table(tempdd$flag))
  cat("  unique N =", length(unique(tempdd[flag == T]$pseudoNHSnumber)), "individuals\n")

  mydata$newcol = 0
  mydata[pseudoNHSnumber %in% tempdd[flag==T]$pseudoNHSnumber, newcol := 1]
  #print(table(mydata$newcol))
  names(mydata)[ncol(mydata)] = whichmed
  
  return(mydata)
}

latestmeasurement = function(mydata, whichmeas, mymeas) {
  # relevant data
  tempdd = merge(mydata[,.(pseudoNHSnumber, QRISK3_Date)],
                 mymeas[,.(pseudoNHSnumber, date, value)], by = "pseudoNHSnumber")
  
  # remove data that are after QRISK3 assessment
  cat("  Remove", nrow(tempdd[QRISK3_Date < date]), "rows - data after QRISK3 assessment date\n")
  tempdd = tempdd[QRISK3_Date >= date]
  
  # most recent value
  tempdd = tempdd[order(date, decreasing = T)]
  tempdd = tempdd[!duplicated(pseudoNHSnumber)]
  
  mydata = merge(mydata, tempdd[,.(pseudoNHSnumber, value, date)], 
                 by = "pseudoNHSnumber", all.x = T)
  #print(table(!is.na(mydata$value)))
  names(mydata)[(ncol(mydata)-1):ncol(mydata)] = paste0(whichmeas, c("","_date"))
  
  return(mydata)
}
#-----


## gender - 1: women & 0: men (ELGH data: male 1, female 2)
forQRISK3[, gender := gender - 1]
table(forQRISK3[CAD == 1]$gender)

## age - Specify the age of the patient in year (e.g. 64 years-old)

## atrial_fibrillation - (0: No, 1:Yes)
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "AF")

## atypical_antipsy - On atypical antipsychotic medication? (0: No, 1:Yes)
forQRISK3 = getmedication(mydata = forQRISK3, whichmed = "antipsy", mymed = antipsy)

## regular_steroid_tablets - On regular steroid tablets? (0: No, 1:Yes)
forQRISK3 = getmedication(mydata = forQRISK3, whichmed = "steorid", mymed = steorid)

## erectile_disfunction - A diagnosis of or treatment for erectile disfunction? (0: No, 1:Yes)
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "ED")

## migraine - Do patients have migraines? (0: No, 1:Yes)
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "Migraine")

## rheumatoid_arthritis - (0: No, 1:Yes)
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "RA")

## chronic_kidney_disease - Chronic kidney disease (stage 3, 4 or 5)? (0: No, 1:Yes)
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "CKD345")

## severe_mental_illness - 0: No, 1:Yes
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "SMI")

## systemic_lupus_erythematosis - 0: No, 1:Yes
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "SLE")

## blood_pressure_treatment - On blood pressure treatment? (0: No, 1:Yes)
forQRISK3 = getmedication(mydata = forQRISK3, whichmed = "BPmed", mymed = BPmed)

## diabetes1
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "T1D")

## diabetes2
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "T2D")

## weight - Weight of patients (kg)
forQRISK3 = latestmeasurement(mydata = forQRISK3, whichmeas = "Weight", mymeas = weight)

## height - Height of patients (cm)
forQRISK3 = latestmeasurement(mydata = forQRISK3, whichmeas = "Height", mymeas = height)

## ethniciy - 1 White or not stated 2 Indian 3 Pakistani 4 Bangladeshi 5 Other Asian
# 6 Black Caribbean 7 Black African 8 Chinese 9 Other ethnic group
forQRISK3[, ethnicity := 5]
forQRISK3[UMAPancestry == "Pakistani", ethnicity := 3]
forQRISK3[UMAPancestry == "Bangladeshi", ethnicity := 4]

## heart_attack_relative - Angina or heart attack in a 1st degree relative < 60? (0: No, 1:Yes)
forQRISK3 = getcondition(mydata = forQRISK3, whichcondition = "FH_CHD")


## cholesterol_HDL_ratio - TC/HDL ratio (range from 1 to 11, e.g. 4)
# latest TC and HDL value
forQRISK3 = latestmeasurement(mydata = forQRISK3, whichmeas = "TC", mymeas = TC)
#*** No HDL data available before 2010
forQRISK3 = latestmeasurement(mydata = forQRISK3, whichmeas = "HDL", mymeas = HDL)
# impute using the female and male averages calculated using later measurements
forQRISK3[gender == 1, TC_HDL_ratio := 3.905]
forQRISK3[gender == 0, TC_HDL_ratio := 4.882]


## systolic_blood_pressure - Systolic blood pressure (mmHg, e.g. 180 mmHg)
forQRISK3 = latestmeasurement(mydata = forQRISK3, whichmeas = "SBP", mymeas = SBP)

## std_systolic_blood_pressure - Standard deviation of at least two most recent systolic blood pressure readings (mmHg)
# recent SBP values before QRISK3 assessment
recentSBP = merge(forQRISK3[,.(pseudoNHSnumber, QRISK3_Date)],
                  SBP[,.(pseudoNHSnumber, date, value)], by = "pseudoNHSnumber")
recentSBP = recentSBP[QRISK3_Date >= date]
recentSBP = recentSBP[order(date, decreasing = T)]
# difference in time compared with the most recent measurement
mostrecentSBP = recentSBP[!duplicated(pseudoNHSnumber)]
mostrecentSBP = mostrecentSBP[, .(pseudoNHSnumber, date_recent = date)]
recentSBP = merge(recentSBP, mostrecentSBP, by = "pseudoNHSnumber")
recentSBP[, timediff := date_recent - date]
# define most recent - within 2 years
recentSBP_2y = recentSBP[timediff <= 2*365*24*3600]
indcount = recentSBP_2y[, .N, by = pseudoNHSnumber]
recentSBP_2y_sd = recentSBP_2y[, sd(value), by = pseudoNHSnumber]
names(recentSBP_2y_sd)[2] = "SBP2ysd"

forQRISK3 = merge(forQRISK3, recentSBP_2y_sd, by = "pseudoNHSnumber", all.x = T)


## smoke - 1 non-smoker 2 ex-smoker 3 light smoker (less than 10) 4 moderate smoker (10 to 19) 5 heavy smoker (20 or over)
recentsmoking = merge(forQRISK3[, .(pseudoNHSnumber, QRISK3_Date)], 
                      smoking[, .(pseudoNHSnumber, DateL, smoking)],
                      by = "pseudoNHSnumber")
# Nonsmoker -> 1
# if the QRISK3 assessment date is before the CurrentSmoker_LatestDate -> 3, 4, or 5
# if CurrentSmoker_LatestDate is NA, or the QRISK3 assessment date is later -> 2
recentsmoking[QRISK3_Date > DateL, smoking := 2]
# treat them as ex-smoker
recentsmoking[is.na(DateL), smoking := 2]
table(recentsmoking$smoking)

forQRISK3 = merge(forQRISK3, recentsmoking[,.(pseudoNHSnumber, smoking)], 
                  by = "pseudoNHSnumber", all.x = T)

## townsend - Townsend deprivation scores
forQRISK3$townsend = 3.307

saveRDS(forQRISK3, "QRISK3_variables_prior_2010-01-01.RDS")




