#-------------------------------------------------------
# 2021-04-30
# 
# Aims:
#  Run association tests for PGS from the PGS Catalog
#-------------------------------------------------------
### Functions are defined in: https://github.com/QinqinHuang/GnH28k_polygenic_scores/blob/main/My_functions_GnH_PGS.R
source("My_functions_GnH_PGS.R")


#----- set up -----
### update with your path file

# "PGS_dir": the dir for calculated PGS

# "PGStrait_pairs": a file with all PGS-trait pairs that need to be tested


#----- test for associations -----
PGSpfm = foreach(ii = 1:nrow(PGStrait_pairs), .combine = rbind) %do% {
  
  myPGSID = PGStrait_pairs$PGS_ID[ii]
  mytrait = PGStrait_pairs$trait[ii]
  mytraittype = PGStrait_pairs$traittype[ii]
  outputfilename = paste0("stats_",mytrait,"_",myPGSID,".RDS")
  
  
  cat("----- working on", mytrait, myPGSID, ";", ii, "/", nrow(PGStrait_pairs) ,"-----\n")
  
  
  # read the PGS
  myPGS = fread(paste0(PGS_dir,"/plink2_score_sum_",myPGSID,".txt"))
  # scale 
  myPGS$PRS_scaled = scale(myPGS$score_sum)
  
  # merge PGS to the phenotype table
  mydd = merge(phenotable_2ndexcl, myPGS[,.(ID, PRS_scaled)], by = "ID")
  
  if(mytraittype == "quantitative") {
    cc = paste0("age_", mytrait)
    mydd[, age := get(cc)]
  }
  
  # trait
  mydd$pheno = mydd[, get(mytrait)]
  # remove missing data
  mydd = mydd[!is.na(pheno)]
  
  
  # test the PGS
  combined_stats = getstatistics(mydd, binary = ifelse(mytraittype == "binary", yes = T, no = F), 
                                 quandata = ifelse(mytraittype == "quantitative", yes = "IRN", no = "raw"),
                                 scale_PRS = F,
                                 bootstrapping = "yes_for_siginicant",
                                 fullmodel = fullmodel, nullmodel = nullmodel) 
  
  returndd = PGStrait_pairs[ii]
  returndd = cbind(returndd, combined_stats)
  
  if(!"bootstrapest" %in% names(returndd)) {
    returndd$bootstrapest = NA
  }
  
  saveRDS(returndd, outputfilename)
  return(returndd)
}






