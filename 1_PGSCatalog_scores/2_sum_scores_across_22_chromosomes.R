#------------------------------------------------------------
# 2021-01-10
# 
# Sum scores across 22 autosomes.
#------------------------------------------------------------
library(data.table)
library(foreach)


## The list of scores from the PGS Catalog. 
PGSmeta = fread("Benchmark_PGS_metadata_updated_NSNPs.csv")

for(ii in 1:nrow(PGSmeta)) {
  PGSID = PGSmeta$PGS_ID[ii]
  
  # ELGH and eMerge
  for(cohort in c("ELGH", "eMerge")) {
    
    # the 22 files for the score
    filelist = grep(paste0("plink2_score_",PGSID,"_chr"), dir(cohort), value = T)
    filelist = grep(".sscore", filelist, value = T)
    
    scores = foreach(ff = filelist, .combine = rbind) %do% {
      score = fread(paste0(cohort, "/", ff))
      names(score)[1] = "ID"
      return(score)
    }
    
    # check sample size
    if(length(filelist) != 22) {
      print(PGSmeta[ii])
      cat("  number of chr:", length(filelist), "\n")
    }
    if(nrow(scores)/length(filelist) != ifelse(cohort == "ELGH", yes = 28022, no = 43877)) {
      cat("not 28022 individuals\n")
      stop(PGSID)
    }
    # check number of SNPs
    if(sum(scores[ID == scores$ID[1]]$ALLELE_CT)/2 != PGSmeta$NSNPs_overlapping[ii]) {
      cat("error with Number of SNPs\n")
      stop(PGSID)
    }
    
    sumscores = scores[, sum(SCORE1_SUM), by = ID]
    names(sumscores)[2] = "score_sum"
    
    fwrite(sumscores, paste0(cohort,"_PGS/plink2_score_sum_",PGSID,".txt"), sep = "\t")
  } 
  
}



