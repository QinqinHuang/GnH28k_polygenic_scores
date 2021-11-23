#--------------------------------------------------
# 2020-05-05
#
# Aim: to match the GWAS summary data with validation
# cohort and using the SNP IDs in the latter
#
#
# "gwas": GWAS summary data, with updated column
# names:
#  SNP, CHR, BP, A1 (effect allele), A2, Beta, SE, Pval
#
# "impsite": SNP site file for the validation cohort 
# should contain:
#  SNP, CHR, BP, REF, ALT, MAF
#--------------------------------------------------
library(data.table)
library(foreach)
library(ggplot2)


#----- My functions -----
# To merge GWAS data and imputed genotype data by CHR, POS, and alleles
# Exclude multi-allelic sites
# Update the variant IDs in the GWAS data
# Remove variants with non-matching alleles

matchedgwas = function(gwas, impsite) {
  
  # Remove multi-allelic SNPs from GWAS dataset
  gwas[,SNP:=paste0(CHR,":",BP)]
  gwas = gwas[!SNP %in% gwas[duplicated(SNP)]$SNP]
  print(paste0(" Number of GWAS variants excluding multi-allelic variants: ", nrow(gwas)))
  gwas$SNP = NULL
  
  # Merge two datasets
  gwas_overlapping = merge(gwas, impsite, by = c("CHR","BP"))
  print(paste0(" Number of variants with matched position: ", nrow(gwas_overlapping)))
  
  # In the GWAS dataset, A1 the effect allele is usually the REF, and A2 is usually the ALT
  # Keep variants if alleles are matching the imputed data
  gwas_overlapping[(A1==REF & A2==ALT), flag:="same"]
  gwas_overlapping[(A1==ALT & A2==REF), flag:="switch"]
  if(sum(is.na(gwas_overlapping$flag))>0) {
    gwas_overlapping[A1=="A", A1_flipped:="T"]
    gwas_overlapping[A1=="T", A1_flipped:="A"]
    gwas_overlapping[A1=="C", A1_flipped:="G"]
    gwas_overlapping[A1=="G", A1_flipped:="C"]
    gwas_overlapping[A2=="A", A2_flipped:="T"]
    gwas_overlapping[A2=="T", A2_flipped:="A"]
    gwas_overlapping[A2=="C", A2_flipped:="G"]
    gwas_overlapping[A2=="G", A2_flipped:="C"]
    gwas_overlapping[(is.na(flag) & A1_flipped==REF & A2_flipped==ALT), flag:="flipped"]
    gwas_overlapping[(is.na(flag) & A1_flipped==ALT & A2_flipped==REF), flag:="flipped_switch"]
  }
  
  print(paste0(" Number of variants with non-matching alleles: ", length(which(is.na(gwas_overlapping[,flag])))))
  gwas_overlapping = gwas_overlapping[!is.na(flag),]
  print(paste0(" Number of overlapping variants: ", nrow(gwas_overlapping)))
  print(paste0("    Among which ", sum(gwas_overlapping$flag == "same"), " same, ", sum(gwas_overlapping$flag == "switch"), " need to update A1/A2 alleles"))
  if("flipped" %in% gwas_overlapping$flag | "flipped_switch" %in% gwas_overlapping$flag) {
    print(paste0("    Also, ", sum(gwas_overlapping$flag == "flipped"), " flipped, ", sum(gwas_overlapping$flag == "flipped_switch"), " flipped and switched"))
  }
  
  # Update beta when A1/A2 is not REF/ALT alleles
  # If change the sign of beta i should use REF instead of A1 as the effect allele
  gwas_overlapping[flag=="switch" | flag=="flipped_switch", `:=` (Beta = -Beta)]
  
  print(paste0(" Number of variants with MAF < 1%: ", sum(gwas_overlapping$MAF <0.01)))
  
  return(gwas_overlapping[,.(SNP,CHR,BP,A1=REF,A2=ALT,Beta,SE,Pval,MAF)])
}


