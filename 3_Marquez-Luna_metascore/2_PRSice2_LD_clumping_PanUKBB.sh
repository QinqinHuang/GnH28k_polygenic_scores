#------------------------------------------------------------
# 2021-03-15
# 
# Aims:
# Run PRSice2 to do LD clumping in 1kG SAS, and calculate
#  PRSs in the validation cohort
#
# LD r2: 0.1, 0.2, 0.5, 0.8
# P-values: 1, 0.8, 0.5, 0.4, 0.3, 0.2, 0.1, 
#			0.08, 0.05, 0.02, 0.01, 0.005, 1e-3, 5e-4, 1e-4,
#			5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8
#
#------------------------------------------------------------
### update your own paths
### "GWASsumm": the path to the GWAS summary data with matched variants and updated variant IDs
### "genodir": the path to the folder where bgen data of the validation cohort are located
### "LDfile": the path to the 1000 Genomes EUR reference plink data
### "LDR2": LD r2 cutoff


PRSICE="PRSice_2.2.11"

/software/hgi/installs/anaconda3/envs/hgi_base/bin/Rscript \
~/software/PRSice2/"$PRSICE"/PRSice.R \
--dir ~/R/x86_64-pc-linux-gnu-library/3.6 \
--prsice ~/software/PRSice2/"$PRSICE"/PRSice_linux \
--base "$GWASsumm" \
--beta \
--A1 A1 --A2 A2 --stat Beta --bp BP --chr CHR --pvalue Pval --snp SNP \
--target "$genodir"/chr#_R2_0.3_MAF_0.001_withCHR,"$genodir"/Feb28k_R2_0.3_MAF_0.001.samples \
--type bgen \
--fastscore \
--bar-levels 1,0.8,0.5,0.4,0.3,0.2,0.1,0.08,0.05,0.02,0.01,0.005,1e-3,5e-4,1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7,5e-8 \
--all-score \
--print-snp \
--ld "$LDdir"/1000G_SAS_filtered \
--allow-inter \
--clump-r2 "$LDR2" \
--out PRSice_output_"$LDR2" \
--no-regress \
--memory 4000
