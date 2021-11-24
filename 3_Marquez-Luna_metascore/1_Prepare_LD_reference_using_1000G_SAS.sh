#------------------------------------------------------------
# Last modified 2021-03-15
#
# Aim: 
# (1) to extract SAS 489 samples from 1000G 
# (2) to convert to plink files: geno 0.05 MAF 0.1%
# (3) to update SNP IDs to CHR:BP
#------------------------------------------------------------
### working directory


# make a list of 489 SAS samples
cat ../integrated_call_samples_v3.20130502.ALL.panel |grep SAS |cut -f 1 > SAS_list_489.txt
awk '{print "1 "$1}' SAS_list_489.txt > SAS_list_489_constfid1.txt


# 22 chromosomes
for cc in {1..22}
do

# convert vcf to plink files: keep ATCG SNPs only
plink --vcf ../ALL.chr"$cc".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
--keep-allele-order --snps-only 'just-acgt' --const-fid 1 --memory 4000 \
--make-bed -out ALL.chr"$cc"_1000G_SNPs

# extract SAS samples
plink --bfile ALL.chr"$cc"_1000G_SNPs --keep SAS_list_489_constfid1.txt \
--keep-allele-order --make-bed -out ALL.chr"$cc"_1000G_SNPs_SAS489 --memory 4000

# geno 0.05, MAF 0.1% (at least one MAC)
plink --bfile ALL.chr"$cc"_1000G_SNPs_SAS489 --geno 0.05 --maf 0.001 --memory 4000 \
--keep-allele-order --make-bed -out ALL.chr"$cc"_1000G_SNPs_SAS489_geno0.05_MAF0.001

done



# Update SNP ID to chr:pos to replace "." 
for ii in {1..22} 
do
mv ALL.chr"$ii"_1000G_SNPs_SAS489_geno0.05_MAF0.001.bim ALL.chr"$ii"_1000G_SNPs_SAS489_geno0.05_MAF0.001.rsID.bim
cat ALL.chr"$ii"_1000G_SNPs_SAS489_geno0.05_MAF0.001.rsID.bim | \
awk '{print$1 "\t" $1 ":" $4 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > \
ALL.chr"$ii"_1000G_SNPs_SAS489_geno0.05_MAF0.001.bim
done


# try merging
plink --bfile ALL.chr1_1000G_SNPs_SAS489_geno0.05_MAF0.001 --merge-list plink_to_merge_21.txt \
--keep-allele-order --make-bed -out plink_1000G_SAS489_SNPs_geno0.05_MAF0.001 \
--memory 10000
# warning: 12 SNPs with 3+ alleles present.

# excluding them
mv plink_1000G_SAS489_SNPs_geno0.05_MAF0.001-merge.missnp SNPs_12_with_morethan2_alleles.txt
for ii in {1..22}
do
plink --bfile ALL.chr${ii}_1000G_SNPs_SAS489_geno0.05_MAF0.001 \
--keep-allele-order --exclude SNPs_12_with_morethan2_alleles.txt \
--make-bed -out ALL.chr${ii}_1000G_SNPs_SAS489_geno0.05_MAF0.001_biallelic --memory 4000
done


# merge 22 chr
for ii in {2..22}; do echo ALL.chr"$ii"_1000G_SNPs_SAS489_geno0.05_MAF0.001_biallelic; done > plink_to_merge_21_biallelic.txt
plink --bfile ALL.chr1_1000G_SNPs_SAS489_geno0.05_MAF0.001_biallelic --keep-allele-order \
--merge-list plink_to_merge_21_biallelic.txt \
--make-bed -out plink_1000G_SAS489_SNPs_geno0.05_MAF0.001_biallelic --memory 10000


# hwe 1e-6
plink --bfile plink_1000G_SAS489_SNPs_geno0.05_MAF0.001_biallelic --keep-allele-order \
--hwe 1e-6 --make-bed -out plink_1000G_SAS489_SNPs_geno0.05_MAF0.001_biallelic_hwe6 \
--memory 2000







	
	



