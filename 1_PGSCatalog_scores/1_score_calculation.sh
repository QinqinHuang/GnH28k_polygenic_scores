#------------------------------------------------------------
# 2021-01-10
# 
# Aim:
#  to calculate polygenic scores using pgen files for 22 chrs
#  for ELGH GAsP and eMerge HRC
#------------------------------------------------------------

## $PGSID is the score ID in the PGS Catalog
## $scorefile is the harmonised score file with weights
## $pgenfile is the prefix of the plink2 pgen files

for ii in {1..22}
do
plink2 --pfile ${pgenfile} --memory 4000 \
--score ${scorefile} 1 3 4 'cols=nallele,dosagesum,scoresums' \
--out plink2_score_${PGSID}_chr${ii}
done








