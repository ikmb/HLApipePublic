#!/bin/bash 
##############################################
# Prepare data for downstream analysis
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################
prefix_orig=$1
prefix=$2
DUPLICATES_R=`which readDataDUPLICATES.R`
N=$3
## 
base_name=$prefix"_6_29_34"
base_name_phase=$prefix".phase"
fam_phase=$base_name_phase".fam"
plink --allow-no-sex --bfile $prefix_orig --make-bed --out $prefix


## PREPARE IMPUTATION                
plink --allow-no-sex --bfile $prefix  --chr 6 --from-mb 29 --to-mb 34 --make-bed --out $base_name

## PREPARE PHASING
Rscript $DUPLICATES_R $prefix
plink --allow-no-sex --bfile $prefix --exclude exclude_duplicates.txt --chr 6 --from-mb 25 --to-mb 35 --make-bed --out $base_name_phase
split -l $N  -d $fam_phase "split_files"
