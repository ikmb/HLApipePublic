#!/bin/bash 
##############################################
# Phase SNPs using SHAPEIT
# author:  Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################

name=$1
SHAPEIT=$2
REFERENCE_DIR_IMPUTE2=$3
CPUS=$4

CERTAINTY_R=`which phaseSNPsCERTAINTY.R`

######################################
# Get HLA types and prepare validation
###################################### 

# EXCLUDE ATCG
awk 'length($5) + length($6) ==2' $name.bim | awk '{ if($5=="A" && $6=="T" || $5=="T" && $6=="A" || $5=="G" && $6=="C" || $5=="C" && $6=="G"){}else{print}}' > $name.atcg
 
plink --bfile $name --extract $name.atcg  --make-bed --out  $name

cp $name.bim $name.old.bim

awk '$2="chr"$1":"$4' $name.old.bim > $name.bim 

#SHAPEIT=shapeit
plink --bfile $name --recode vcf --out $name 
$SHAPEIT --input $name.vcf --map  chr20.b37.gmap.gz --region 6 --mcmc-iterations 5b,1p,1b,1p,1b,1p,5m --window 2 --pbwt-mdr 0.2 --output phased.vcf.gz 


