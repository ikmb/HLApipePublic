#!/bin/bash

##############################################
# Impute missing data from the dataset using
# dataset itself as reference to fill gaps
# author:  Mareike Wendorff, Frauke Degenhardt
# contact: f.degenhardt@ikmb.uni-kiel.de
# March 2020
##############################################

bedbimfam_ref=$1
bedbimfam_beagle=$2
beagle_map=$3
beagle_jar=$4

fam=${bedbimfam_ref}.fam

### run beagle if defined

# code file to .vcf
plink --bfile  $bedbimfam_ref --list-duplicate-vars ids-only suppress-first
cat  $bedbimfam_ref.bim | grep [DI] | cut -f2 >> plink.dupvar
plink --bfile $bedbimfam_ref --exclude plink.dupvar --recode vcf --allow-no-sex --out tmp



# run beagle
java -jar $beagle_jar impute=true gt=tmp.vcf map=$beagle_map out=$bedbimfam_beagle'.tmpimp'

## unzip file
	gunzip $bedbimfam_beagle'.tmpimp.vcf'

## vcf back to bim/bed/fam
	plink --vcf ${bedbimfam_beagle}.tmpimp.vcf --const-fid --make-bed --allow-no-sex --out ${bedbimfam_beagle}

## replace new fam file by old one, because of problems with _
	cp ${fam} ${bedbimfam_beagle}.fam
##fi

