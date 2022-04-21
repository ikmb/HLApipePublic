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

CERTAINTY_R=`which certainty.R`

######################################
# Get HLA types and prepare validation
###################################### 

# EXCLUDE ATCG
awk 'length($5) + length($6) ==2' $name.bim | awk '{ if($5=="A" && $6=="T" || $5=="T" && $6=="A" || $5=="G" && $6=="C" || $5=="C" && $6=="G"){}else{print}}' > $name.atcg
 
plink --bfile $name --extract $name.atcg  --make-bed --out  $name

cp $name.bim $name.old.bim

awk '$2="chr"$1":"$4' $name.old.bim > $name.bim 

#SHAPEIT=shapeit

$SHAPEIT \
	-check \
	-B $name  \
	-M $REFERENCE_DIR_IMPUTE2/genetic_map_chr6"_combined_b37.txt" \
	--input-ref $REFERENCE_DIR_IMPUTE2/1000GP_Phase3_chr6.hap.gz $REFERENCE_DIR_IMPUTE2/1000GP_Phase3_chr6.legend.gz $REFERENCE_DIR_IMPUTE2/1000GP_Phase3.sample \
	--output-log $name.check \
	--thread 1

cat $name.*strand | grep Strand | awk '{print $4}' |sort | uniq > $name.flip

plink --bfile $name --flip  $name.flip --make-bed --out $name



$SHAPEIT \
	-check \
	-B $name  \
	-M $REFERENCE_DIR_IMPUTE2/genetic_map_chr6"_combined_b37.txt" \
	--input-ref $REFERENCE_DIR_IMPUTE2/1000GP_Phase3_chr6.hap.gz $REFERENCE_DIR_IMPUTE2/1000GP_Phase3_chr6.legend.gz $REFERENCE_DIR_IMPUTE2/1000GP_Phase3.sample \
	--output-log $name.check \
	--thread 1



if [[ $(wc -l $name.fam | awk '{print $1}') -gt 100 ]]
then

cat $name*.strand | grep Strand | awk '{print $4}' |sort| uniq> $name.exclude
plink --bfile $name --exclude  $name.exclude --make-bed --out $name

$SHAPEIT \
--input-bed $name \
--input-map $REFERENCE_DIR_IMPUTE2/genetic_map_chr6"_combined_b37.txt" \
--input-thr 0.9 \
--missing-code 0 \
--states 100 \
--window 2 \
--effective-size 18000 \
--burn 7 \
--prune 8 \
--main 20 \
--seed 123456789 \
--output-graph $name.hgraph \
--thread $CPUS \
--output-log $name.hgraph


else

cat $name*.strand | awk '{print $4}' |sort| uniq> $name.exclude
plink --bfile $name --exclude  $name.exclude --make-bed --out $name


$SHAPEIT \
--input-ref  $REFERENCE_DIR_IMPUTE2/1000GP_Phase3_chr6.hap.gz  $REFERENCE_DIR_IMPUTE2/1000GP_Phase3_chr6.legend.gz  $REFERENCE_DIR_IMPUTE2/1000GP_Phase3.sample \
--input-bed $name \
--input-map $REFERENCE_DIR_IMPUTE2/genetic_map_chr6"_combined_b37.txt" \
--input-thr 0.9 \
--missing-code 0 \
--force \
--states 100 \
--window 2 \
--effective-size 18000 \
--burn 7 \
--prune 8 \
--main 20 \
--seed 123456789 \
--output-graph $name.hgraph \
--thread $CPUS \
--output-log $name.hgraph


fi
#######################
# Phasing certainty
#######################

#SHAPEIT=shapeit

for r in $(seq 1 100); do
    $SHAPEIT -convert \
            --input-graph $name.hgraph \
            --output-sample $name.S$r \
   	    --output-log $name.S$r.log \
            --seed $r;
done

Rscript $CERTAINTY_R $name $(ls *haps)


sed -i 's/chr6://g' $name.certainty.exclude  
sed -i 's/"//g' $name.certainty.exclude  

$SHAPEIT -convert \
--exclude-snp $name.certainty.exclude \
--input-graph $name.hgraph  \
--output-max $name \
--output-log $name.convert



