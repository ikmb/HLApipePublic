![](images/ikmb_bfx_logo.png)

# HLApipePublic


This pipeline offers a workflow for HLA imputation using HIBAG (Zheng et al., 2014) and phasing for small datasets (Degenhardt et al., 2020) including utility scripts to evaluate the accuracy of the imputation.

Preprocessing:
- Alignment of a study to the imputation reference

Main:
- HLA, Amino acid & SNP imputation
- Phasing of HLA alleles

Utility: 
- Calculation of marginal probabilitites per allele of imputation results
- Calculation of alleles that have similar SNP haplotypes given the positions in your input data
- Calculation of alleles that are difficult to phase given your input data

Additional code includes: 
- Preparation of SNP/AA database from IMGT


## Dependencies 

- Nextflow
- SHAPEIT2 [https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#download]
- IMPUTE2 reference. [https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#reference]
- BEAGLE (optional) [https://faculty.washington.edu/browning/beagle/b4_1.html]
- Beagle map plink.chr6.GRCh37.map [http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/]
- PLINK 1.9 [https://www.cog-genomics.org/plink/1.9/]
- R (version 3.5+)
- R packages: parallel, ggplot2, data.table, reshape, reshape2, HIBAG, SNPRelate, grid, gridExtra, knitR, R.utils

## Installation

Change the following lines in the files to fit your requirements: 

nextflow.config:

```
params.impute2_reference_dir = /path/to/impute2/reference/files
params.ref_1000G =  /path/to/1000G/.bim/.bed/.fam annotation
params.shapeit = /path/to/shapeit2/executable
```

```
params {
  // Defaults only, expecting to be overwritten
  max_memory = 120.GB
  max_cpus = 8
  max_time = 36.h
  maxMultiqcEmailFileSize = 25.MB 
}
```

- R packages can be installed using anaconda and `conda env create -f environment.yml`
- Download the plink.chr6.GRCh37.map and place it into assets/beagle_map/


# Usage

## Basic execution
```
Usage:  nextflow run HLApipePublic --prefix FILE --reference_name REFERENCE --run_name NAME --shapeit SHAPEIT --impute2_reference_dir IMPUTE2_REF_DIR 
```
## Parameters

```--prefix```		An input prefix referencing a set of PLINK files \

```--reference_name``` 	Name of the reference imputation panel (see below for details) 

## Optional parameters:

General:

```--loci```     	Loci that should be imputed. Default: As specified in conf/resources.config.\


Software/References: 

```--shapeit```		Path to the SHAPEIT2 executable. \
```--impute2_ref_dir``` Path to the IMPUTE2 reference. \
```--beagle```        Location of the .jar file of BEAGLE4.1. \
```--do_beagle```    	This flag is optional and enables phasing using Beagle. Default: false. \
```--ref_1000G ```    Path to population used for PCA. PLINK files (hg19). Default: 1000G Phase 3 population. \
```--sample```	Path to sample file used for PCA. Default: 1000G Phase 3 population. \
```--subpop```    	Name of a sub population to use. Valid options are: AA, AFR, AMR, CHN, EAS, EUR, GER, IND, IRN, JPN, KOR, MLT. Can be used together with the IKMB reference.\

Others: 

```--email```         Email address to send reports to (enclosed in '') \
```--outdir```        Path to output directory. Default: results.


## Input
 - PLINK file in .bed/.bim/.fam format. Assembly must be the same as the assembly of the reference. Default here hg19.

## Outputs
<details> <summary>
alignToReference: .*.summary.refchecked.txt</summary>

| Name  | Description |
| ------------- | ------------- |
| chr | chromosome in input file|
| id | SNP id in input file  |
| pos | position (bp) in input file  |
| A1| minor allele in input file  |
| A2 | major allele in input file  |
| MAF | minor allele frequency in input file |
| ref.pos | position (bp) in reference model  |
| ref.A1 |  minor allele in reference model|
| ref.A2 | major allele in reference model  |
| ref.MAF | minor allele frequency in reference model  |
| action | action taken |
| type | type of genotype AT/CG or not |
</details>
<details> <summary>
imputeHLA: .*marginal_prob.*.txt</summary>

| Name  | Description |
| ------------- | ------------- |
| locus| HLA locus|
| id | HLA allele code |
| prob | marginal posterior probability cauclated from the HIBAG model |
| freq | frequency of the HLA allele  |
| digits | calculated for X digits  |

</details>

<details> <summary>
imputeHLA: .*overlap_alleles.*.txt</summary>
| Name  | Description |
| ------------- | ------------- |
| gene| HLA locus|
| alleles | HLA allele code; compared allele combination |
| 0%-100%| min, max and quartiles of the % position overlap for N classifiers|
| NClassifiers | Number of classifiers |
</details>
## HLA reference panels

###  Add a reference

To add a reference make an entry into the conf/resources.config


```
	'IKMB' {
		model = "${baseDir}/assets/models/multiethnic_IKMB/model_multiethnic.RData"
		loci = ["A","B","C","DRB1","DQA1","DQB1","DPA1","DPB1","DRB3","DRB4","DRB5"]
	}
  
```

## 
## REFERENCES
Zheng X, Shen J, Cox C, Wakefield JC, Ehm MG, Nelson MR, Weir BS. HIBAG--HLA genotype imputation with attribute bagging. Pharmacogenomics J. 2014 Apr;14(2):192-200. doi: 10.1038/tpj.2013.18. Epub 2013 May 28. PMID: 23712092; PMCID: PMC3772955.









