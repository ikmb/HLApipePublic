![](../images/ikmb_bfx_logo.png)

# Usage

## Basic execution

The typical way to launch this pipeline would be:

```
nextflow run ikmb/HLApipe --run_name <some_name> --prefix '/path/to/datafile' --beagle
```

### `--run_name`

A custom name for this analysis, will be used to name output files.

May be deprecated to allow for processing of multiple samples!

### `--prefix`

This option points to the root name of a set of bed/bim/fam files. So if you have the following data:

```
ls /path/to/data
  ..my_data.bed
  ..my_data.bim
  ..my_data.fam
```

then the way to run the pipeline on these data would be:

```
nextflow run ikmb/HLApipe --prefix '/path/to/data/my_data'
```

### `--beagle`

This flag is optional and enables phasing using Beagle. 

Default: false


### `--subpop`

Name of a sub population to use. Valid options are: AA, AFR, AMR, CHN, EAS, EUR, GER, IND, IRN, JPN, KOR, MLT

Default: false

### `--chip`

The name of a genotyping chip the input data is based on. This option is overriden for certain reference names where a pre-defined chip
configuration has been set. 

Currently configured options are:

- immunochip

- affymetrix_axiom

- gsa

### `--reference_name`

The name of the reference data set. This option overrides the `--chip` option for certain reference names.

Currently configured options are:

- IKMB

- IKMB_1KG

- IKMB_g

- 1KG






