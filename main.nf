#!/usr/bin/env nextflow

/**
===============================
HLA Pipeline
===============================

### Homepage / git
git@github.com:ikmb/HLApipe.git
### Implementation
Implemented in Q1 2020

Author: Frauke Degenhardt, f.degenhardt@ikmb.uni-kiel.de
Author: Mareike Wendorff, m.wendorff@ikmb.uni-kiel.de

**/

params.version = workflow.manifest.version

// Help message
helpMessage = """
===============================================================================
IKMB HLA pipeline | version ${params.version}
===============================================================================
Usage:  nextflow run HLApipePublic --prefix FILE --reference_name REFERENCE --run_name NAME --shapeit SHAPEIT --impute2_reference_dir IMPUTE2_REF_DIR 


--prefix		An input prefix referencing a set of PLINK files
--reference_name 	Name of the reference imputation panel (see below for details)
--run_name		Name prefix + "_" + run_name will be used for the output

Optional parameters:
General: 

--loci      	Loci that should be imputed. Default: As specified in conf/resources.config.

Software/References: 
--shapeit	Path to the SHAPEIT2 executable.
--impute2_ref_dir Path to the IMPUTE2 reference.
--beagle        Location of the .jar file of BEAGLE4.1.
--do_beagle    	This flag is optional and enables phasing using Beagle. Default: false
--ref_1000G     Path to population used for PCA. PLINK files (hg19). Default: 1000G Phase 3 population.
--sample	Path to sample file used for PCA. Default: 1000G Phase 3 population.
--subpop    	Name of a sub population to use. Valid options are: AA, AFR, AMR, CHN, EAS, EUR, GER, IND, IRN, JPN, KOR, MLT. Can be used together with the IKMB reference.

Others: 

--email         Email address to send reports to (enclosed in '')
--outdir        Path to output directory. Default: results.

"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

def summary = [:]

// *****************************
// Define a run name
// *****************************
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

// *****************************
// Define a location for files
// *****************************

params.impute2_reference_dir = "/path/to/impute2/reference/files"
params.shapeit = "/path/to/shapeit2/executable"
params.beagle = "/path/to/beagle/jar"

// *****************************
// Location of preloaded files
// *****************************
params.supplementary_dir = baseDir + "/assets/supplementary"
params.tex_dir = baseDir + "/assets/tex"
params.sample = baseDir + "/assets/samples/1000GP_Phase3_global.cluster"
params.do_beagle = false
params.bin_dir= baseDir + "/bin"

// *****************************
// Initialize other parameters
// *****************************
params.loci = false
params.dict = false
params.splitlnumber=100
params.frqs=""
params.subpop=""
// *****************************
// Sanity checks and validations
// *****************************
if (!params.shapeit || params.shapeit == "/path/to/shapeit2/executable" ) {
	exit 1, "Must provide a path to the shapeit executable. Edit in netflow.config or supply by using (--shapeit)."
} 

if (!params.impute2_reference_dir || params.impute2_reference_dir == "/path/to/impute2/reference/files") {
	exit 1, "Must provide a path to the IMPUTE2 reference. Edit in netflow.config or supply by using (--impute2_reference_dir)."
} 

if (!params.references) {
	exit 1, "No hashmap for references defined for this execution profile, please see documentation for how to create such information."
}


if (!params.references.containsKey(params.reference_name)) {
	exit 1, "Provided an unknown reference name for model selection (--reference_name)."
} else {
	MODEL_NAME = params.references[params.reference_name].model ?: false
}


DICTIONARY = params.dict ?:  params.references[params.reference_name].dict

if (!DICTIONARY) {
	exit 1, "Must provide dictionary for AA and SNP translation. Edit in conf/resources.config."
}


if (!params.sample){
	exit 1, "Must provide a path to the 1000G .sample file. Edit in netflow.config or supply by using (--sample)."
}


if (!params.valid_pops.contains(params.subpop) && params.subpop!=""  && (params.reference_name == "IKMB" || params.reference_name == "IKMB_g" || params.reference_name == "IKMB_1KG") ) {
	exit 1, "Requested for an unknown subpopulation to be analysed (--subpop)"
}
if (params.subpop &&  !(params.reference_name== "IKMB" || params.reference_name == "IKMB_g" || params.reference_name == "IKMB_1KG") ) {
		exit 1, "This reference does not support subpopulations (do not put --subpop)."
}
 
if (params.do_beagle) {
BEAGLE_MAP = file("${baseDir}/assets/beagle_map/plink.chr6.GRCh37.map")
if (!BEAGLE_MAP.exists() ){
	exit 1, "Could not find the beagle map that is supposed to be included with this code base."
}
}

if (params.do_beagle && params.beagle == "/path/to/beagle/jar") {
	exit 1, "Must provide a path to the Beagle jar-file or supply by using (--beagle)."
}

// *****************************
//  Assign loci
// *****************************


if (!params.loci) {
LOCI=params.references[params.reference_name].loci}else{
LOCI= params.loci.tokenize(',')
}


// *******************************
// Define empty channels if needed
// ******************************
if (!params.do_beagle) {
	bimbedfam_beagle = Channel.empty()
}

// ****************************************
// Record analysis parameters for reporting
// ****************************************

summary['runName'] = run_name
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Reference'] = params.reference_name
summary['Model'] = MODEL_NAME
summary['Beagle'] =  params.do_beagle
summary['Analysed loci'] =  LOCI
summary['Dictionary'] = DICTIONARY

if (params.subpop && (params.references == "IKMB" || params.references == "IKMB_g" || params.references == "IKMB_1KG")) {
	summary['SubPopulation'] = params.subpop
}
// ****************************************
// WORKFLOW STARTS HERE
// ****************************************


Channel.fromPath("${params.prefix}.{bed,bim,fam}")
	.ifEmpty { exit 1, "Could not find the specified input..." }
	.toSortedList()
	.set { inputFiles }


// **********
// This process reads data & excludes duplicated sites
// ***********


process readData {

//	publishDir "${params.outdir}/preprocess", mode: 'copy'
	scratch params.scratch
	input:
	set file(bed),file(bim),file(fam) from inputFiles

	output:
	set file(bed_mod),file(bim_mod),file(fam_mod) into (bimbedfam_in, bimbedfam_report)
        set file(bed_phase),file(bim_phase),file(fam_phase) into (bimbedfam_phase)
    	val(base_name) into baseNameTMP
	val(prefix) into basenameRunname
        file("split_files*") into chunks

	script:
	prefix_orig = bed.getBaseName()
	prefix = prefix_orig + "_" + run_name
	base_name= prefix +  "_6_29_34"
	bed_mod = base_name + ".bed"
	bim_mod = base_name + ".bim"
	fam_mod = base_name + ".fam"
	base_name_phase = prefix + ".phase"
	bed_phase = base_name_phase + ".bed"
        bim_phase = base_name_phase + ".bim"
        fam_phase = base_name_phase + ".fam"
        
	"""	
		readData.sh $prefix_orig $prefix ${params.splitlnumber}      
	"""

}



// **********
// This process uses Beagle to fill up missing data in the data set.
// **********

// We either run Beagle, or we just use the bed/bim/fam files as-is
if (params.do_beagle) {

	process runBeagle {
		scratch params.scratch
//		publishDir "${params.outdir}/Beagle", mode: 'copy'

		input:
		set file(bed),file(bim),file(fam) from bimbedfam_in
                val(base_name_old) from baseNameTMP

		output:
		set file(bed_mod),file(bim_mod),file(fam_mod) into bimbedfam_beagle
	        val(base_name) into baseName

		script:
                base_name = base_name_old + "_beagle"
		bed_mod = base_name + ".bed"
		bim_mod = base_name  + ".bim"
		fam_mod = base_name + ".fam"

		"""
			beagle.sh $base_name_old $base_name $BEAGLE_MAP ${params.beagle}
		"""

	}

}  else {

        bimbedfam_in.set{bimbedfam_beagle}
        baseNameTMP.set{baseName}

        
}


process alignToReference {
	scratch params.scratch
	publishDir "${params.outdir}/alignToReference", mode: 'copy', pattern: '*summary*'

    input:
        set file(bed),file(bim),file(fam) from bimbedfam_beagle
        val(base_name) from baseName
	val(basenameRunname) from basenameRunname

    output:
        file(out)
        set file(bed_mod),file(bim_mod),file(fam_mod) into (bimbedfam_impute, bimbedfam_impcombine, bimbedfam_phasecombine)
        val(base_checked) into checkedName
        set file(refA), file(refB), file(ref)  into report_datacheck
      
        
	script:
        base_checked = base_name +  ".refchecked"
        bed_mod = base_checked + ".bed"
        bim_mod = base_checked + ".bim"
        fam_mod = base_checked + ".fam"
	refA = "refchecked_A_" + basenameRunname + ".pdf"
	refB = "refchecked_B_" + basenameRunname + ".pdf"
	ref = "report_refchecked_" + basenameRunname + ".txt"
        out = base_name + ".summary.refchecked.txt"

        def options = ""
        if (params.subpop && (params.references == "IKMB" || params.references == "IKMB_g" || params.references == "IKMB_1KG")) {
                options = params.frqs

        }
        """
        Rscript ${baseDir}/bin/alignToReference.R $MODEL_NAME $base_name ${options}
		mv "refchecked_A.pdf" $refA
		mv "refchecked_B.pdf" $refB
		mv "report.refchecked.txt" $ref
        """
}


// **********
// This process performs HLA imputation
// **********
process imputeHLA {
	scratch params.scratch
	input:
	set file(bed),file(bim),file(fam) from bimbedfam_impute
	each locus from LOCI
	val(checked_name) from checkedName
	
	output:
	file(impute_result) into imputation
	set file(impute_result),val(locus) into imputation_tmp
	file(same_result) into imputation_same

	script:	
	impute_result = "imputation_" + checked_name + "_" + locus + ".RData"
	same_result = locus + "_" + checked_name + ".txt"

	"""
		Rscript ${baseDir}/bin/imputeHLA.R $checked_name $MODEL_NAME $locus 
        Rscript ${baseDir}/bin/imputeHLAhaplo.R $locus $MODEL_NAME $checked_name $checked_name
	"""
}

// **********
// This process combines all HLA loci after imputation
// **********
process imputeHLACombine {
	scratch params.scratch
	publishDir "${params.outdir}/imputeHLACombine", mode: 'copy', pattern: 'imputation*'
	
	input:
	file(singles) from imputation.collect()
	file(singles_same) from imputation_same.collect()
        val(checked_name) from checkedName
	set file(bed),file(bim),file(fam) from bimbedfam_impcombine
	val(basenameRunname) from basenameRunname
	

	output:
        set file(postprob_pdf), file(info), file(unsure_pdf), file(unsure), file(same_result) into report_postprob
	set file(csv_out), file(data), file(info), file(map), file(ped), file(bim), file(bed), file(fam), file(same_result)
        
	script:	
	unsure = "imputation_marginal_prob_" +  checked_name +  ".txt"
	csv_out = "imputation_" + checked_name +  ".csv"
        data = "imputation_" + checked_name + ".RData" 
        info = "imputation_" + checked_name + ".info"
	map =  "imputation_" + checked_name + ".map"
	ped =  "imputation_" + checked_name + ".ped"
	bim =  "imputation_" + checked_name + ".bim"
	bed =  "imputation_" + checked_name + ".bed"
	fam =  "imputation_" + checked_name + ".fam"
	same_result = "imputation_overlap_alleles_" + checked_name + ".txt"
	postprob_pdf = "postprob_" + basenameRunname +".pdf"
	unsure_pdf = "unsure_" + basenameRunname + ".pdf"
	"""
        Rscript ${baseDir}/bin/imputeHLACombineRData.R $singles $checked_name $params.supplementary_dir/location_HLA_genes_hg19.txt $DICTIONARY
		Rscript ${baseDir}/bin/imputeHLACombineCSV.R $singles $checked_name
		Rscript ${baseDir}/bin/imputeHLACombinePLINK.R $checked_name $params.bin_dir/utilityFUNCTIONS.R
                
		cat $singles_same > tmp.txt
                cat <(head -1 tmp.txt) <(grep -v gene tmp.txt) > $same_result
                imputeHLACombineINFO.R tmp.RData $checked_name $data $params.bin_dir/utilityFUNCTIONS.R

		mv "postprob.pdf" $postprob_pdf
		mv "unsure.pdf" $unsure_pdf
	"""	
}




// **********
// This process performs SNP phasing
// **********
process phaseSNPs {

	scratch params.scratch
//	publishDir "${params.outdir}/phasedSNPs", mode: 'copy'
        errorStrategy 'ignore' 
	input:
	set file(bed), file(bim), file(fam) from bimbedfam_phase
	val(checked_name) from checkedName
	val(base_name) from baseName
	each fam_chunk from chunks

	output:
//	set file(haps), file(sample), file(correlation) into phased
        val(phased_name_val) into phased_name 
        file(haps) into phased_haps
 	file(sample) into phased_sample
        file(correlation) into phased_correlation

	script:
	prefix = bed.getBaseName()
        val_chunk = fam_chunk.toString().replaceAll(/.*\//,"")
	haps = checked_name + "_" + val_chunk + ".haps"
	sample = checked_name + "_" + val_chunk + ".sample"
	correlation = checked_name + "_" + val_chunk + ".certainty.all"
	phased_name_val = checked_name + "_" + val_chunk
        """                  
		plink --bfile ${prefix} --keep ${fam_chunk} --geno 0.05 --make-bed --out ${checked_name}"_"${val_chunk}
		phaseSNPs.sh ${checked_name}"_"${val_chunk} ${params.shapeit} ${params.impute2_reference_dir} ${task.cpus}
	#	rm *.log *.sample *.haps
	"""
}

// **********
// This process performs HLA phasing
// **********

process phaseHLA {

	scratch params.scratch
	input:
        each name from phased_name.collect()        
//	set file(haps), file(sample), file(certainty) from phased
        file(haps_tmp) from phased_haps.collect()
        file(sample_tmp) from phased_sample.collect()
        file(correlation_tmp) from phased_correlation.collect()
	set file(imputed), val(locus) from imputation_tmp

	output:
	file(phase_txt) into phasedHLA
	file(phase_info)

	script:
        checked_name=name
	phase_txt = checked_name +  "." + locus +  ".HLA.phased.txt"
	phase_info = checked_name + "." + locus + ".HLA.info.txt"
        haps=name + ".haps"
        sample = name + ".sample"
        certainty = name + ".certainty.all"
	"""
                  
		Rscript ${baseDir}/bin/phaseHLA.R $checked_name $imputed $haps $sample $certainty $MODEL_NAME $locus
	"""
}


// **********
// This process combines loci from HLA phasing
// **********

process phaseHLACombine {

	scratch params.scratch
	publishDir "${params.outdir}/phaseHLAcombine", mode: 'copy', pattern: 'imputation*'

	input:
	file singles from phasedHLA.collect()
	set file(bed),file(bim),file(fam) from bimbedfam_phasecombine
 	val(basenameRunname) from basenameRunname
        val(checked_name) from checkedName

	output:
	file(phased_result)
	set file(phased_pdf), file(phased_result) into report_phasing
        set file(data), file(info), file(map), file(ped), file(bim), file(bed), file(fam)
        file(phased_comb)

	script:
        phased_comb = "imputation_" + checked_name + ".phased.csv"
        phased_result = "imputation_" + checked_name + ".META.PHASING.txt"
        data = "imputation_" + checked_name + ".haplotypes.RData"
        info = "imputation_" + checked_name + ".haplotypes.info"
	phased_pdf = "phased_" + basenameRunname + ".pdf"
        map =  "imputation_" + checked_name + ".haplotypes.map"
        ped =  "imputation_" + checked_name + ".haplotypes.ped"
        bim =  "imputation_" + checked_name + ".haplotypes.bim"
        bed =  "imputation_" + checked_name + ".haplotypes.bed"
        fam =  "imputation_" + checked_name + ".haplotypes.fam"
	"""
		Rscript ${baseDir}/bin/phaseHLACombineMETA.R ${singles} $phased_result
		Rscript ${baseDir}/bin/phaseHLACombineRData.R $checked_name
		Rscript ${baseDir}/bin/phaseHLACombinePLINK.R $checked_name $params.bin_dir/utilityFUNCTIONS.R
        Rscript ${baseDir}/bin/phaseHLACombineCSV.R $data $phased_comb
        Rscript ${baseDir}/bin/phaseHLACombineINFO.R $checked_name $data $params.bin_dir/utilityFUNCTIONS.R
		mv "phased.pdf" $phased_pdf
	"""
}



// **********
// This process writes an automated report about the data and the imputation quality
// **********

process report {

	scratch params.scratch
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    set file(figrefcheckA),file(figrefcheckAB), file(tabrefcheck) from report_datacheck
    set file(postprob), file(info), file(unsure_pdf), file(unsure_txt), file(haplo) from report_postprob
    file(phased) from report_phasing
    set file(bim), file(bed), file(fam) from bimbedfam_report
    val(checked_name) from checkedName
    val(base_name) from baseName
    val(basenameRunname) from basenameRunname
	
    output:
    file output
 
    script:
    rootname = bed.getBaseName()
    output =  "report_" + basenameRunname +  ".html"
    """
		plink --bfile $rootname --missing --out $rootname
        cp -r ${params.tex_dir}/* .
        R -e 'rmarkdown::render("report.Rmd", output_file="${output}", params = list(rootname="${rootname}", checked_name="${checked_name}",  pop="${params.subpop}", model="${MODEL_NAME}", shapeit="${params.shapeit}", modules="${LOADEDMODULES}", basenamerunname="${basenameRunname}"))'
    """
}
