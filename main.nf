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
Usage:  nextflow run HLApipePublic --prefix FILE --reference_name REFERENCE --run_name NAME--shapeit SHAPEIT --impute2_reference_dir IMPUTE2_REF_DIR 


--prefix		An input prefix referencing a set of PLINK files
--reference_name 	Name of the reference imputation panel (see below for details)

Optional parameters:
General: 

--assembly  	Build of your data (default hg19; also supported hg18, data are then lifted)
--loci      	Loci that should be imputed. Default: As specified in conf/resources.config.
--splitlnnumber Optional split of data for the phasing step. Number of samples in each split.
 
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
params.liftover_file = baseDir + "/assets/liftover/hg18ToHg19.over.chain.gz"
params.tex_dir = baseDir + "/assets/tex"

params.sample = baseDir + "/assets/samples/1000GP_Phase3_global.cluster"
params.ref_1000G = baseDir + "/assets/variants/chr1_22_0.05"
params.do_beagle = false


// *****************************
// Initialize other parameters
// *****************************
params.assembly = "hg19"
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

SHAPEIT=file("${params.shapeit}")
if (!SHAPEIT.exists() ) {
        exit 1, "Shapit executable not found at path given by --shapeit. Please check the path."}


if (!params.impute2_reference_dir || params.impute2_reference_dir == "/path/to/impute2/reference/files") {
	exit 1, "Must provide a path to the IMPUTE2 reference. Edit in netflow.config or supply by using (--impute2_reference_dir)."
} 
IMPUTEFILE=file("${params.impute2_reference_dir}/1000GP_Phase3_chr6.hap.gz")
if(!IMPUTEFILE.exists()){
    exit 1, "Reference file 1000GP_Phase3_chr6.hap.gz not found at PATH declared with --impute2_reference_dir. Please check the path and existence of the file."}

if (!params.valid_assembly.contains(params.assembly) ) {
		exit 1, "Requested for an unknown assembly (--assembly)"
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

if (!params.ref_1000G){
	exit 1, "Must provide a path to the 1000G .bim/.bed/.fam. Edit in netflow.config or supply by using (--ref_1000G)."
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



LOCI = params.loci ?:  params.references[params.reference_name].loci

if (!LOCI) {
	exit 1, "Must provide loci to be analysed. Default should be in the conf/rescources.config or specified using (--loci)"
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
summary['Assembly'] = params.assembly
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
// This process performs a liftover or prepares input data.
// **********
if (params.assembly != "hg19") {

	process liftover {
		publishDir "${params.outdir}/liftover", mode: 'copy'

		input:
		set file(bed),file(bim),file(fam) from inputFiles

		output:
		set file(bed_lifted),file(bim_lifted),file(fam_lifted) into (inputPrep)
		file(unmSNPs)
		file(unmCHR) into report_liftover
	 
		script:

		prefix = bed.getBaseName()
		base_name = prefix + "_" + run_name + "_hg19"
		bed_lifted = base_name + ".bed"
		bim_lifted = base_name + ".bim"
		fam_lifted = base_name + ".fam"
		
		unmSNPs = "unmappedSNPs_" + base_name + ".txt"
		unmCHR = "unmappedCHR_" + base_name + ".txt"

		"""
		liftover_pos.sh $bim 
		liftOver prelift.pos ${params.liftover_file} postlift.pos unMapped.txt -positions
		liftover.R postlift.pos unMapped.txt $prefix $base_name
		"""

	}
}  else {


	process noliftover {

		input:
		set file(bed),file(bim),file(fam) from inputFiles
        
		output:
		file("noliftover.txt") into report_liftover
	        set file(bed_out), file(bim_out), file(fam_out) into inputPrep
		script:

		prefix = bed.getBaseName()
                base_name = prefix + "_" + run_name 
                bed_out = base_name + ".bed"
                bim_out = base_name + ".bim"
                fam_out = base_name + ".fam"

	       	"""
		plink --bfile $prefix --make-bed --out $base_name
		echo "" > noliftover.txt
		
		"""
	}	        
}


// **********
// This process reads data from liftover, excludes duplicated sites (the duplicate site with less data present)
// and extracts chromosome 6.
//
// **********


process readData {

//	publishDir "${params.outdir}/preprocess", mode: 'copy'

	input:
	set file(bed),file(bim),file(fam) from inputPrep

	output:
	set file(bed_mod),file(bim_mod),file(fam_mod) into (bimbedfam_in, bimbedfam_report)
        set file(bed_qc),file(bim_qc),file(fam_qc) into (bimbedfam_qc)
        set file(bed_phase),file(bim_phase),file(fam_phase) into (bimbedfam_phase)
    	val(base_name) into baseNameTMP
	val(prefix) into basenameRunname
        file("split_files*") into chunks
	script:

	prefix = bed.getBaseName()
	base_name= prefix +  "_6_29_34"
	bed_mod = base_name + ".bed"
	bim_mod = base_name + ".bim"
	fam_mod = base_name + ".fam"
	base_name_qc = prefix + ".qc"
	bed_qc = base_name_qc + ".bed"
	bim_qc = base_name_qc + ".bim"
	fam_qc = base_name_qc + ".fam"
	base_name_phase = prefix + ".phase"
	bed_phase = base_name_phase + ".bed"
        bim_phase = base_name_phase + ".bim"
        fam_phase = base_name_phase + ".fam"

	"""
#		check_build.R ${params.ref_1000G} $prefix
		duplicates.R $prefix
		plink --bfile $prefix --exclude exclude_duplicates.txt --make-bed --out $base_name_qc
         

                
		plink --allow-no-sex --bfile $base_name_qc  --chr 6 --from-mb 25 --to-mb 35 --make-bed --out $base_name_phase
              

                split -l ${params.splitlnumber}  -d $fam_phase "split_files"
                


   
           	plink --allow-no-sex --bfile $base_name_qc --make-bed --out $base_name_qc
		plink --allow-no-sex --bfile $prefix  --chr 6 --from-mb 29 --to-mb 34 --make-bed --out $base_name

               
	"""

}



// **********
// This process uses Beagle to fill up missing data in the data set.
// **********

// We either run Beagle, or we just use the bed/bim/fam files as-is
if (params.do_beagle) {

	process runBeagle {

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

        input:
        set file(bed),file(bim),file(fam) from bimbedfam_beagle
        val(base_name) from baseName
	val(basenameRunname) from basenameRunname

        output:
        set file(bed_mod),file(bim_mod),file(fam_mod) into (bimbedfam_impute, bimbedfam_impcombine, bimbedfam_phasecombine)
        val(base_checked) into checkedName
        set file(refA), file(refB), file(ref)  into report_datacheck

        script:
        base_checked = base_name +  ".refchecked"
        bed_mod = base_checked + ".bed"
        bim_mod = base_checked + ".bim"
        fam_mod = base_checked + ".fam"
	refA = "refchecked_A_" + basenameRunname + ".png"
	refB = "refchecked_B_" + basenameRunname + ".png"
	ref = "report_refchecked_" + basenameRunname + ".txt"


         def options = ""
        if (params.subpop && (params.references == "IKMB" || params.references == "IKMB_g" || params.references == "IKMB_1KG")) {
                options = params.frqs

        }
        """
                match_study_to_reference.R $MODEL_NAME $base_name ${options}
		mv "refchecked_A.png" $refA
		mv "refchecked_B.png" $refB
		mv "report.refchecked.txt" $ref
        """
}


// **********
// This process performs HLA imputation
// **********
process imputeHLA {

	input:
	set file(bed),file(bim),file(fam) from bimbedfam_impute
	each locus from LOCI
	val(checked_name) from checkedName
	
	output:
	file(impute_result) into imputation
	set file(impute_result),val(locus) into imp_HLA

	script:	
	impute_result = "imputation_" + checked_name + "_" + locus + ".RData"

	"""
		impute.R $checked_name ${params.assembly} $MODEL_NAME $locus
	"""
}

// **********
// This process combines all HLA loci after imputation
// **********
process imputeHLACombine {

	publishDir "${params.outdir}/imputeHLA", mode: 'copy'
	
	input:
	file(singles) from imputation.collect()
        val(checked_name) from checkedName
	set file(bed),file(bim),file(fam) from bimbedfam_impcombine
	val(basenameRunname) from basenameRunname
	

	output:
	file(result) into imputation_all
        set file(postprob_png), file(unsure_png), file(unsure) into report_postprob
	set file(csv_out), file(assoc_data), file(assoc_af), file(map), file(ped), file(bim), file(bed), file(fam)

	script:	
	result = "imputation_" + checked_name + ".RData"
	unsure = "marginal_prob_" +  checked_name +  ".txt"
	csv_out = "imputation_" + checked_name +  ".csv" 
        assoc_data = "imputation_" + checked_name + ".data"
        assoc_af = "imputation_" + checked_name + ".af"
        nuc_data = "imputation_" + checked_name + ".nuc.data"
        prot_data =  "imputation_" + checked_name + ".prot.data"
	map =  "imputation_" + checked_name + ".map"
	ped =  "imputation_" + checked_name + ".ped"
	bim =  "imputation_" + checked_name + ".bim"
	bed =  "imputation_" + checked_name + ".bed"
	fam =  "imputation_" + checked_name + ".fam"
	postprob_png = "postprob_" + basenameRunname +".png"
	unsure_png = "unsure_" + basenameRunname + ".png"
	"""
		combine_loci_imputation.R $singles $result
		get_marginal_prob.R $result $checked_name
		postprob.R $result
		make_output.R $result $checked_name
		prepare_single.R $result $checked_name $params.supplementary_dir/location_HLA_genes_hg19.txt $DICTIONARY
		util.R $checked_name
		mv "postprob.png" $postprob_png
		mv "unsure.png" $unsure_png
	"""	
}




// **********
// This process performs SNP phasing
// **********
process phaseSNPs {

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
		phase_global_SNPs_parallel.sh ${checked_name}"_"${val_chunk} ${params.shapeit} ${params.impute2_reference_dir} ${task.cpus}
	#	rm *.log *.sample *.haps
	"""
}

// **********
// This process performs HLA phasing
// **********

process phaseHLA {

//	publishDir "${params.outdir}/phasedHLA/Loci", mode: 'copy'

	input:
        each name from phased_name.collect()        
//	set file(haps), file(sample), file(certainty) from phased
        file(haps_tmp) from phased_haps.collect()
        file(sample_tmp) from phased_sample.collect()
        file(correlation_tmp) from phased_correlation.collect()
	set file(imputed), val(locus) from imp_HLA

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
                  
		phase.R $checked_name $imputed $haps $sample $certainty $MODEL_NAME $locus
	"""
}


// **********
// This process combines loci from HLA phasing
// **********

process phaseHLACombine {

	publishDir "${params.outdir}/phasedHLA", mode: 'copy'

	input:
	file singles from phasedHLA.collect()
	set file(bed),file(bim),file(fam) from bimbedfam_phasecombine
 	val(basenameRunname) from basenameRunname
        val(checked_name) from checkedName
	output:
	file(phased_result)
	file(phased_png) into report_phasing
        file(assoc_data_hapl)
        file(assoc_af_hapl)
        file(phased_comb)
	script:
        phased_comb = checked_name + ".phased_overview.txt"
        phased_result = checked_name + ".HLA.all.phased.txt"
        assoc_data_hapl = "imputation_" + checked_name + ".hapl.data"
        assoc_af_hapl = "imputation_" + checked_name + ".hapl.af"
	phased_png = "phased_" + basenameRunname + ".png"
	"""
		combine_loci_phasing.R ${singles} $phased_result
		prepare_haplo.R $checked_name
                util.R $checked_name
		mv "phased.png" $phased_png
	"""
}


// **********
// This process writes an automated report about the data and the imputation quality
// **********

process report {
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    file liftover from report_liftover
    set file(figrefcheckA),file(figrefcheckAB), file(tabrefcheck) from report_datacheck
    set file(postprob), file(unsurepng), file(unsuretxt) from report_postprob
    file(phased) from report_phasing
    set file(bim), file(bed), file(fam) from inputPrep
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
        R -e 'rmarkdown::render("report.Rmd", output_file="${output}", params = list(rootname="${rootname}", checked_name="${checked_name}",  pop="${params.subpop}", model="${MODEL_NAME}", shapeit="${params.shapeit}", modules="${LOADEDMODULES}", liftover="${liftover}", assembly="${params.assembly}", basenamerunname="${basenameRunname}"))'
    """
}
