
params.bfile = "pgc_gwas/cobg_dir_genome_wide/*.hg19.ch.fl.bgn.{bed,bim,fam}"
params.pheno = "pgc_gwas/mddw2_symptoms.pheno"
params.covar = "pgc_gwas/mddw2.covar"
params.keep = "pgc_gwas/cases.id"
params.remove = "pgc_gwas/959_PGC_UKB_overlap.id"
params.minCount = "50"

workflow {

    BFILE_CH = Channel
        .fromFilePairs(params.bfile, size: 3) { file -> file.simpleName }
        
    PHENO_CH = Channel
        .fromPath(params.pheno, checkIfExists: true)
		
	SYMPTOMS_CH = Channel
		.value(["MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD5a", "MDD6", "MDD7", "MDD8", "MDD9"])
		.collect()

    COVAR_CH = Channel
        .fromPath(params.covar, checkIfExists: true)

    KEEP_CH = Channel
        .fromPath(params.keep, checkIfExists: true)

    REMOVE_CH = Channel
        .fromPath(params.remove, checkIfExists: true)

    INPUTS_CH = BFILE_CH
        .combine(PHENO_CH)
        .combine(COVAR_CH)
        .combine(KEEP_CH)
        .combine(REMOVE_CH)
		
	minCount = params.minCount.toInteger()
	
	// calculate minimum case/control counts, then filter symptoms with more than params.minCount
	COUNTS_CH = COUNTS(INPUTS_CH)
		.splitCsv(header: true)
		.filter { it -> it[1].case.toInteger() >= minCount & it[1].control.toInteger() >= minCount }
		.map { it -> [it[0], it[1].symptom] }
	
	// merge inputs with symptoms to analyse
	ANALYZE_CH = INPUTS_CH
		.cross(COUNTS_CH)
		.map { it -> [it[0][0], it[0][1], it[0][2], it[0][3], it[0][4], it[0][5], it[1][1]] }

    GWAS_CH = GWAS(ANALYZE_CH)
        .multiMap { it ->
        glm: it[0]
        log: it[1]
    }

	GLM_CH = GWAS_CH.glm
        .flatten()
	
    DANER_CH = DANER(GLM_CH)
	
}

// check symptom present/absent counts
process COUNTS {
	tag "${cohort}"
	
	executor 'local'
	
	module '2022:R/4.2.1-foss-2022a'
	
	input:
	tuple val(cohort), path(bfile), path(pheno), path(covar), path(keep), path(remove)
	
	output:
	tuple val(cohort), path("counts.csv")
	
	script:
	"""
	#!Rscript
	library(dplyr)
	library(readr)
	library(tidyr)
	
	fam <- read_table("${bfile[0].baseName}.fam", col_names=c("FID", "IID", "father", "mother", "sex", "pheno"))
	
	pheno <- read_table("${pheno}")
	covar <- read_table("${covar}")
	keep <- read_table("${keep}", col_names=c("FID", "IID"))
	remove <- read_table("${remove}", col_names=c("FID", "IID"))
	
	symptom_counts <- 
	pheno |>
	filter(FID %in% pull(fam, FID),
		   FID %in% pull(covar, FID),
	       FID %in% pull(keep, FID),
		   !FID %in% pull(remove, FID)) |>
	pivot_longer(cols=MDD1:MDD9, names_to = "symptom", values_to = "status") |>
	filter(status %in% 1:2) |>
	mutate(status=recode(status, "1"="case", "2"="control")) |>
	group_by(symptom, status) |>
	count() |>
	pivot_wider(names_from = "status", values_from="n", values_fill=0)
	
	write_csv(symptom_counts, "counts.csv")
	"""
}

// Symptom GWAS
process GWAS {
    tag "${cohort}"
	
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

    publishDir "sumstats/PGC/cohorts", pattern: "*.log", mode: 'copy'

    cpus = 16
    memory = 28.GB
    time = '30m'

    module '2022:PLINK/2.00a3.6-GCC-11.3.0'

    input:
    tuple val(cohort), path(bfile), path(pheno), path(covar), path(keep), path(remove), val(symptom)

    output:
    tuple path("*.gz"), path("*.log")

    script:
    """
    plink2 \
    --bfile ${bfile[0].baseName} \
    --pheno ${pheno} \
	--pheno-name ${symptom} \
    --covar ${covar} \
    --covar-name C1,C2,C3 \
    --keep ${keep} \
    --remove ${remove} \
    --mac 20 \
    --glm 'hide-covar' cols=chrom,pos,ax,a1freqcc,totallelecc,machr2,firth,test,nobs,orbeta,se,p \
    --covar-variance-standardize \
    --out ${cohort}_${symptom} \
    --threads ${task.cpus - 1} \
    --memory ${task.memory.bytes.intdiv(1000000)}
	
	gzip *.hybrid
    """
}

// convert Plink output to daner format
process DANER {
    tag "${sumstats}"
	
	publishDir "sumstats/PGC/cohorts", mode: 'copy'

    cpus = 1
    memory = 16.GB
    time = '30m'

    module '2022:R/4.2.1-foss-2022a'

    input:
    path(sumstats)

    output:
    path("*.gz")

    script:
    """
    #!Rscript

    library(readr)
    library(dplyr)
    library(stringr)

    sumstats_gz <- "${sumstats}"

    analysis <- str_split(sumstats_gz, pattern=fixed("."))[[1]]
    cohort <- str_split(analysis[1], pattern=fixed("_"))[[1]][2]
    symptom <- analysis[2]

    daner_gz <- paste(paste("daner", symptom, cohort, sep="_"), "gz", sep=".")

    sumstats <- read_table(sumstats_gz)

    n_cases <- max(pull(sumstats, CASE_ALLELE_CT))/2
    n_controls <- max(pull(sumstats, CTRL_ALLELE_CT))/2

    frq_a_col <- paste('FRQ', 'A', n_cases, sep='_')
    frq_u_col <- paste('FRQ', 'U', n_controls, sep='_')

    daner <- sumstats |>
    transmute(
    CHR=`#CHROM`,
    SNP=ID,
    BP=POS,
    A1,
    A2=AX,
    !!frq_a_col:=A1_CASE_FREQ,
    !!frq_u_col:=A1_CTRL_FREQ,
    INFO=MACH_R2,
    OR,
    SE=`LOG(OR)_SE`,
    P,
    NCAS=CASE_ALLELE_CT/2,
    NCON=CTRL_ALLELE_CT/2
    )

    write_tsv(daner, daner_gz)
    """
}