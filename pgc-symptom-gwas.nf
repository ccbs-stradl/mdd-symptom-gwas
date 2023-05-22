
params.bfile = "pgc_gwas/cobg_dir_genome_wide/*.hg19.ch.fl.bgn.{bed,bim,fam}"
params.pheno = "pgc_gwas/mddw2_symptoms_cases.pheno"
params.covar = "pgc_gwas/mddw2.covar"
params.keep = "pgc_gwas/cases.id"
params.remove = "pgc_gwas/959_PGC_UKB_overlap.id"

workflow {

    BFILE_CH = Channel
        .fromFilePairs(params.bfile, size: 3) { file -> file.simpleName }
        
    PHENO_CH = Channel
        .fromPath(params.pheno)

    COVAR_CH = Channel
        .fromPath(params.covar)

    KEEP_CH = Channel
        .fromPath(params.keep)

    REMOVE_CH = Channel
        .fromPath(params.remove)

    INPUTS_CH = BFILE_CH
        .combine(PHENO_CH)
        .combine(COVAR_CH)
        .combine(KEEP_CH)
        .combine(REMOVE_CH)

    GWAS_CH = GWAS(INPUTS_CH)
        .multiMap { it ->
        glm: it[0]
        log: it[1]
    }

    GLM_CH = GWAS_CH.glm
        .flatten()

    DANER_CH = DANER(GLM_CH)


}

process GWAS {
    tag "${cohort}"
	
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

    publishDir "sumstats/PGC/cohorts", pattern: "*.log"

    cpus = 8
    memory = 16.GB
    time = '6h'

    module '2022:PLINK/2.00a3.6-GCC-11.3.0'

    input:
    tuple val(cohort), path(bfile), path(pheno), path(covar), path(keep), path(remove)

    output:
    tuple path("*.gz"), path("*.log")

    script:
    """
    plink2 \
    --bfile ${bfile[0].baseName} \
    --pheno ${pheno} \
    --covar ${covar} \
    --covar-name C1,C2,C3,C4,C5,C6 \
    --keep ${keep} \
    --remove ${remove} \
    --mac 20 \
    --glm 'hide-covar' 'skip-invalid-pheno' cols=chrom,pos,ax,a1freqcc,totallelecc,machr2,firth,test,nobs,orbeta,se,p \
    --covar-variance-standardize \
    --out ${cohort} \
    --threads ${task.cpus} \
    --memory ${task.memory.bytes.intdiv(1000000)}
	
	gzip *.hybrid
    """
}

process DANER {
    tag "${sumstats}"
	
	publishDir "sumstats/PGC/cohorts"

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