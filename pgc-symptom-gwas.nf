
params.bfile = "pgc_gwas/cobg_dir_genome_wide/*.hg19.ch.fl.bgn.{bed,bim,fam}"
params.pheno = "pgc_gwas/mddw2_symptoms.pheno"
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
}

process GWAS {
    tag "${cohort}"

    cpus = 8
    memory = 16.GB
    time = '3h'

    module '2022:PLINK/2.00a3.6-GCC-11.3.0'

    input:
    tuple val(cohort), path(bfile), path(pheno), path(covar), path(keep), path(remove)

    output:
    tuple path("*.zst"), path("*.log")

    script:
    """
    plink2 \
    --bfile ${bfile[0].baseName} \
    --pheno ${pheno} \
    --covar ${covar} \
    --covar-name C1,C2,C3,C4,C5,C6 \
    --keep ${keep} \
    --mac 20 \
    --remove ${remove} \
    --glm 'zs' 'hide-covar' 'skip-invalid-pheno' cols=chrom,pos,ax,a1freqcc,totallelecc,machr2,firth,test,nobs,orbeta,se,p \
    --covar-variance-standardize \
    --out ${cohort} \
    --threads ${task.cpus} \
    --memory ${task.memory.bytes.intdiv(1000000)}
    """
}