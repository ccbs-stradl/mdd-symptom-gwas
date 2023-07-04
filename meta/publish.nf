/* PGC sumstats format for posting */

nextflow.enable.dsl = 2

params.sumstats = "inputs/*.{gz,xls}"
params.glue = "pgc.glue"
params.cff = "../CITATION.cff"
params.contig = "http://fileserve.mrcieu.ac.uk/ref/2.8/b37/human_g1k_v37.fasta.fai"
params.cohorts = "cohort_alignment.txt"
params.prefix = "mdd_symptoms_2023"
params.acknowledgments = "The PGC has received funding from the US National Institute of Mental Health (5 U01MH109528-04). Statistical analyses were carried out on the Genetic Cluster Computer (http://www.geneticcluster.org) hosted by SURFsara and financially supported by the Netherlands Scientific Organization (NWO 480-05-003) along with a supplement from the Dutch Brain Foundation and the VU University Amsterdam."
params.analyst = "Mark J Adams"

workflow {
    
    SUMSTATS_CH = Channel
        .fromFilePairs(params.sumstats, size: 2)

    GLUE_CH = Channel
        .fromPath(params.glue)

    CFF_CH = Channel
        .fromPath(params.cff)

    CONTIG_CH = Channel
        .fromPath(params.contig)

    COHORT_CH = Channel
        .fromPath(params.cohorts)
    
    PREFIX_CH = Channel
        .of(params.prefix)

    ACKNOWLEDGMENTS_CH = Channel
        .of(params.acknowledgments)

    ANALYST_CH = Channel
        .of(params.analyst)

    INPUTS_CH = SUMSTATS_CH
        .combine(GLUE_CH)
        .combine(CFF_CH)
        .combine(COHORT_CH)
        .combine(CONTIG_CH)
        .combine(PREFIX_CH)
        .combine(ACKNOWLEDGMENTS_CH)
        .combine(ANALYST_CH)

    PGC_CH = PGC(INPUTS_CH)
    
    GZIP(PGC_CH)

}

process PGC {
    tag "${dataset}"

    cpus = 1
    memory = 6.GB
    time = '10m'

    input:
    tuple val(dataset), path(sumstats), path(glue), path(cff), path(cohorts), path(contig), val(prefix), val(acknowledgments), val(analyst)

    output:
    path("${prefix}-${dataset}.txt")

    script:
    """
    #!Rscript

    library(dplyr)
    library(readr)
    library(readxl)
    library(stringr)
    library(lubridate)
    library(yaml)

    daner <- read_tsv("${dataset}.gz")
    basic <- read_excel("${dataset}.xls")

    filedate <- format(now(), "%Y-%Om-%d") 

    analyst <- "${analyst}"

    # read citation information
    cff <- read_yaml("${cff}")

    # cohorts filename information
    cohorts <- read_tsv("${cohorts}")

    # header template
    header_glues <- readLines("${glue}")

    # fasta contig
    fasta_fai <- read_table("${contig}", col_names=c('NAME', 'LENGTH', 'OFFSET', 'LINEBASES', 'LINEWIDTH' ))
    
    doi <- cff[['identifiers']][[1]][['value']]
    sumstats_url = cff[['url']]
    code_url = cff[['repository-code']]
    abstract = cff[['abstract']]

    acknowledgments <- "${acknowledgments}"
    analysis_version <- "${prefix}"
    dataset <- "${dataset}"
    symptom <- str_split_1(dataset, pattern='-')[2]

    # variants
    variants <- nrow(daner)

    # total sample sizes
    ncase <- basic |> filter(Dataset == "SUM") |> pull(N_cases)
    ncontrol <- basic |> filter(Dataset == "SUM") |> pull(N_controls)
    neff <- basic |> filter(Dataset == "SUM") |> mutate(N_eff=2*N_eff_half) |> pull(N_eff)
    ntrio <- 0

    # per-cohort-info
    ##cohortList="{}"
    ##sumstatsVersionByCohort="{}"
    basic_cohorts <- basic |>
        mutate(fileset=str_replace(Dataset, ".align", "")) |>
        inner_join(mutate(cohorts, fileset=str_replace(filename, "daner_", "")), by="fileset") |>
        mutate(cohort=if_else(cohort == "PGC",
                             true=str_match(fileset, pattern="MDD.+_([a-z0-9]+)")[,2],
                             false=cohort),
                core=if_else(cohort == "PGC", true=TRUE, false=FALSE),
                neff=2*N_eff_half)

    ncohort <- nrow(basic_cohorts)
    cohort_list <- basic_cohorts |> pull(cohort) |> paste(collapse=",")
    versions_list <- basic_cohorts |> pull(filename) |> paste(collapse=",")
    cases_by_cohort <- basic_cohorts |> pull(N_cases) |> paste(collapse=",")
    controls_by_cohort <- basic_cohorts |> pull(N_controls) |> paste(collapse=",")
    neff_by_cohort <- basic_cohorts |> pull(neff) |> paste(collapse=",")
    trios_by_cohort <- rep(0, ncohort) |> paste(collapse=",")
    snps_by_cohort <- basic_cohorts |> pull(`N-SNPs`) |> paste(collapse=",")
    processed_by_core  <- basic_cohorts |> pull(core) |> paste(collapse=",")

    # format contig
    daner_chr <- daner %>% select(CHR) |>
     distinct(CHR) |>
     mutate(CHR=as.character(CHR)) |>
     pull(CHR)

    # map fai "X" -> "23"
    contig_lengths <- fasta_fai %>%
    mutate(CHR=case_when(NAME == "X" ~ "23",
                        NAME == "XY" ~ "24",
                        NAME == "MT" ~ "25",
                        TRUE ~ NAME)) %>%
    filter(CHR %in% daner_chr) %>%
    mutate(contig=str_glue("##contig=<ID={CHR},length={LENGTH}>"))

    contig <- paste(pull(contig_lengths, contig), collapse='\\n')

    # header
    headers <- sapply(header_glues, str_glue)
    header <- paste(headers, collapse='\\n')

    # format sumstats
    pgc_sumstats <- daner %>%
        arrange(CHR, BP) |>
        select(CHR, BP, SNP, A1, A2, OR, SE,
            FRQ_A=starts_with('FRQ_A'), FRQ_U=starts_with('FRQ_U'), P,
            INFO, Neff_half, Nca, Nco) %>%
        transmute(`#CHROM`=CHR, POS=BP, ID=SNP, A1=A1, A2=A2,
                BETA=log(OR), SE, PVAL=P, FCAS=FRQ_A, FCON=FRQ_U,
                IMPINFO=INFO, NEFF=2*Neff_half, NCAS=Nca, NCON=Nco)

    out <- "${prefix}-${dataset}.txt"

    cat(header, "\\n", file=out)

    write_tsv(pgc_sumstats, file=out, append=TRUE, col_names=TRUE)
    """

}

process GZIP {
    tag "${sumstats}"

    publishDir "publish", mode: 'copy'

    cpus = 1
    memory = 1.GB
    time = '10m'

    input:
    path(sumstats)

    output:
    path("*.gz")

    script:
    """
    gzip -c ${sumstats} > ${sumstats}.gz
    """

}