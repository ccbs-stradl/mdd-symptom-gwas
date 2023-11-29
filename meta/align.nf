/* align sumstats to reference panel */

params.out = "aligned"
params.frq = "*.frq2.gz"
params.daner = "daner*.gz"

workflow {

    FRQ_CH = Channel
        .fromPath(params.frq)
        .collect()

    REF_CH = REF(FRQ_CH)

    DANER_CH = Channel
        .fromPath(params.daner)

    QC_CH = ALIGN(DANER_CH, REF_CH)


}

process REF {

	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

    cpus = 1
    memory = 28.GB
    time = '30m'

    //module '2022:R/4.2.1-foss-2022a'

    input:
    path frq

    output:
    path("ref.rds")

    script:
    """
    #!Rscript
    library(dplyr)
    library(readr)
    library(stringr)

    impute_frq2_files <- str_split("${frq}", pattern=" ")[[1]]

    impute_frq2 <-
    bind_rows(
    lapply(impute_frq2_files,
        function(frq2_file)
            read_table(frq2_file,
                        col_types=cols(SNP = col_character(),
                                        CHR = col_integer(),
                                        POS = col_integer(),
                                        A1 = col_character(),
                                        A2 = col_character(),
                                        FA1 = col_double(),
                                        NCHROBS = col_integer()
    )))) |>
    select(-NCHROBS) |>
    arrange(CHR, POS)

    saveRDS(impute_frq2, "ref.rds")
    """

}

process ALIGN {
    tag "${daner}"
	
	scratch true
	stageInMode 'copy'
	stageOutMode 'copy'

    publishDir params.out, mode: 'copy'

    cpus = 1
    memory = { 16.GB * task.attempt }
    time = '30m'

    //module '2022:R/4.2.1-foss-2022a'

    input:
    each path(daner)
    path(ref)

    output:
    tuple path("*.gz"), path("*.snpcount")

    script:
    """
    #!Rscript
    library(dplyr)
    library(readr)
    library(stringr)
    library(QCGWAS)

    daner <- read_table("${daner}", na = c("NA", "Inf", "nan"))
    impute_frq2 <- readRDS("${ref}")

    frq_a_col <- str_subset(names(daner), "FRQ_A")
    frq_u_col <- str_subset(names(daner), "FRQ_U")

    # check for NCAS/NCON columns
    if(!"NCAS" %in% names(daner)) {
        Nca <- as.numeric(str_extract(frq_a_col, "[[:digit:]]+"))
        Nco <- as.numeric(str_extract(frq_u_col, "[[:digit:]]+"))
        daner_n <- daner |>
            mutate(NCAS = Nca, NCON = Nco)
    } else {
        daner_n <- daner
    }

    # positions that are matchable to the reference
    daner_preharmonise <- daner_n |>
        na.omit() |>
        mutate(BETA = log(OR)) |>
        # merge on chromosome and position
        inner_join(impute_frq2 ,
            by=c('CHR'='CHR', 'BP'='POS'), suffix=c('', '.ref'),
            multiple = "all") |>
        filter((A1 == A1.ref & A2 == A2.ref) | (A1 == A2.ref & A2 == A1.ref)) |>
        mutate(CPID = str_c(CHR, BP, sep = ":")) |>
        filter(!duplicated(CPID))
        
    # QC input
    daner_set <- daner_preharmonise |>
        select(MARKER = SNP.ref, CHR, POSITION = BP, EFFECT_ALL = A1, OTHER_ALL = A2, EFFECT = BETA, STDERR = SE, PVALUE = P, EFF_ALL_FREQ = starts_with("FRQ_U"))

    ref_set <- impute_frq2 |>
        select(SNP, CHR, POSITION = POS, MINOR = A1, MAJOR = A2, MAF = FA1)

    # match and harmonise alleles and effects
    matched_set <- match_alleles(daner_set, ref_set, return_SNPs = TRUE,
                   check_FRQ = TRUE, check_ambiguous = TRUE,
                   delete_mismatches = FALSE, delete_diffEAF = TRUE,
                   threshold_diffEAF = 0.20
                   )

    # harmonise daner to matches
    gwas_matched <- as_tibble(matched_set[c("MARKER", "EFFECT_ALL", "OTHER_ALL", "EFFECT", "EFF_ALL_FREQ")]) |>
        filter(!is.na(EFFECT_ALL), !is.na(OTHER_ALL))

    daner_harmonised <- gwas_matched |>
        inner_join(daner_preharmonise, by=c("MARKER" = "SNP.ref")) |>
        mutate(flip = sign(EFFECT) != sign(BETA)) |>
        transmute(CHR, SNP = MARKER, BP, A1 = EFFECT_ALL, A2 = OTHER_ALL,
        !!frq_a_col := if_else(flip, true = 1 - .data[[frq_a_col]], false = .data[[frq_a_col]]),
        !!frq_u_col := if_else(flip, true = 1 - .data[[frq_u_col]], false = .data[[frq_u_col]]),
        INFO, OR = exp(EFFECT), SE, P, NCAS, NCON)

    write_tsv(daner_harmonised, "${daner.baseName}.align.gz")

    aligned_size <- data.frame(gwas="${daner}", pre_align=nrow(daner), post_align=nrow(daner_harmonised))

    write_tsv(aligned_size, "${daner.baseName}.snpcount")
    """
}