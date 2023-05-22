/* align sumstats to reference panel */

params.out = "meta/aligned"
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
    memory = 16.GB
    time = '30m'

    module '2022:R/4.2.1-foss-2022a'

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
    select(-FA1, -NCHROBS) |>
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
    memory = 16.GB
    time = '30m'

    module '2022:R/4.2.1-foss-2022a'

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

    daner <- read_table("${daner}")
    impute_frq2 <- readRDS("${ref}")

    # merge on chromosome and position
    daner_ref <- daner |>
    left_join(impute_frq2 ,
            by=c('CHR'='CHR', 'BP'='POS'), suffix=c('', '.ref')) |>
    # keep SNPs where alleles match
    filter((A1 == A1.ref & A2 == A2.ref) | (A1 == A2.ref & A2 == A1.ref)) |>
    # user ref markername
    mutate(SNP=SNP.ref) |>
    select(-ends_with('.ref')) |>
    arrange(CHR, BP)

    write_tsv(daner_ref, "${daner.baseName}.align.gz")

    aligned_size <- data.frame(gwas="${daner}", pre_align=nrow(daner), post_align=nrow(daner_ref))

    write_tsv(aligned_size, "${daner.baseName}.snpcount")
    """
}