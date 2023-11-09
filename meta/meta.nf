/* LOO meta-analysis */

nextflow.enable.dsl = 2

params.sumstats = "aligned/*.gz"
params.datasets = "*.meta"

workflow {

    // sumstats in daner format
    DANER_CH = Channel
        .fromPath(params.sumstats)
        .map { it -> [it.name, it]}

    // lists of datasets to meta analyse
    DATASETS_CH = Channel
        .fromPath(params.datasets)
        .map { it -> [it, it.readLines()]}
        .transpose()
        .map { it -> [it[1], it[0]]}

    // unzip daner files so they can be streamed
    // harmonise allele frequency column names
    SUMSTATS_CH = SUMSTATS(DANER_CH)
    
    // match up datasets with the sumstats files they specify
    DATASETS_SUMSTATS_CH =
    DATASETS_CH
        .join(SUMSTATS_CH)
        .map { it -> [it[1], it[2]]}
        .groupTuple()

    META_CH = META(DATASETS_SUMSTATS_CH)
        .view()
    
}

process SUMSTATS {
    tag "${dataset}"

    cpus = 1
    memory = 16.GB
    time = '10m'

    input:
    tuple val(dataset), path(daner)

    output:
    tuple val(dataset), path("${daner.baseName}")

    script:
    """
    #!Rscript

    library(dplyr)
    library(readr)
    library(stringr)

    daner <- read_tsv("${daner}")

    col_names <- names(daner)

    if('NCAS' %in% col_names) {
        daner_n <- daner
    } else {
        frq_a_col <- str_subset(col_names, "FRQ_A")
	    Nca <- as.numeric(str_extract(frq_a_col, "[[:digit:]]+"))
	    frq_u_col <- str_subset(col_names, "FRQ_U")
	    Nco <- as.numeric(str_extract(frq_u_col, "[[:digit:]]+"))

        daner_n <- daner |>
            mutate(NCAS = Nca, NCON = Nco)
    }

    sumstats <- daner_n |>
        select(CHR, SNP, BP, A1, A2,
               FRQ_A = starts_with("FRQ_A"), FRQ_U = starts_with("FRQ_U"),
               INFO, OR, SE, NCAS, NCON)

    write_tsv(sumstats, "${daner.baseName}")
    """
}

process META {
    tag "${dataset.baseName}"

    cpus = 4
    memory = 16.GB
    time = '30m'

    input:
    tuple path(dataset), path(sumstats)

    output:
    tuple val(dataset.baseName), path("${dataset.baseName}.tsv")

    script:
    """
    #!python3
    import os
    
    os.environ["POLARS_MAX_THREADS"] = "${task.cpus}"
    import polars as pl

    sumstats_paths = "${sumstats}"

    sumstats = pl.concat(
        [pl.scan_csv(path, separator= "\\t", null_values = "NA", dtypes = {"INFO": pl.Float64}) for path in sumstats_paths.split()],
        how = "vertical"
    )

    meta = (
        sumstats.select(
            pl.col("*"),
            (pl.col("NCAS") + pl.col("NCON")).alias("N"),
            pl.col("OR").log().alias("BETA"),
            (1 / (pl.col("SE") ** 2)).alias("ivw"),
            (4 * pl.col("NCAS") * pl.col("NCON") / (pl.col("NCAS") + pl.col("NCON"))).alias("NEFF")
        )
        .select(
            pl.col("*"),
            (pl.col("BETA") * pl.col("ivw")).alias("wBETA"),
            (pl.col("INFO") * pl.col("N")).alias("wINFO"),
            (pl.col("FRQ_A") * pl.col("NCAS")).alias("wAFCAS"),
            (pl.col("FRQ_U") * pl.col("NCON")).alias("wAFCON")
        )
        .group_by("CHR", "BP", "SNP", "A1", "A2")
        .agg(
            (pl.sum("wBETA") /  pl.sum("ivw")).exp().alias("OR"),
            (1 / pl.sum("ivw")).sqrt().alias("SE"),
            (pl.sum("wINFO") / pl.sum("N")).alias("INFO"),
            (pl.sum("wAFCAS") / pl.sum("NCAS")).alias("AFCAS"),
            (pl.sum("wAFCON") / pl.sum("NCON")).alias("AFCON"),
            pl.sum("NCAS"),
            pl.sum("NCON"),
            pl.sum("NEFF")
        )
        .collect()
    )

    meta.write_csv("${dataset.baseName}.tsv", separator = "\\t")
    """
}