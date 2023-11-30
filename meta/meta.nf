/* LOO meta-analysis */

nextflow.enable.dsl = 2

params.sumstats = "aligned/*.gz"
params.datasets = "datasets/*.meta"

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
        .map { it -> [it[1], it[0]] }
    
    // List meta analyses that each dataset is specified for
    DATASETS_META_CH =
    DATASETS_CH
        .groupTuple()
    
    // Line up files listed in datasets with available daner files
    // Keep files listed in datasets, marking those where the file is not found
    DATASETS_DANERS_CH =
    DATASETS_META_CH
        .join(DANER_CH, remainder: true)
        .filter { it[1] != null }
    
    
    // Check whether any of the listed datasets do not have an associated daner file
    DATASETS_DANERS_CHECK_CH =
    DATASETS_DANERS_CH
        .branch { missing: it[2] == null
                 complete: true
                }
                
    // Print missing daner files
    DATASETS_DANERS_CHECK_CH
        .missing
        .transpose()
        .subscribe { println "Meta-analysis ${it[1].baseName} missing input dataset ${it[0]}" }
        
    // Unique daner files that are in the datasets
    DANER_SETS_CH =
    DATASETS_DANERS_CHECK_CH
        .complete
        .map { it -> [it[0], it[2]] }
    
    // Unzip daner files so they can be streamed
    // Harmonise allele frequency column names
    SUMSTATS_CH = SUMSTATS(DANER_SETS_CH)
    
    // Match sumstat files back up with meta-analysis dataset lists
    // Regroup sumstats into meta-analyses
    DATASETS_SUMSTATS_CH =
    DATASETS_DANERS_CHECK_CH
        .complete
        .map { it -> [it[0], it[1]] }
        .join(SUMSTATS_CH)
        .map { [it[2], it[1]] }
        .transpose()
        .map { [it[1], it[0]] }
        .groupTuple()
    
    // Perform meta analysis and post process
    META_CH = META(DATASETS_SUMSTATS_CH)
    POST_CH = POST(META_CH)
}

// Prepare danerfiles for input to meta-analysis
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

// Perform meta analysis
process META {
    tag "${meta.baseName}"

    cpus = 4
    memory = 16.GB
    time = '30m'

    input:
    tuple path(meta), path(sumstats)

    output:
    tuple val(meta.baseName), path("${meta.baseName}.tsv")

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

    meta.write_csv("${meta.baseName}.tsv", separator = "\\t")
    """
}

// Postprocess meta-analysis
process POST {
    tag "${meta}"
    
    publishDir "output", mode: "copy" 
    
    cpus = 1
    memory = 16.GB
    time = '30m'
    
    input:
    tuple val(meta), path(sumstats)
    
    output:
    tuple val(meta), path("${sumstats.baseName}.gz")
    
    script:
    """ 
    #!Rscript
    
    library(dplyr)
    library(readr)
    library(plyranges)
    
    sumstats <- read_tsv("${sumstats}")
    
    # filter based on 80% max(Neff)
    max_neff <- sumstats |> summarize(neff = max(NEFF)) |> pull(neff)
    
    sumstats_neff <- sumstats |>
        filter(NEFF >= 0.8 * max_neff)
        
    # remove duplicate positions
    sumstats_gr <- sumstats_neff |>
        transmute(seqnames = CHR, start = BP, width = 1, SNP, index = seq_along(SNP)) |>
        as_granges()
    
    duplicate_positions <- 
    join_overlap_self(sumstats_gr) |>
        filter(index != index.overlap) |>
        as_tibble()
        
    duplicate_snps <- c(
        pull(duplicate_positions, SNP),
        pull(duplicate_positions, SNP.overlap)    
    )
    
    sumstats_post <- sumstats_neff |> 
        filter(!SNP %in% duplicate_snps) |>
        mutate(P = pchisq(log(OR)^2 / SE^2, df = 1, lower.tail = FALSE))
        
    write_tsv(sumstats_post, "${sumstats.baseName}.gz")
    """
}