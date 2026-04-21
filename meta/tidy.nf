/* tidy sumstats */
import groovy.json.JsonOutput

params.sumstats = "sumstats/cohorts/**/daner*.gz"
params.dbsnp = "sumstats/reference/dbSNP155"
params.cohorts = "meta/cohort_alignment.txt"

workflow {

  // sumstats files keyed to their filenames
  SUMSTATS_CH = Channel.fromPath(params.sumstats)
    .map { it -> [it.name, it]}

  // reference directory
  DBSNP_CH = Channel.fromPath(params.dbsnp, type: "dir")

  // list of datasets with metadata, keyed to filename
  COHORTS_CH = Channel.fromPath(params.cohorts)
    .splitCsv(header: true, sep: '\t')
    .map { it -> [it.filename, it]}

  // merge sumstats with metadata
  SUMSTATS_INFO_CH = COHORTS_CH
    .join(SUMSTATS_CH)
  
  // clean and preformat
  FORMAT_CH = FORMAT(SUMSTATS_INFO_CH) 
    .combine(DBSNP_CH)

  TIDY_CH = TIDY(FORMAT_CH)
  ALIGN_CH = ALIGN(TIDY_CH)

}

// Prepare sumstats for tidying
process FORMAT {
  tag "${filename}"

  memory 8.GB

  input:
  tuple val(filename), val(info), path(daner)

  output:
  tuple val("${info.study}-${info.group}-${info.reference}"), val(info), path("${info.study}-${info.group}-${info.reference}.gz")

  script:
  """
  #!Rscript
  library(dplyr)
  library(readr)
  library(stringr)

  daner <- read_table("${daner}", na = c("NA", "Inf", "nan"))

  # extract case/control counts from header
  frq_a_col <- str_subset(names(daner), "FRQ_A")
  frq_u_col <- str_subset(names(daner), "FRQ_U")
  Nca <- as.numeric(str_extract(frq_a_col, "[[:digit:]]+"))
  Nco <- as.numeric(str_extract(frq_u_col, "[[:digit:]]+"))

  # check for NCAS/NCON columns
  if(!"NCAS" %in% names(daner)) {
      daner_n <- daner |>
        mutate(NCAS = Nca, NCON = Nco)
  } else {
      daner_n <- daner
  }

  # format to tidyGWAS expectations
  pretidy <- daner_n |>
    na.omit() |>
    mutate(B = log(OR)) |>
    select(
      CHR = CHR,
      POS = BP,
      RSID = SNP,
      EffectAllele = A1,
      OtherAllele = A2,
      B = B,
      SE = SE,
      EAF = starts_with("FRQ_A"),
      P = P,
      CaseN = NCAS,
      ControlN = NCON,
      INFO
    )

  write_tsv(pretidy, "${info.study}-${info.group}-${info.reference}.gz")
  """
}

// Tidy sumstats
process TIDY {
  tag "${dataset}"

  cpus 8
  memory 8.GB
  beforeScript "OMP_NUM_THREADS=8"

  input:
  tuple val(dataset), val(info), path(daner), path(dbsnp)

  output:
  tuple val(dataset), val(info), path("${dataset}")

  script:
  """
  #!Rscript 
  library(tidyGWAS)

  cleaned <- tidyGWAS(
    "${daner}",
    dbsnp = "${dbsnp}",
    output_dir = "${dataset}"
  )
  """
}

// Align and output with metadata
process ALIGN {
  tag "${dataset}"

  publishDir "sumstats/tidy"

  cpus 4
  memory 8.GB
  beforeScript "OMP_NUM_THREADS=4"

  input:
  tuple val(dataset), val(info), path(sumstats)

  output:
  tuple val(dataset), val(info), path("${dataset}-ref.parquet"), path("${dataset}-ref.json")

  script:
  """
  #!Rscript 
  library(tidyGWAS)
  library(arrow)
  library(dplyr)

  sumstats <- open_dataset("${sumstats}/tidyGWAS_hivestyle")
  metadata_str <- '${JsonOutput.prettyPrint(JsonOutput.toJson(info))}'

  aligned <- align_to_ref(sumstats, "REF_37") |>
    group_by(CHR)

  write_dataset(aligned, "${dataset}-ref.parquet", format = "parquet", hive_style = TRUE)
  write(metadata_str, "${dataset}-ref.json")
  """
}