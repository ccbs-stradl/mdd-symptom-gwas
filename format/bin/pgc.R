# Create VCF-like PGC sumstats file with header

# Rscript bin/pgc.R --name=MDD1 --daner=../meta/distribution/ALSPAC_UKB.MDD1_depressed/daner_ALSPAC_UKB.MDD1_depressed.gz --basic=../meta/distribution/ALSPAC_UKB.MDD1_depressed/basic.ALSPAC_UKB.MDD1_depressed.num.xls --analyst="Mark James Adams, UoE" --fai=human_g1k_v37.fasta.fai --cff=../CITATION.CFF --template=pgc.glue --out=pgc-mdd1-comm --name="MDD symptom: low mood (MDD1)" --pop=EUR

library(dplyr)
library(readr)
library(stringr)
library(yaml)
library(lubridate)
library(readxl)
library(optparse)

cat("Formatting sumstats to pgc\n")

parser <- OptionParser()
parser <- add_option(parser, c("-n", "--name"), type = "character",
help = "Phenotype ID")
parser <- add_option(parser, c("-d", "--daner"), type = "character",
    help = "Daner sumstats file")
parser <- add_option(parser, c("-b", "--basic"), type = "character",
    help = "Summary spreadsheet ('basic') file")
parser <- add_option(parser, c("-t", "--template"), type = "character",
    help = "Header template")
parser <- add_option(parser, c("-a", "--analyst"), type = "character",
    help = "Name of analyst")
parser <- add_option(parser, c("-f", "--fai"), type = "character",
    help = "fasta.fai contig file")
parser <- add_option(parser, c("-c", "--cff"), type = "character",
    help = "Citation CFF file")
parser <- add_option(parser, c("-p", "--pop"), type = "character",
    help = "Reference population ID")
parser <- add_option(parser, c("-o", "--out"), type = "character",
    help = "Output filename")

ARGS <- commandArgs(TRUE)

args <- parse_args(parser, args = ARGS)

filedate <- format(now(), "%Y-%Om-%d")

# read analyst from config file
analyst <- args$analyst

# read citation information
cff <- read_yaml(args$cff)

# manuscript and DOI information

doi_idx <- which(sapply(cff$identifiers, function(i) i$type) == "doi")
doi <- cff$identifiers[[doi_idx]]$value
sumstats_url <- cff$url
code_url <- cff$`repository-code`

methods <- ""
acknowledgments <- "The PGC has received funding from the US National Institute of Mental Health (5 U01MH109528-04). Statistical analyses were carried out on the Genetic Cluster Computer (http://www.geneticcluster.org) hosted by SURFsara and financially supported by the Netherlands Scientific Organization (NWO 480-05-003) along with a supplement from the Dutch Brain Foundation and the VU University Amsterdam."
abstract <- cff$abstract

# Genome build contig
cat(str_glue("Reading contig {args$fai}\n"))
fasta_fai <-
    read_table(args$fai,
        col_names = c("NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"))

# open daner file
cat(str_glue("Reading daner {args$daner}\n"))
daner <- read_tsv(args$daner,
    col_types = cols(SNP = col_character()))

# Remove problematic rows
if (nrow(problems(daner)) > 0) {
    daner_qc <- daner |> slice(-problems(daner)$row)
} else {
    daner_qc <- daner
}

# format contig
daner_chr <- daner_qc |>
    select(CHR) |>
    distinct(CHR) |>
    mutate(CHR = as.character(CHR)) |>
    pull(CHR)

# map fai "X" -> "23"
contig_lengths <- fasta_fai |>
mutate(CHR = case_when(NAME == "X" ~ "23",
                      NAME == "XY" ~ "24",
                      NAME == "MT" ~ "25",
                      TRUE ~ NAME)) |>
filter(CHR %in% daner_chr) |>
mutate(contig = str_glue("##contig=<ID={CHR},length={LENGTH}>"))

contig <- paste(pull(contig_lengths, contig), collapse = "\n")

# get number of cases and controls from the header
frq_cols_split <- str_split(str_subset(names(daner), "FRQ"), "_")
ncase <- as.numeric(last(frq_cols_split[[1]]))
ncontrol <- as.numeric(last(frq_cols_split[[2]]))
ntrio <- 0

# open cohort files

# Analysed cohorts
basic <- read_excel(args$basic)

# cohort counts
cohort_counts <-
    basic |>
    filter(Dataset != "SUM")

analysis_datasets <- cohort_counts$Dataset
analysis_ncases <- cohort_counts$N_cases
analysis_ncontrols <- cohort_counts$N_controls
analysis_neffhalf <- cohort_counts$N_eff_half
analysis_nsnps <- cohort_counts$`N-SNPs`

name <- args$name 

cat(str_glue("Making sumstats file for: {name}\n"))

# dataset lists

dataset_list <- paste(analysis_datasets, collapse = ",")
neff <- 2 * sum(analysis_neffhalf)

# reference population
pop <- toupper(args$pop)

# number of variants
variants <- nrow(daner_qc)

# output prefix
id <- args$out

# header template
cat("Preparing header\n")

header_glues <- readLines(args$template)
headers <- sapply(header_glues, str_glue)
header <- paste(headers, collapse = "\n")

cat("#########################################################\n")
cat("#########################################################\n")
cat("#########################################################\n")
cat("\n\n\n")
cat(header)

cat("\n\n\n")
cat("#########################################################\n")
cat("#########################################################\n")
cat("#########################################################\n")

# format table
pgc_sumstats <- daner_qc |>
select(CHR, BP, SNP, A1, A2, OR, SE,
       FRQ_A = starts_with("FRQ_A"), FRQ_U = starts_with("FRQ_U"), P,
       INFO, ngt, Neff_half, Nca, Nco, Direction, HetISqt, HetDf, HetPVa) |>
transmute(`#CHROM`=CHR, POS=BP, ID=SNP, A1=A1, A2=A2,
          BETA = log(OR), SE, PVAL = P, NGT = ngt, FCAS = FRQ_A, FCON = FRQ_U,
          IMPINFO = INFO, NEFF = 2*Neff_half,
          NCAS = Nca, NCON = Nco, DIRE = Direction,
          HETI = HetISqt, HETDF = HetDf, HETPVAL = HetPVa)

out <- paste(id, 'tsv', sep='.')

cat(str_glue("Writing sumstats to {out}"))

cat(header, "\n", file = out)

write_tsv(pgc_sumstats, file = out, append = TRUE, col_names = TRUE)
