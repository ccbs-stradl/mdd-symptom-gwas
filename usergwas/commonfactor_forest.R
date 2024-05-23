library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

sumstats <- readRDS(here::here("sumstats", "clin.comm.ukb.common.sumstats.rds"))
gwas <- read_tsv(here::here("meta", "usergwas", "comm.ukb.common.sumstats.txt.gz"))
loci <- read_tsv(here::here("usergwas", "comm.ukb.common.loci.txt"))

loci_keep <- loci |>
  select(CHR, BP, SNP)

sumstats |>
  filter(SNP %in% pull(loci, SNP))