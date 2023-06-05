# Prepare LDSC and sumstats of well-powered symptoms across both samples


library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(GenomicSEM)

symptoms_h2 <- read_tsv(here::here("ldsc", "symptoms.h2.txt"))

# symptoms with positive heritabilities, sufficient chisq
symptoms_analyse <-
symptoms_h2 |>
    filter(cohorts %in% c("Clin", "Comm", "UKBt")) |>
    filter(h2 > 0, h2 < 1, Intercept > 1, MeanChiSq > 1.02)
    # "Comm.MDD1", "Comm.MDD4a", "Comm.MDD4b", "Comm.MDD6",
    # "Comm.MDD7", "Comm.MDD8", "Comm.MDD9", "UKBt.MDD1"


covstruct_prefix <- "comm.ukb.common.covstruct"
sumstats_prefix <- "comm.ukb.common.sumstats"
covstruct_r <- here::here("ldsc", paste(covstruct_prefix, "deparse.R", sep = "."))
covstruct_rds <- here::here("ldsc", paste(covstruct_prefix, "rds", sep = "."))
sumstats_rds <- here::here("sumstats", paste(sumstats_prefix, "rds", sep = "."))

# find the daner files
meta_daner_files <- list.files(here::here("meta/distribution/"),
    pattern = "^daner_[A-Za-z_]+\\.MDD[1-9][a-z]*_[A-Za-z]+\\.gz$",
    full.name = TRUE, recursive = TRUE)

# non-meta-analysed files
other_daner_files <- list.files("sumstats/UKB/Touchscreen",
    pattern = "^daner_[A-Za-z_]+\\.MDD[1-9][a-z]*_[A-Za-z]+\\.gz$",
    full.name = TRUE, recursive = TRUE)

daner_files <- c(meta_daner_files, other_daner_files)
# parse out info from filenames
daner_info <- str_match(basename(daner_files),
    "daner_([A-Za-z_]+)\\.(MDD[0-9a-b]+)")

symptoms_daner_analyse <- symptoms_analyse |>
left_join(tibble(daner=daner_files, cohorts = daner_info[,2], ref = daner_info[,3]),
          by = c("cohorts", "ref")) |>
select(Sample, sample_symptom, ref, samp_prev, pop_prev, filename, daner)

 write_tsv(symptoms_daner_analyse, here::here("ldsc", paste(covstruct_prefix, "prevs", "txt", sep = ".")))
	
if (!file.exists(covstruct_r)) {

  symptoms_covstruct <- ldsc(traits = symptoms_daner_analyse$filename,
    sample.prev = rep(0.5, times = length(symptoms_daner_analyse$samp_prev)),
    population.prev = symptoms_daner_analyse$pop_prev,
    ld = here::here("sumstats/reference/eur_w_ld_chr/"),
    wld = here::here("sumstats/reference/eur_w_ld_chr/"),
    trait.names = symptoms_daner_analyse$sample_symptom)

  dput(symptoms_covstruct, covstruct_r, control = c("exact"))
  saveRDS(symptoms_covstruct, covstruct_rds)
  # check for exact match of deparsed object
  identical(dget(covstruct_r), symptoms_covstruct)

} else {

  symptoms_covstruct <- dget(covstruct_r)

}

# download 1KG reference file from https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v

if (!file.exists(sumstats_rds)) {
    symptoms_sumstats <- 
    sumstats(files = symptoms_daner_analyse$daner,
    ref = here::here("sumstats", "reference", "reference.1000G.maf.0.005.txt"),
    trait.names = symptoms_daner_analyse$sample_symptom,
    se.logit = rep(TRUE, nrow(symptoms_daner_analyse)),
    OLS = NULL,
    linprob = rep(FALSE, nrow(symptoms_daner_analyse)),
    info.filter = 0.6,
    maf.filter = 0.01,
    keep.indel = FALSE,
    parallel = TRUE,
    cores = 2)
    saveRDS(symptoms_sumstats, sumstats_rds)
}

# test for traits that don"t require large modifications to the matrix
# symptoms_covstruct <- dget(here::here("ldsc/clin.comm.covstruct.deparse.R"))
# cc.model <- "
# MDD =~ NA*Comm.MDD1 + Comm.MDD4a + Comm.MDD4b + Comm.MDD6 + Comm.MDD7 + Comm.MDD8 + Comm.MDD9 + UKBt.MDD1
# MDD ~~ 1*MDD
# "
#  [1] Clin.MDD5b
#  [2] Clin.MDD7     
#  [3] Comm.MDD1       
#  [4] Comm.MDD4a     
#  [5] Comm.MDD4b    
#  [6] Comm.MDD5b
#  [7] Comm.MDD6        
#  [8] Comm.MDD7       
#  [9] Comm.MDD8   
# [10] Comm.MDD9          
# [11] UKBt.MDD1
