# Prepare LDSC and sumstats of well-powered symptoms across both samples


library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(GenomicSEM)

symptoms_h2 <- read_tsv(here::here('ldsc', 'symptoms.h2.txt'))

# symptoms with positive heritabilities, sufficient chisq
symptoms_analyse <-
symptoms_h2 %>%
    filter(cohorts %in% c('ALL', 'UKBt')) %>%
    filter(h2 > 0) %>%
    filter(sample_symptom %in% c("AllAnh", "AllDep", "AllSui", "AllAppInc", "AllGuilt", "AllSleInc", "AllConc", "AllAppDec", "AllMotoInc"))

covstruct_prefix <- 'all.common.covstruct'
sumstats_prefix <- 'all.common.sumstats'
covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))
sumstats_rds <- here::here('sumstats', paste(sumstats_prefix, 'rds', sep='.'))

# find the daner files
meta_daner_files <- list.files("meta/distribution", pattern="^daner_[A-Z_]+\\.MDD[1-9][a-z]*_[A-Za-z]+\\.gz$", full.name=TRUE, recursive=TRUE)

# non-meta-analysed files
other_daner_files <- list.files("sumstats/UKB/Touchscreen", pattern="^daner_[A-Za-z_]+\\.MDD[1-9][a-z]*_[A-Za-z]+\\.gz$", full.name=TRUE, recursive=TRUE)

daner_files <- c(meta_daner_files, other_daner_files)
# parse out info from filenames
daner_info <- str_match(basename(daner_files), "daner_([A-Za-z_]+)\\.(MDD[0-9a-b]+)")

symptoms_daner_analyse <- symptoms_analyse %>%
left_join(tibble(daner=daner_files, cohorts=daner_info[,2], ref=daner_info[,3]),
          by=c('cohorts', 'ref')) %>%
select(Sample, sample_symptom, ref, samp_prev, pop_prev, filename, daner)

 write_tsv(symptoms_daner_analyse, here::here('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
	
if(!file.exists(covstruct_r)) {

  symptoms_covstruct <- ldsc(traits=symptoms_daner_analyse$filename,
							 sample.prev=rep(0.5, times=length(symptoms_daner_analyse$samp_prev)),
							 population.prev=symptoms_daner_analyse$pop_prev,
							 ld=here::here('sumstats/reference/eur_w_ld_chr/'),
							 wld=here::here('sumstats/reference/eur_w_ld_chr/'),
							 trait.names=symptoms_daner_analyse$sample_symptom)

  dput(symptoms_covstruct, covstruct_r, control=c('exact'))
  saveRDS(symptoms_covstruct, covstruct_rds)
  
  # check for exact match of deparsed object
  identical(dget(covstruct_r), symptoms_covstruct)

} else {

  symptoms_covstruct <- dget(covstruct_r)

}

# download 1KG reference file from https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v

if(!file.exists(sumstats_rds)){
	
	symptoms_sumstats <- sumstats(files=symptoms_daner_analyse$daner,
	                              ref=here::here('sumstats', 'reference', 'reference.1000G.maf.0.005.txt'),
								  trait.names=symptoms_daner_analyse$sample_symptom,
								  se.logit=rep(TRUE, nrow(symptoms_daner_analyse)),
								  OLS=NULL,
								  linprob=rep(FALSE, nrow(symptoms_daner_analyse)),
								  info.filter=0.6,
								  maf.filter=0.01,
								  keep.indel=FALSE,
								  parallel=TRUE,
								  cores=8)
						
	saveRDS(symptoms_sumstats, sumstats_rds)
}

# test for traits that don't require large modifications to the matrix
# cc.model <- "
# MDD =~ NA*AllAnh + AllDep + AllSui + AllAppInc + AllGuilt + AllSleInc + AllConc + AllAppDec + AllMotoInc
# MDD ~~ 1*MDD
# "
# 
# x <- usermodel(symptoms_covstruct, estimation='DWLS', model=cc.model)

#  1 AllAnh        
#  2 AllDep        
#  3 AllSui        
#  4 AllAppInc     
#  5 UkbDep        
#  6 AllGuilt      
#  7 UkbAnh        
#  8 AllFatig      
#  9 AllSleInc     
# 10 AllConc       
# 11 AllAppDec     
# 12 AllSleDec     
# 13 AllMotoInc   