# Prepare LDSC and sumstats of well-powered symptoms across from all meta-analysis

library(GenomicSEM)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)

dsm_mdd_symptoms_labels <-
read_delim("
MDD1;Mood;Mood;Dep
MDD2;Interest;Interest;Anh
MDD3;Weight⇅;Weight⇆;App
MDD3a;Weight⇊;Weight⇇;AppDec
MDD3b;Weight⇈;Weight⇉;AppInc
MDD4;Sleep⇅;Sleep⇆;Sle
MDD4a;Sleep⇊;Sleep⇇;SleDec
MDD4b;Sleep⇈;Sleep⇉;SleInc
MDD5;Motor⇅;Motor⇆;Moto
MDD5a;Motor⇈;Motor⇉;MotoInc
MDD5b;Motor⇊;Motor⇇;MotoDec
MDD6;Fatigue;Fatigue;Fatig
MDD7;Guilt;Guilt;Guilt
MDD8;Concentrate;Concentrate;Conc
MDD9;Suicidality;Suicidality;Sui
", col_names=c('ref', 'h', 'v', 'abbv'), delim=';')

symptoms_sample_prev_file <- here::here('meta', 'symptoms_prev.txt')
symptoms_sample_prev <- read_tsv(symptoms_sample_prev_file)

pop_prevs_w <-
symptoms_sample_prev %>%
mutate(w=case_when(cohorts == 'AGDS_PGC' ~ 0.15,
				   symptom %in% c('MDD1', 'MDD2') ~ 1.0,
				   TRUE ~ 0.57)) %>%
mutate(pop_prev=samp_prev*w) %>%
select(symptom, cohorts, pop_prev) %>%
pivot_wider(names_from=cohorts, values_from=pop_prev) %>%
group_by(symptom) %>%
mutate(pop_prev=mean(c(AGDS_PGC, ALSPAC_UKB))) %>%
select(symptom, pop_prev)

covstruct_prefix <- 'all.common.covstruct'
sumstats_prefix <- 'all.common.sumstats'
covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- here::here('ldsc', paste(covstruct_prefix, 'rds', sep='.'))
sumstats_rds <- here::here('sumstats', paste(sumstats_prefix, 'rds', sep='.'))

# list sumstats distribution directories
sumstats_files <- list.files(here::here('meta', 'munged'), 'ALL.+gz$', full.names=TRUE)

# pull out which cohorts and symptom 'x' this is from the filename (COHORTS_MDDx_*)
cohorts_symptoms <- str_match(basename(sumstats_files), '([A-Z_]+).(MDD[:digit:](a|b)?)')[,1]

sumstats_paths <- data.frame(filename=sumstats_files, sumstats=str_remove(basename(sumstats_files), '.sumstats.gz'))

sumstats_prevs <- 
sumstats_paths %>%
mutate(ref=str_match(sumstats, "(MDD[1-9a-b]+)")[,2]) %>%
left_join(pop_prevs_w, by=c('ref'='symptom')) %>%
left_join(dsm_mdd_symptoms_labels, by='ref') %>%
mutate(daner=file.path('meta', 'distribution', sumstats, paste0('daner_', sumstats, '.gz')))

sumstats_keep <- sumstats_prevs %>%
filter(ref %in% c('MDD1', 'MDD2', 'MDD3b', 'MDD4b', 'MDD6', 'MDD7', 'MDD8'))

write_tsv(sumstats_prevs, here::here('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
	
if(!file.exists(covstruct_r)) {

  symptoms_covstruct <- ldsc(traits=sumstats_keep$filename,
							 sample.prev=rep(0.5, length=sumstats_keep$pop_prev),
							 population.prev=sumstats_keep$pop_prev,
							 ld=here::here('sumstats/reference/eur_w_ld_chr/'),
							 wld=here::here('sumstats/reference/eur_w_ld_chr/'),
							 trait.names=sumstats_keep$trait_name)

  dput(symptoms_covstruct, covstruct_r, control=c('exact'))
  saveRDS(symptoms_covstruct, covstruct_rds)
  
  # check for exact match of deparsed object
  identical(dget(covstruct_r), symptoms_covstruct)

} else {

  symptoms_covstruct <- dget(covstruct_r)

}

# download 1KG reference file from https://utexas.box.com/s/vkd36n197m8klbaio3yzoxsee6sxo11v

if(!file.exists(sumstats_rds)){
	
	symptoms_sumstats <- sumstats(files=here::here(sumstats_keep$daner),
	                              ref=here::here('sumstats', 'reference', 'reference.1000G.maf.0.005.txt'),
								  trait.names=sumstats_keep$trait_name,
								  se.logit=rep(TRUE, nrow(sumstats_keep)),
								  OLS=NULL,
								  linprob=rep(FALSE, nrow(sumstats_keep)),
								  info.filter=0.6,
								  maf.filter=0.01,
								  keep.indel=FALSE,
								  parallel=TRUE,
								  cores=4)
						
	saveRDS(symptoms_sumstats, sumstats_rds)
}