# Convert MTAG'd sumstats to daner


library(GenomicSEM)
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(readxl)

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
MDD5;Motor⇅;Motor⇆;Psyc
MDD5a;Motor⇈;Motor⇉;PsycInc
MDD5b;Motor⇊;Motor⇇;PsycDec
MDD6;Fatigue;Fatigue;Fatig
MDD7;Guilt;Guilt;Guilt
MDD8;Concentrate;Concentrate;Conc
MDD9;Suicidality;Suicidality;Sui
", col_names=c('ref', 'h', 'v', 'abbv'), delim=';')

# Check which symptoms have genetic variance > 0
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
sumstats_prefix <- 'agds_pgc.alspac_ukb.sumstats'
covstruct_r <- here::here('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

symptoms <- dimnames(symptoms_covstruct$S)[[2]]

symptoms_pos <- symptoms[which(diag(symptoms_covstruct$S) > 0)]
print(symptoms_pos)

# Clin/Pop naming scheme
cohort_symptom_names <- str_split(symptoms_pos, '\\.')
cohort_names <- sapply(cohort_symptom_names, first)
ref_names <- sapply(cohort_symptom_names, last)

sumstats_names <- 
tibble(cohorts=cohort_names, ref=ref_names) %>%
mutate(Sample=case_when(cohorts %in% 'AGDS_PGC' ~ 'Clin',
                        cohorts %in% 'ALSPAC_UKB' ~ 'Pop',
                        TRUE ~ NA_character_)) %>%
left_join(dsm_mdd_symptoms_labels, by='ref') %>%
mutate(SampleSymptom=paste(Sample, abbv, sep='_'))

#  1	Clin_AppInc.txt
#  2	Clin_SleDec.txt
#  3	Clin_SleInc.txt
#  4	Pop_Dep.txt
#  5	Pop_Anh.txt
#  6	Pop_AppDec.txt
#  7	Pop_AppInc.txt
#  8	Pop_SleDec.txt
#  9	Pop_SleInc.txt
# 10	Pop_Fatig.txt
# 11	Pop_Guilt.txt
# 12	Pop_Sui.txt

mtag_traits <- 
c("Clin_AppInc",
"Clin_SleDec",
"Clin_SleInc",
"Pop_Dep",
"Pop_Anh",
"Pop_AppDec",
"Pop_AppInc",
"Pop_SleDec",
"Pop_SleInc",
"Pop_Fatig",
"Pop_Guilt",
"Pop_Sui")

dir.create(here::here('mtag', 'daner'))

for(i in seq.int(length(mtag_traits))) {
    # open the MTAG sumstats
    trait <- mtag_traits[i]
    trait_txt <- here::here('mtag', 'clin_pop', paste0('clinpop_trait_', i, '.txt'))
    
    mtag_daner_gz <- here::here('mtag', 'daner', paste0('daner_', trait, '.mtag', '.gz'))
    
    if(!file.exists(mtag_daner_gz)) {
    
        cat(paste('Reading', basename(trait_txt), '\n'))
        mtag <- read_tsv(trait_txt, col_types=cols(SNP=col_character()))
        
        # check if trait has any genome-wide significant SNPs
        gw <- mtag %>% filter(mtag_pval <= 5e-8)
        
        if(nrow(gw) > 1) {
            # open the original daner file
            sumstats_info <- sumstats_names %>% filter(SampleSymptom == trait)
            cohorts = sumstats_info %>% pull(cohorts)
            ref <- sumstats_info %>% pull(ref)
            cohorts_ref <- paste(cohorts, ref, sep='.')
            prefix <- list.files(here::here('meta', 'distribution'), cohorts_ref)
            daner_gz <- paste0('daner_', prefix, '.gz')
            
            cat(paste('Reading', daner_gz, '\n'))
            daner <- read_tsv(here::here('meta', 'distribution', prefix, daner_gz), col_types=cols(SNP=col_character()))
            
            mtag_daner <- 
            mtag %>% select(SNP, A1, A2, starts_with('mtag')) %>%
            inner_join(daner, by='SNP', suffix=c('.mtag', '.daner')) %>%
            # get log odds to match original coding
            mutate(logOR=case_when(A1.mtag == A1.daner & A2.mtag == A2.daner ~ mtag_beta,
                                   A1.mtag == A2.daner & A2.mtag == A1.daner ~ -mtag_beta,
                                   TRUE ~ NA_real_)) %>%
            filter(!is.na(logOR)) %>%
            mutate(or=exp(logOR)) %>%
            select(CHR, SNP, BP, A1=A1.daner, A2=A2.daner, starts_with('FRQ'), INFO,
                      OR=or, SE=mtag_se, P=mtag_pval, ngt, Direction, HetISqt, HetDf, HetPVa, Nca, Nco, Neff_half) %>%
            arrange(CHR, BP)
            
            cat(paste('Writing', basename(mtag_daner_gz), '\n'))
            write_tsv(mtag_daner, mtag_daner_gz)
            
        }
    }
}
    
    