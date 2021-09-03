# Format sumstats for MTAG

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

# open and format each set of sumstats
for(i in seq.int(nrow(sumstats_names))) {
    cohorts <- sumstats_names$cohorts[i]
    ref <- sumstats_names$ref[i]
    sample_symptom <- sumstats_names$SampleSymptom[i]
    
    cohorts_ref <- paste(cohorts, ref, sep='.')
    prefix <- list.files(here::here('meta', 'distribution'), cohorts_ref)
    daner_gz <- paste0('daner_', prefix, '.gz')
    sumstats_txt <- paste0(sample_symptom, '.txt')
    basic_xls <- list.files(here::here('meta', 'distribution', prefix), 'basic.+\\.xls', full.names=TRUE)
    
    if(length(basic_xls) > 0) {
        basic <- read_excel(basic_xls)
        
        lambda_gc <- basic %>% filter(Dataset == 'SUM') %>% pull(`LAMBDA-GC`)
        
        if(lambda_gc > 1.02) {
            cat(paste0(sumstats_txt, ' GC=', lambda_gc, '\n'))
        }
    }
    
    if(!file.exists(sumstats_txt)) {
        cat(paste("Reading", daner_gz, "\n"))
        daner <- read_tsv(here::here('meta', 'distribution', prefix, daner_gz), col_types=cols(SNP=col_character()))
        
        daner_patch <- daner#[-problems(daner)$row,]
        
        mtag_sumstats <- daner_patch %>%
            select(snpid=SNP, chr=CHR, bpos=BP, a1=A1, a2=A2, freq=starts_with('FRQ_A'), or=OR, se=SE, pval=P, Nca, Nco) %>%
            transmute(snpid, chr, bpos, a1, a2, freq, beta=log(or), se, pval, n=4*Nca*Nco / (Nca + Nco))
            
        
        cat(paste("Writing", sumstats_txt, "\n"))
        write_tsv(mtag_sumstats, sumstats_txt)
    }
    
}

