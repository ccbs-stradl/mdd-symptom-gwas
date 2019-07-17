# Symptom phenotypes in cases-only and cases-controls samples in the PGC MDD cohorts

library(yaml)
library(gdata)
library(stringr)
library(dplyr)
library(tidyr)

# load configuration file (see config-example.yaml')
config <- yaml.load_file('config.yaml')

# input and output directories and files
mdd_v1_dir <- config$data$pgc$mdd$v1
output_dir <- config$data$pgc$mdd$output
ukb_remove <- config$data$pgc$mdd$ukb_remove


##################################################################################################
#Part 1

#Extract 9 symptom phenotypes
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")
studies <- c("RADIANT", "STR", "ROTTERDAM", "QIMR", "NESDANTR", "GSK", "GENRED", "GENRED2", "EDINBURGH", "CoLaus", "COFAMS", "BonnMH")
study_noupdate <- c("STARD", "SHIP0", "MARS", "Janssen", "GENPOD")

dir.create(output_dir, recursive=TRUE)

# find spreadsheets with symptom-level data
updated_sheets <- list.files(path=file.path(mdd_v1_dir, 'secondary_phenotypes/1216_update'), pattern='_PGC2updated.xls', full.names=T)
no_update_sheets <- list.files(path=file.path(mdd_v1_dir, 'secondary_phenotypes/1216_update'), pattern='_PGC2.xls', full.names=T)
sheets <- c(updated_sheets, no_update_sheets)

# extract cohort studies names from filesname
xls_filesnames <- sapply(sheets, basename)

studies <- sapply(str_split(xls_filesnames, pattern='_', n=2), first)
names(sheets) <- studies


# open each spreadsheet
mdd_list <- lapply(sheets, read.xls, stringsAsFactors=FALSE)

# coerce column types so that they are mergable
# convert columns that are reported with a mix of integer and character responses to character
# if Study name is missing, copy over sub-study name
mdd_list_chars <- lapply(mdd_list, function(df) {
        df %>% mutate_at(.funs=as.character,
                         .vars=c('Sub.study', 'Instrument', 'MDD0', 'ID2', 'SuicideDef', 'Bulimia',
                           'FH_MDD_1st', 'Study', 'Sex', 'NumSeriousSuicideAttempts', 'SuicideAttempt',
                           'P1WithinMDD', 'P2CrossDisorder', 'P3SubstanceUse', 'P4AnyPGC')) %>%
               mutate(Study=if_else(is.na(Study), true=Sub.study, false=Study))
        })

# bind together
mdd <- bind_rows(mdd_list_chars) %>% as_tibble() 

# check coding of DSMIVMDD categories in each study
mdd %>% group_by(Study, DSMIVMDD) %>% tally() %>% as.data.frame()

# code symptom phenotypes
# QC symptoms: replace values that are not 0/1 with missing
#     mutate_at() accepts "tidy dot" formula syntax for the function to apply, i.e.,
#         .funs=~some expression(.)
#     where "." is being replaced by each column names listed in .vars.
#     this saves us having type "some expression" 12 times, like
#         mutate(MDD1=some expression(MDD1), MDD2=some expression(MDD2), ...
# 
mdd_symptoms <- mdd %>%
  select(ID1, ID2, Study, Sub.study, DSMIVMDD, starts_with('MDD'), -MDDscore, -MDD0) %>%
  mutate_at(.vars=c('MDD1',  'MDD2',  'MDD3a', 'MDD3b', 'MDD4a', 'MDD4b',
                    'MDD5a', 'MDD5b', 'MDD6',  'MDD7',  'MDD8', 'MDD9'),
            .funs=~if_else(. %in% c(0, 1), true=., false=NA_integer_))
  

#Part 2a
##Subset by studies that have symptom data for cases and controls or cases-only

# determine whether each case/control group has symptom data
# check that there is variation (e.g., controls have data, not just all 0s
mdd_symptoms %>%
  gather(key='symptom', value='value', MDD1:MDD9) %>%
  filter(!is.na(value)) %>%
  group_by(Study, DSMIVMDD, symptom) %>%
  summarize(varies=0 %in% value & 1 %in% value) %>%
  group_by(Study, symptom) %>%
  tally() %>%
  ungroup() %>%
  transmute(Study, symptom, grouping=case_when(n == 1 ~ 'Ca',
                                               n == 2 ~ 'CaCo',
                                               TRUE   ~ 'None')) %>%
  spread(symptom, grouping)

# TODO: remove overlap with UKB (see below for necessary ID1 munging

# TODO: save phenotypes out to *_CasesControls_PGCMDD2_recoded.tsv and *_CasesOnly_PGCMDD2_final.tsv
# columns: ID1, ID2, PHENO
# coding missing as -9


#Part 3
#Remove overlapping samples with UKB - remove according to ID1 and ID2 

UKB <- read.table(ukb_remove, sep="\t", header=F)
UKB$ID1adj <- gsub("\\*.*","",UKB$V2) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
for (i in symptoms) {
MDD <-read.table(file.path(output_dir, paste0(i,"_CasesControls_PGCMDD2_recoded.tsv", sep="")), sep="\t", header=T)
MDD$ID1adj <- gsub("\\*.*","",MDD$ID1) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
print(dim(MDD))
#filter MDD for IDs not in the UKB sample
MDD2 <- MDD %>%
filter(!MDD$ID2 %in% UKB$V3 & !MDD$ID1adj %in% UKB$V2)
print(dim(MDD2))
# check how many overlapped that had phenotypes
MDD %>%
  filter(MDD$ID2 %in% UKB$V3 & MDD$ID1adj %in% UKB$ID1adj) %>%
  filter_at(i, any_vars(. != -9)) %>%
  group_by(Study) %>%
  tally()
#293 removed due to overlap with UKB - 205 from RADIANT which makes sense
MDD2 <- MDD2[,-c(3,4,6)]
write.table(MDD2, file.path(output_dir, paste0(i,"_CasesControls_PGCMDD2_final.tsv", sep="")), sep="\t", row.names=F, quote=F)
}

for (i in symptoms) {
MDD <-read.table(file.path(output_dir, paste0(i,"_CasesOnly_PGCMDD2_recoded.tsv", sep="")), sep="\t", header=T)
MDD$ID1adj <- gsub("\\*.*","",MDD$ID1) #Recode ID1 to remove IID component as there may be inconsistency between two datasets
print(dim(MDD))
#filter MDD for IDs not in the UKB sample
MDD2 <- MDD %>%
filter(!MDD$ID2 %in% UKB$V3 & !MDD$ID1adj %in% UKB$V2)
print(dim(MDD2))
MDD %>%
  filter(MDD$ID2 %in% UKB$V3 & MDD$ID1adj %in% UKB$ID1adj) %>%
  filter_at(i, any_vars(. != -9)) %>%
  group_by(Study) %>%
  tally()
#293 removed due to overlap with UKB - 205 from RADIANT which makes sense
MDD2 <- MDD2[,-c(3,4,6)]
write.table(MDD2, file.path(output_dir, paste0(i,"_CasesOnly_PGCMDD2_final.tsv", sep="")), sep="\t", row.names=F, quote=F)
}

#COMPLETE

