# calculate symptom prevalences in PGC cohorts

library(readxl)
library(tidyr)
library(readr)
library(dplyr)
library(yaml)
library(stringr)
library(psych)


# load configuration file (see config-example.yaml')
config <- yaml.load_file('config.yaml')

# input and output directories and files
mdd_v1_dir <- config$data$pgc$mdd$v1
output_dir <- config$data$pgc$mdd$output
ukb_remove <- config$data$pgc$mdd$ukb_remove

# numbering scheme for symptoms
symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")

# subdirectory of Excel files containing system data
secondary_phenotypes_files <- list.files(file.path(mdd_v1_dir, 'secondary_phenotypes', '1216_update'), pattern='xls', full.names=TRUE)

# coerce columns across files so they are of the same type
secondary_phenotypes_cohorts <-
  lapply(secondary_phenotypes_files,
         function(filename)
          read_excel(filename) %>%
            mutate_at(.funs=as.character,
                     .vars=c('Sub.study', 'Instrument', 'MDD0', 'ID2',
                             'SuicideDef', 'Bulimia', 'FH_MDD_1st',
                             'Study', 'Sex', 'NumSeriousSuicideAttempts',
                             'SuicideAttempt', 'P1WithinMDD',
                             'P2CrossDisorder', 'P3SubstanceUse',
                             'P4AnyPGC')) %>%
            mutate_at(.funs=as.numeric,
                      .vars=c('YrsEduc', 'AgeAtOnset', 'Early_Onset',
                              'NumEpisodes', 'Recurrent',
                              'DurationLongestEpisode',
                              'Hallucinations', 'Delusions', 'GenAnxDis',
                              'SocialPhobia', 'PanicDis', 'OCD',
                              'PostTraumaticDis', 'Anorexia',
                              'AlcoholAbuse', 'AlcoholDepend',
                              'MaxStandardDrinks', 'NicotineDepend',
                              'RegularSmoker', 'NumberCigarettes',
                              'CurrentSmoker', 'Cig100',
                              'CannabisAbuseDepend', 'Subs_dep',
                              'FH_BIP_1st', 'BMI', 'Waist_circumference',
                              'PostND', 'Migraine', 'Migraine_aura',
                              'PeriND', 'AgeLastMetCriteria', 'MDDscore',
                              'MDD1', 'MDD2', 'MDD3a', 'MDD3b', 'MDD4a',
                              'MDD4b', 'MDD5a', 'MDD5b', 'MDD6',
                              'MDD7', 'MDD8', 'MDD9', 'DimLibido',
                              'FH_Method', 'CommClin', 'Year_of_birth',
                              'Month_of_birth', 'AgeInterview')) %>%
               # get study name from filename
               mutate(STUDY=str_split(basename(filename), '_')[[1]][1],
                      # harmonise 0/1 and 1/2 coding of case/control status
                      DSMIVMDD=if_else(rep(all(DSMIVMDD %in% c(0, 1)), n()),
                                       true=DSMIVMDD,
                                       false=DSMIVMDD-1)))

# merge all study phenotype data together
secondary_phenotypes <- bind_rows(secondary_phenotypes_cohorts)


# extract symptom data and reshape to long format
cohort_status_symptoms <- 
secondary_phenotypes %>%
  select(STUDY, DSMIVMDD, starts_with('MDD'), -MDD0, -MDDscore) %>%
  gather(key='symptom', value='value', MDD1:MDD9) %>%
  mutate(value=if_else(value %in% 0:1,
                       true=value, false=NA_real_))

# count presence/absence of symptoms across MDD status
cohort_symptom_counts <- 
cohort_status_symptoms %>%
  group_by(STUDY, DSMIVMDD, symptom, value) %>%
  tally() %>%
  ungroup() %>%
  mutate(MDD=recode(DSMIVMDD, `1`='Case', `0`='Control'),
         status=recode(value, `1`='Present', `0`='Absent')) %>%
  select(Study=STUDY, MDD, Symptom=symptom, Status=status, N=n) 

# convert to wide for display
cohort_symptom_counts_wide <- 
cohort_symptom_counts %>%
  filter(!is.na(Status)) %>%
  spread(key=Status, value=n, fill='---') %>%
  unite(AbsentPresent, Absent, Present, sep=':') %>%
  mutate(Symptom=paste(mdd, symptom, sep=': ')) %>%
  select(STUDY, Symptom, AbsentPresent) %>%
  spread(key=Symptom, value=AbsentPresent)

# counts for symptom presence/absence within cases and controls
# exclude cohorts that don't have positive data within one category
symptom_counts <- 
cohort_symptom_counts %>%
group_by(Study, MDD) %>%
mutate(positive=rep('Present' %in% Status, n())) %>%
filter(positive, !is.na(Status))  %>%
group_by(MDD, Symptom, Status) %>%
summarize(N=sum(N))


## sample_count <- 
## status_symptoms %>%
## filter(!is.na(value)) %>%
## group_by(DSMIVMDD, symptom) %>%
## tally()
## 
## k <- 0.15
## sample_signif <- nchar(min(sample_count$n))
##   
## prevalences <- 
## status_symptoms %>%
##   group_by(DSMIVMDD, symptom) %>%
##   summarise(prev=mean(value, na.rm=T)) %>%
##   ungroup() %>%
##   mutate(MDD=recode(DSMIVMDD, `0`='control_prev', `1`='case_prev')) %>%
##   select(-DSMIVMDD) %>%
##   spread(key=MDD, value=prev) %>%
##   mutate(pop_prev=k*case_prev + (1-k)*control_prev) %>%
##   mutate(case_prev=round(case_prev, digits=sample_signif),
##          control_prev=round(control_prev, digits=sample_signif),
##          pop_prev=round(pop_prev, digits=sample_signif))
## 
## write_tsv(prevalences, 'ldsc/pgc_dsm_prev.txt')

dir.create('sumstats/PGC/CasesAllCohorts', recursive=TRUE, showWarnings=FALSE)

write_tsv(cohort_symptom_counts, 'sumstats/PGC/CasesAllCohorts/pgc_dsm_cohort_symptom_counts.txt')
write_tsv(symptom_counts, 'sumstats/PGC/CasesAllCohorts/pgc_dsm_symptom_counts.txt')

# Symptom correlations

# get cases and symptoms

cases_symptoms <- 
secondary_phenotypes %>%
filter(DSMIVMDD == 1) %>%
select(ID1, ID2, Study, starts_with('MDD'), -MDD0, -MDDscore) %>%
gather(key='symptom', value='value', MDD1:MDD9) %>%
filter(value %in% 0:1) %>%
group_by(Study) %>%
mutate(group_var=var(value)) %>%
ungroup() %>%
filter(group_var != 0) %>%
select(-group_var) %>%
spread(symptom, value) %>%
select(-ID1, -ID2, -Study)

cases_symptoms_cor <- tetrachoric(cases_symptoms)

dput(cases_symptoms_cor, 'sumstats/PGC/CasesAllCohorts/pgc_cases_tetra_cor.deprase.R', control=c('all', 'digits17'))
