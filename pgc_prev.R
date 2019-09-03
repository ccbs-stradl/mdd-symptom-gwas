# calculate symptom prevalences in PGC cohorts

library(readxl)
library(tidyr)
library(readr)
library(dplyr)
library(yaml)


# load configuration file (see config-example.yaml')
config <- yaml.load_file('config.yaml')

# input and output directories and files
mdd_v1_dir <- config$data$pgc$mdd$v1
output_dir <- config$data$pgc$mdd$output
ukb_remove <- config$data$pgc$mdd$ukb_remove

symptoms <- c("MDD1", "MDD2", "MDD3a", "MDD3b", "MDD4a", "MDD4b", "MDD5a", "MDD5b", "MDD6", "MDD7", "MDD8", "MDD9")

secondary_phenotypes_files <- list.files(file.path(mdd_v1_dir, 'secondary_phenotypes', '1216_update'), pattern='xls', full.names=TRUE)

secondary_phenotypes_cohorts <-
  lapply(secondary_phenotypes_files,
         function(file)
          read_excel(file) %>%
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
               mutate(Study=if_else(is.na(Study),
                                    true=Sub.study,
                                    false=Study),
                      # harmonise 0/1 and 1/2 coding of case/control status
                      DSMIVMDD=if_else(rep(all(DSMIVMDD %in% c(0, 1)), n()),
                                       true=DSMIVMDD,
                                       false=DSMIVMDD-1)))

secondary_phenotypes <- bind_rows(secondary_phenotypes_cohorts)


status_symptoms <- 
secondary_phenotypes %>%
  select(DSMIVMDD, starts_with('MDD'), -MDD0, -MDDscore) %>%
  gather(key='symptom', value='value', MDD1:MDD9) %>%
  mutate(value=if_else(value %in% 0:1,
                       true=value, false=NA_real_))

sample_count <- 
status_symptoms %>%
filter(!is.na(value)) %>%
group_by(DSMIVMDD, symptom) %>%
tally()

k <- 0.15
sample_signif <- nchar(min(sample_count$n))
  
prevalences <- 
status_symptoms %>%
  group_by(DSMIVMDD, symptom) %>%
  summarise(prev=mean(value, na.rm=T)) %>%
  ungroup() %>%
  mutate(MDD=recode(DSMIVMDD, `0`='control_prev', `1`='case_prev')) %>%
  select(-DSMIVMDD) %>%
  spread(key=MDD, value=prev) %>%
  mutate(pop_prev=k*case_prev + (1-k)*control_prev) %>%
  mutate(case_prev=round(case_prev, digits=sample_signif),
         control_prev=round(control_prev, digits=sample_signif),
         pop_prev=round(pop_prev, digits=sample_signif))

write_tsv(prevalences, 'ldsc/pgc_dsm_prev.txt')

