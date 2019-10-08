---
title: UKB MDD symptom phenotypes 
author: Mark Adams 
output:
  html_document:
    toc: TRUE
    code_folding: hide
    number_sections: TRUE
    df_print: kable
    keep_md: true
  md_document:
    variant: markdown_github
---

# Data extraction

Extract data from the UKB release files


```bash

cd sumstats/UKB
ukb_unpack ukb12345.enc
wget  -nd  biobank.ndph.ox.ac.uk/showcase/util/encoding.ukb
ukb_conv ukb12345.enc_ukb r -oukb_mdd -iukb_mdd.field 

```


```r
source('sumstats/UKB/ukb_mdd.r')
```


```r
library(dplyr)
library(readr)
library(psych)
```

# Fields

## DSM lifetime symptoms (CIDI)



```r
mhq_cidi_fields <-
read_delim("
MDD1;f.20446.0.0 ;Ever had prolonged feelings of sadness or depression
MDD2;f.20441.0.0 ;Ever had prolonged loss of interest in normal activities
MDD3;f.20536.0.0 ;Weight change during worst episode of depression
MDD4;f.20532.0.0 ;Did your sleep change
MDD4a;f.20533.0.0;Trouble falling asleep
MDD4a;f.20535.0.0;Waking too early
MDD4b;f.20534.0.0;Sleeping too much
MDD6;f.20449.0.0 ;Feelings of tiredness during worst episode of depression
MDD7;f.20450.0.0 ;Feelings of worthlessness during worst period of depression
MDD8;f.20435.0.0 ;Difficulty concentrating during worst depression
MDD9;f.20437.0.0 ;Thoughts of death during worst depression
", col_names=c('Reference', 'Field', 'Description'), delim=';')

mhq_cidi_fields
```

<div class="kable-table">

Reference   Field         Description                                                 
----------  ------------  ------------------------------------------------------------
MDD1        f.20446.0.0   Ever had prolonged feelings of sadness or depression        
MDD2        f.20441.0.0   Ever had prolonged loss of interest in normal activities    
MDD3        f.20536.0.0   Weight change during worst episode of depression            
MDD4        f.20532.0.0   Did your sleep change                                       
MDD4a       f.20533.0.0   Trouble falling asleep                                      
MDD4a       f.20535.0.0   Waking too early                                            
MDD4b       f.20534.0.0   Sleeping too much                                           
MDD6        f.20449.0.0   Feelings of tiredness during worst episode of depression    
MDD7        f.20450.0.0   Feelings of worthlessness during worst period of depression 
MDD8        f.20435.0.0   Difficulty concentrating during worst depression            
MDD9        f.20437.0.0   Thoughts of death during worst depression                   

</div>


## Current symptoms (PHQ)


```r
mhq_phq_fields <- 
read_delim("
MDD1;f.20510.0.0;Recent feelings of depression
MDD2;f.20514.0.0;Recent lack of interest or pleasure in doing things
MDD3;f.20511.0.0;Recent poor appetite or overeating
MDD4;f.20517.0.0;Trouble falling or staying asleep or sleeping too much
MDD5;f.20518.0.0;Recent changes in speed amount of moving or speaking
MDD6;f.20519.0.0;Recent feelings of tiredness or low energy
MDD7;f.20507.0.0;Recent feelings of inadequacy
MDD8;f.20508.0.0;Recent trouble concentrating on things
MDD9;f.20513.0.0;Recent thoughts of suicide or self harm
", col_names=c('Reference', 'Field', 'Description'), delim=';')

mhq_phq_fields
```

<div class="kable-table">

Reference   Field         Description                                            
----------  ------------  -------------------------------------------------------
MDD1        f.20510.0.0   Recent feelings of depression                          
MDD2        f.20514.0.0   Recent lack of interest or pleasure in doing things    
MDD3        f.20511.0.0   Recent poor appetite or overeating                     
MDD4        f.20517.0.0   Trouble falling or staying asleep or sleeping too much 
MDD5        f.20518.0.0   Recent changes in speed amount of moving or speaking   
MDD6        f.20519.0.0   Recent feelings of tiredness or low energy             
MDD7        f.20507.0.0   Recent feelings of inadequacy                          
MDD8        f.20508.0.0   Recent trouble concentrating on things                 
MDD9        f.20513.0.0   Recent thoughts of suicide or self harm                

</div>
   
   
   
# Input coding

## CIDI
   

Code up different kinds of missingness (did not take MHQ vs declind screening symptoms vs answered negative to screening symptoms  vs declined symptom)


```r
# Symptoms coded with codings 502 and 503
# Coding 503 for Criterion A1/A2
## -818	Prefer not to answer
##    0	No
##    1	Yes
#
# 'Not assessed' = not assessed by MHQ
# 'Declined'     = preferred not to answer
# 'Absent' = reports not having symptom 
# 'Present'  = reports having symptom 

# Coding 502 for other Criteria (except weight change)
## -818	Prefer not to answer
## -121	Do not know
##    0	No
##    1	Yes

# 'Not assessed' = not assessed by MHQ
# 'Missing'      = did not respond to depression/anhedonia items
# 'Screened'     = answered 0 for depression and anhedonia
# 'Declined'     = did not know or preferred not to answer
# 'Absent' = reports not having symptom 
# 'Present'  = reports having symptom 

resp_prefer <- 'Prefer not to answer'
resp_know <- 'Do not know'

screen_symptom_502 <- function(symptom_field, depression_field, anhedonia_field) {
        case_when(symptom_field %in% c(resp_prefer, resp_know) ~ 'Declined',
                  symptom_field %in% 'No' ~ 'Absent',
                  symptom_field %in% 'Yes' ~ 'Present',
                  depression_field %in% 'No' | anhedonia_field %in% 'No' ~ 'Screened',
                  depression_field %in% resp_prefer & anhedonia_field %in% resp_prefer ~ 'Missing',
                  TRUE ~ 'Not assessed')
}

# Coding for sleep criteria

### f.20534.0.0		Sleeping.too.much									###
### f.20533.0.0		Trouble.falling.asleep									###
### f.20535.0.0		Waking.too.early									###
# Coding 508
## 0	No
## 1	Yes

# Coding 507 for weight change criterion
# 
## -818	Prefer not to answer
## -121	Do not know
##    0	Stayed about the same or was on a diet
##    1	Gained weight
##    2	Lost weight
##    3	Both gained and lost some weight during the episode

# 'Not assessed' = not assessed by MHQ
# 'Missing'      = did not respond to depression/anhedonia items
# 'Screened'     = answered 0 for depression and anhedonia
# 'Declined'     = did not know or preferred not to answer
# 'Absent'       = reports not having symptom 
# 'Present'      = reports having symptom 

resp_same <- 'Stayed about the same or was on a diet'
resp_gain <- 'Gained weight'
resp_lost <- 'Lost weight'
resp_both <- 'Both gained and lost some weight during the episode'


screen_symptom_507 <- function(symptom_field, depression_field, anhedonia_field) {
        case_when(symptom_field %in% c(resp_prefer, resp_know) ~ 'Declined',
                  symptom_field %in% 'No' ~ 'Absent',
                  symptom_field %in% c(resp_gain, resp_lost, resp_both) ~ 'Present',
                  depression_field %in% 'No' | anhedonia_field %in% 'No' ~ 'Screened',
                  depression_field %in% resp_prefer & anhedonia_field %in% resp_prefer ~ 'Missing',
                  TRUE ~ 'Not assessed')
}

# code weight 'gain'/'loss' as a symptom
screen_weight_507 <- function(symptom_field, present_levels, absent_levels, depression_field, anhedonia_field) {
        case_when(symptom_field %in% c(resp_prefer, resp_know) ~ 'Declined',
                  symptom_field %in% present_levels ~ 'Absent',
                  symptom_field %in% absent_levels ~ 'Present',
                  depression_field %in% 'No' | anhedonia_field %in% 'No' ~ 'Screened',
                  depression_field %in% resp_prefer & anhedonia_field %in% resp_prefer ~ 'Missing',
                  TRUE ~ 'Not assessed')
}
```


```r
CIDI_symptoms <- bd %>%
select(f.eid, f.20446.0.0, f.20441.0.0, f.20532.0.0, f.20534.0.0, f.20533.0.0, f.20517.0.0, f.20535.0.0, f.20435.0.0, f.20449.0.0, f.20450.0.0, f.20437.0.0, f.20536.0.0) %>%
mutate(
# Depression and anhedonia act as screening items for the rest of the symptoms
       cidi1=case_when(is.na(f.20446.0.0) ~ 'Not assessed',
                                 f.20446.0.0 == resp_prefer ~ 'Missing',
                                 f.20446.0.0 == 'No' ~ 'Absent',
                                 f.20446.0.0 == 'Yes' ~ 'Present'),
       cidi2=case_when(is.na(f.20441.0.0) ~ 'Not assessed',
                                 f.20441.0.0 == resp_prefer ~ 'Missing',
                                 f.20441.0.0 == 'No' ~ 'Absent',
                                 f.20441.0.0 == 'Yes' ~ 'Present'),
# Weight change
       cidi3=screen_weight_507(symptom_field=f.20536.0.0,
                               present_levels=c(resp_gain, resp_lost, resp_both),
                               absent_levels=c(resp_same),
                               depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0),
# Weight gain
       cidi3b=screen_weight_507(symptom_field=f.20536.0.0,
                                    present_levels=c(resp_gain, resp_both),
                                    absent_levels=c(resp_same, resp_lost),
                                    depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0),
# Weight loss
       cidi3a=screen_weight_507(symptom_field=f.20536.0.0,
                                    present_levels=c(resp_lost, resp_both),
                                    absent_levels=c(resp_same, resp_gain),
                                    depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0),
# Sleep (insomnia or hypersomnia)
       cidi4=screen_symptom_502(symptom_field=f.20532.0.0,
                                     depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0),
### f.20534.0.0		Sleeping.too.much					 
### f.20533.0.0		Trouble.falling.asleep		 
### f.20535.0.0		Waking.too.early					 
# Coding 508
## 0	No
## 1	Yes
# Insomnia
       cidi4a=case_when(f.20532.0.0 %in% c(resp_prefer, resp_know) ~ 'Declined',
                               f.20533.0.0 %in% 'Yes' | f.20535.0.0 %in% 'Yes' ~ 'Present',
                               f.20534.0.0 %in% 'Yes' ~ 'Absent',
                               f.20532.0.0 %in% 'No' ~ 'Absent',
                               f.20446.0.0 %in% 'No' | f.20441.0.0 %in% 'No' ~ 'Screened',
                               f.20446.0.0 %in% resp_prefer & f.20441.0.0 %in% resp_prefer ~ 'Missing',
                               TRUE ~ 'Not assessed'),
# Hypersomnia
       cidi4b=case_when(f.20532.0.0 %in% c(resp_prefer, resp_know) ~ 'Declined',
                                  f.20534.0.0 %in% 'Yes' ~ 'Present',
                                  f.20533.0.0 %in% 'Yes' | f.20535.0.0 %in% 'Yes' ~ 'Absent',
                                  f.20532.0.0 %in% 'No' ~ 'Absent',
                                  f.20446.0.0 %in% 'No' | f.20441.0.0 %in% 'No' ~ 'Screened',
                                  f.20446.0.0 %in% resp_prefer & f.20441.0.0 %in% resp_prefer ~ 'Missing',
                                  TRUE ~ 'Not assessed'),
# Tiredness
       cidi6=screen_symptom_502(symptom_field=f.20449.0.0,
                                     depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0),
# Worthlessness
       cidi7=screen_symptom_502(symptom_field=f.20450.0.0,
                                     depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0),
# Concentration
       cidi8=screen_symptom_502(symptom_field=f.20435.0.0,
                                     depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0),
# Thoughts of death
       cidi9=screen_symptom_502(symptom_field=f.20437.0.0,
                                     depression_field=f.20446.0.0, anhedonia_field=f.20441.0.0))
```

Code to binary (absent = 0, present = 1)


```r
Ab0Pr1 <- function(symptom) case_when(symptom == 'Absent' ~ 0L,
                                    symptom == 'Present' ~ 1L,
                                    TRUE ~ NA_integer_)

CIDI_binary <- CIDI_symptoms %>%
mutate(CIDI1=Ab0Pr1(cidi1),
       CIDI2=Ab0Pr1(cidi2),
       CIDI3=Ab0Pr1(cidi3),
       CIDI3a=Ab0Pr1(cidi3a),
       CIDI3b=Ab0Pr1(cidi3b),
       CIDI4=Ab0Pr1(cidi4),
       CIDI4a=Ab0Pr1(cidi4a),
       CIDI4b=Ab0Pr1(cidi4b),
       CIDI6=Ab0Pr1(cidi6),
       CIDI7=Ab0Pr1(cidi7),
       CIDI8=Ab0Pr1(cidi8),
       CIDI9=Ab0Pr1(cidi9))
```

# Phenotypic correlations

Calculate tetrachoric correlations for binary data


```r
options('mc.cores'=4)
cidi_cor <- tetrachoric(CIDI_binary %>% select(starts_with('CIDI', ignore.case=FALSE)))
```

```
## Warning in cor.smooth(mat): Matrix was not positive definite, smoothing was
## done
```

```r
cidi_cor
```

```
## Call: tetrachoric(x = CIDI_binary %>% select(starts_with("CIDI", ignore.case = FALSE)))
## tetrachoric correlation 
##        CIDI1 CIDI2 CIDI3 CIDI3a CIDI3b CIDI4 CIDI4a CIDI4b CIDI6 CIDI7
## CIDI1   1.00                                                          
## CIDI2   0.90  1.00                                                    
## CIDI3  -0.02 -0.22  1.00                                              
## CIDI3a -0.03 -0.13  0.67  1.00                                        
## CIDI3b  0.01 -0.18  0.65 -0.13   1.00                                 
## CIDI4   0.04  0.25 -0.35 -0.28  -0.19   1.00                          
## CIDI4a  0.05  0.17 -0.31 -0.28  -0.13   0.76  1.00                    
## CIDI4b  0.01  0.25 -0.18 -0.01  -0.25   0.66  0.04   1.00             
## CIDI6   0.00  0.39 -0.33 -0.16  -0.33   0.47  0.35   0.46   1.00      
## CIDI7   0.05  0.33 -0.21 -0.06  -0.25   0.26  0.19   0.30   0.41  1.00
## CIDI8   0.02  0.39 -0.33 -0.24  -0.23   0.51  0.46   0.34   0.63  0.44
## CIDI9   0.05  0.17 -0.21 -0.16  -0.13   0.18  0.17   0.12   0.19  0.25
##       CIDI8 CIDI9
## CIDI8  1.00      
## CIDI9  0.23  1.00
## 
##  with tau of 
##  CIDI1  CIDI2  CIDI3 CIDI3a CIDI3b  CIDI4 CIDI4a CIDI4b  CIDI6  CIDI7 
## -0.120  0.265  0.247 -0.195 -0.708 -0.834 -0.678  0.967 -0.912 -0.024 
##  CIDI8  CIDI9 
## -0.795 -0.052
```
