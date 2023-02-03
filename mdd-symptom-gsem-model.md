GenomicSEM of MDD symptoms
================
Mark Adams, Bradley Jermy, Jackson Thorp, Andrew Grotzinger, Michel
Nivard

# Setup

## R packages

R version

``` r
R.version
```

    ##                _                           
    ## platform       x86_64-pc-linux-gnu         
    ## arch           x86_64                      
    ## os             linux-gnu                   
    ## system         x86_64, linux-gnu           
    ## status                                     
    ## major          4                           
    ## minor          1.2                         
    ## year           2021                        
    ## month          11                          
    ## day            01                          
    ## svn rev        81115                       
    ## language       R                           
    ## version.string R version 4.1.2 (2021-11-01)
    ## nickname       Bird Hippie

Package installation

``` r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr', 'corrplot')
for(pack in required_packages) if(!require(pack, character.only=TRUE)) install.packages(pack)

if(!require(GenomicSEM)) remotes::install_github("MichelNivard/GenomicSEM")
```

GenomicSEM version

``` r
require(readr)
require(tidyr)
require(stringr)
require(dplyr)
require(ggplot2)
require(corrplot)
require(GenomicSEM)

packageVersion("GenomicSEM")
```

    ## [1] '0.0.5'

# Symptom labels

MDD DSM symptoms are numbered 1-9:

``` r
# plot labels

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
```

    ## Rows: 15 Columns: 4
    ## ── Column specification ──────────────────────────────────────────────────────────────────
    ## Delimiter: ";"
    ## chr (4): ref, h, v, abbv
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
dsm_mdd_symptoms_reference <-
read_delim("
MDD1;Depressed mood most of the day, nearly every day
MDD2;Markedly diminished interest or pleasure in all, or almost all, activities most of the day, nearly every day
MDD3;Significant change in weight or appetite
MDD3a;Significant weight loss or decrease in appetite
MDD3b;Significant weight gain or increase in appetite
MDD4;Sleeping too much or not sleeping enough
MDD4a;Insomnia nearly every day
MDD4b;Hypersomnia nearly every day
MDD5;Changes in speed/amount of moving or speaking
MDD5a;Psychomotor agitation nearly every day
MDD5b;Psychomotor slowing nearly every day
MDD6;Fatigue or loss of energy nearly every day
MDD7;Feelings of worthlessness or excessive or inappropriate guilt
MDD8;Diminished ability to think or concentrate, or indecisiveness
MDD9;Recurrent thoughts of death or suicide or a suicide attempt or a specific plan for attempting suicide
", col_names=c('Reference', 'Description'), delim=';')
```

    ## Rows: 15 Columns: 2
    ## ── Column specification ──────────────────────────────────────────────────────────────────
    ## Delimiter: ";"
    ## chr (2): Reference, Description
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
dsm_mdd_symptoms_reference %>%
left_join(dsm_mdd_symptoms_labels, by=c('Reference'='ref')) %>%
select(Reference, Abbreviation=abbv, Label=h, Description)
```

    ## # A tibble: 15 × 4
    ##    Reference Abbreviation Label       Description                                         
    ##    <chr>     <chr>        <chr>       <chr>                                               
    ##  1 MDD1      Dep          Mood        Depressed mood most of the day, nearly every day    
    ##  2 MDD2      Anh          Interest    Markedly diminished interest or pleasure in all, or…
    ##  3 MDD3      App          Weight⇅     Significant change in weight or appetite            
    ##  4 MDD3a     AppDec       Weight⇊     Significant weight loss or decrease in appetite     
    ##  5 MDD3b     AppInc       Weight⇈     Significant weight gain or increase in appetite     
    ##  6 MDD4      Sle          Sleep⇅      Sleeping too much or not sleeping enough            
    ##  7 MDD4a     SleDec       Sleep⇊      Insomnia nearly every day                           
    ##  8 MDD4b     SleInc       Sleep⇈      Hypersomnia nearly every day                        
    ##  9 MDD5      Moto         Motor⇅      Changes in speed/amount of moving or speaking       
    ## 10 MDD5a     MotoInc      Motor⇈      Psychomotor agitation nearly every day              
    ## 11 MDD5b     MotoDec      Motor⇊      Psychomotor slowing nearly every day                
    ## 12 MDD6      Fatig        Fatigue     Fatigue or loss of energy nearly every day          
    ## 13 MDD7      Guilt        Guilt       Feelings of worthlessness or excessive or inappropr…
    ## 14 MDD8      Conc         Concentrate Diminished ability to think or concentrate, or inde…
    ## 15 MDD9      Sui          Suicidality Recurrent thoughts of death or suicide or a suicide…

# GenomicSEM covariance structure

``` r
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 19 Columns: 9
    ## ── Column specification ──────────────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): cohorts, symptom, sumstats, filename, trait_name
    ## dbl (4): Nca, Nco, samp_prev, pop_prev
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Rename samples: AGDS/PGC is the **Clin**ical meta-analysis sample
(`Clin`) and ALSPAC/UKB is the **Comm**unity meta-analysis sample
(`Comm`); and rename symptoms numbers (`MDD1`, `MDD2`) to abbreviations
(`Dep`, `Anh`). There are also extra measures of `MDD1` and `MDD2` from
**UKB** Baseline data.

``` r
cohorts_sample_symptoms <-
sumstats_prevs %>%
left_join(dsm_mdd_symptoms_labels, by=c('symptom'='ref')) %>%
select(cohorts, symptom, trait_name, abbv) %>%
mutate(Sample=case_when(cohorts %in% 'AGDS_PGC' ~ 'Clin',
                        cohorts %in% 'ALSPAC_UKB' ~ 'Comm',
                        cohorts %in% 'UKBt' ~ 'Ukb',
                        TRUE ~ NA_character_)) %>%
mutate(sample_symptom=paste0(Sample, abbv))

sample_symptoms <- cohorts_sample_symptoms$sample_symptom
names(sample_symptoms) <- cohorts_sample_symptoms$trait_name

# rename traits in covstruct
dimnames(symptoms_covstruct$S)[[2]] <-
as.vector(sample_symptoms[dimnames(symptoms_covstruct$S)[[2]]])
```

# Structural models

Symptoms with positive variances

``` r
symptoms_S_var <- diag(symptoms_covstruct$S)
names(symptoms_S_var) <- dimnames(symptoms_covstruct$S)[[2]]

symptoms_S_var[which(symptoms_S_var > 0)]
```

    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec     ClinSui 
    ## 0.128871019 0.051543416 0.013663704 0.013924961 0.039896861 0.002375007 0.103892914 
    ##     CommDep     CommAnh  CommAppDec  CommAppInc  CommSleDec  CommSleInc   CommFatig 
    ## 0.082198484 0.087343514 0.035741853 0.079200630 0.043111949 0.047553540 0.057539570 
    ##   CommGuilt    CommConc     CommSui      UkbDep      UkbAnh 
    ## 0.060012455 0.056887041 0.033166683 0.042463388 0.047031487

## Common factor

Common factor across symptoms from both cohorts, as a general MDD factor

``` r
commonfactor.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

MDD ~~ 1*MDD
"

commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.331 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
commonfactor.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC      CFI      SRMR
    ## df 1336.652 152 1.909226e-188 1412.652 0.950653 0.1751983

``` r
commonfactor.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs == 'MDD') %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op         rhs STD_Genotype STD_Genotype_SE
    ## 1  MDD =~  ClinAppDec        0.031           0.075
    ## 2  MDD =~  ClinAppInc       -0.056           0.104
    ## 3  MDD =~  ClinSleDec       -0.095           0.127
    ## 4  MDD =~  ClinSleInc       -0.135           0.143
    ## 5  MDD =~ ClinMotoInc       -0.178           0.107
    ## 6  MDD =~ ClinMotoDec       -0.169           0.165
    ## 7  MDD =~     ClinSui       -0.483           0.094
    ## 8  MDD =~     CommDep       -0.903           0.047
    ## 9  MDD =~     CommAnh       -0.981           0.039
    ## 10 MDD =~  CommAppDec       -0.161           0.085
    ## 11 MDD =~  CommAppInc       -0.321           0.076
    ## 12 MDD =~  CommSleDec       -0.530           0.092
    ## 13 MDD =~  CommSleInc       -0.463           0.102
    ## 14 MDD =~   CommFatig       -0.575           0.091
    ## 15 MDD =~   CommGuilt       -0.584           0.080
    ## 16 MDD =~    CommConc       -0.682           0.085
    ## 17 MDD =~     CommSui       -0.550           0.089
    ## 18 MDD =~      UkbDep       -0.861           0.063
    ## 19 MDD =~      UkbAnh       -0.889           0.056
    ## 20 MDD ~~         MDD        1.000              NA

## Directional symptoms

Common factor with residual correlations among paired directional
symptoms.

``` r
common_dir.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

ClinAppDec ~~ ClinAppInc
ClinSleDec ~~ ClinSleInc
ClinMotoInc ~~ ClinMotoDec
CommAppDec ~~ CommAppInc
CommSleDec ~~ CommSleInc

MDD ~~ 1*MDD
"

common_dir.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=common_dir.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.388 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = common_dir.model):
    ## A difference greater than .025 was observed pre- and post-smoothing in the genetic
    ## covariance matrix. This reflects a large difference and results should be interpreted with
    ## caution!! This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = common_dir.model):
    ## A difference greater than .025 was observed pre- and post-smoothing for Z-statistics
    ## in the genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered traits,
    ## and you might consider removing those traits from the model. If you are going to run a
    ## multivariate GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
common_dir.fit$modelfit
```

    ##     chisq  df       p_chisq    AIC       CFI      SRMR
    ## df 1222.6 147 7.177833e-169 1308.6 0.9551956 0.1664624

``` r
common_dir.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs == 'MDD' | (lhs != 'MDD' & lhs != rhs)) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##            lhs op         rhs STD_Genotype STD_Genotype_SE
    ## 1          MDD =~  ClinAppDec        0.029           0.075
    ## 2          MDD =~  ClinAppInc       -0.056           0.104
    ## 3          MDD =~  ClinSleDec       -0.095           0.126
    ## 4          MDD =~  ClinSleInc       -0.134           0.143
    ## 5          MDD =~ ClinMotoInc       -0.180           0.107
    ## 6          MDD =~ ClinMotoDec       -0.173           0.165
    ## 7          MDD =~     ClinSui       -0.483           0.094
    ## 8          MDD =~     CommDep       -0.902           0.047
    ## 9          MDD =~     CommAnh       -0.980           0.038
    ## 10         MDD =~  CommAppDec       -0.171           0.085
    ## 11         MDD =~  CommAppInc       -0.325           0.076
    ## 12         MDD =~  CommSleDec       -0.542           0.092
    ## 13         MDD =~  CommSleInc       -0.477           0.102
    ## 14         MDD =~   CommFatig       -0.575           0.091
    ## 15         MDD =~   CommGuilt       -0.583           0.080
    ## 16         MDD =~    CommConc       -0.682           0.085
    ## 17         MDD =~     CommSui       -0.550           0.089
    ## 18         MDD =~      UkbDep       -0.860           0.063
    ## 19         MDD =~      UkbAnh       -0.888           0.056
    ## 20  ClinAppDec ~~  ClinAppInc       -0.423           0.185
    ## 21  ClinSleDec ~~  ClinSleInc        0.383           0.459
    ## 22 ClinMotoInc ~~ ClinMotoDec       -0.340           0.424
    ## 23  CommAppDec ~~  CommAppInc       -0.191           0.134
    ## 24  CommSleDec ~~  CommSleInc       -0.275           0.185
    ## 25         MDD ~~         MDD        1.000              NA

## Ascertainment-specific factors

### Clinical-Community

Symptoms were assessed in Clinical and Community samples. Make factors
representing sample type.

``` r
clin_comm.model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

CLIN ~~ 1*CLIN
COMM ~~ 1*COMM

"

clin_comm.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_comm.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.531 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_comm.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1327.448 151 3.816421e-187 1405.448 0.9509947 0.1749241

``` r
clin_comm.fit$results[c(1,2,3,6,7,9)] %>%
     filter(lhs %in% c('CLIN', 'COMM'), rhs %in% c('CLIN', 'COMM'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE p_value
    ## 1 CLIN ~~ COMM         0.72            0.37   0.052

### Gating symptoms

Then consider the community cardinal symptoms (Depression and
Anhedonia), which are gating items to the rest of the community-sample
symptoms.

``` r
gate.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

MDD ~~ 1*MDD
GATE ~~ 1*GATE
MDD ~~ 0*GATE
"

gate.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=gate.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.384 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = gate.model):
    ## A difference greater than .025 was observed pre- and post-smoothing in the
    ## genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered
    ## traits, and you might consider removing those traits from the model. If you are
    ## going to run a multivariate GWAS we strongly recommend setting the smooth_check
    ## argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = gate.model):
    ## A difference greater than .025 was observed pre- and post-smoothing for
    ## Z-statistics in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

``` r
gate.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI     SRMR
    ## df 1150.775 148 9.954511e-155 1234.775 0.9582291 0.162555

``` r
gate.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('GATE')) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op     rhs STD_Genotype STD_Genotype_SE
    ## 1 GATE =~ CommDep         0.65           0.090
    ## 2 GATE =~ CommAnh         0.54           0.089
    ## 3 GATE =~  UkbDep         0.86           0.112
    ## 4 GATE =~  UkbAnh         0.59           0.110
    ## 5 GATE ~~    GATE         1.00              NA

### Gate-Community-Clinical

Comparison of ascertainment and measurement. Gating symptoms from
community sample (distinguish controls from subthreshold), other
symptoms from community sample symptoms (subthreshold from cases), and
clinical cohort symptoms (distinguish cases from each other)

``` r
measure.model <- "
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh
COMM =~ NA*CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui 

CLIN ~~ 1*CLIN
GATE ~~ 1*GATE
COMM ~~ 1*COMM
cl_co < 1
CLIN ~~ cl_co*COMM
"

measure.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=measure.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   9.436 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## measure.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## measure.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
measure.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1480.445 149 7.669439e-218 1562.445 0.9445383 0.1595514

``` r
measure.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('CLIN', 'COMM', 'GATE'), rhs %in% c('CLIN', 'COMM', 'GATE'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 CLIN ~~ COMM         1.00           0.348
    ## 2 CLIN ~~ GATE         0.44           0.177
    ## 3 GATE ~~ COMM         0.79           0.073

### MDD Subtypes

Group symptoms from the clinical cohorts (that had positive genetic
variance) together with the same symptoms from the community cohorts (to
create a dimension of MDD subtypes), then a separate factor for all
other community cohort symptoms (that separates cases from controls)

``` r
subtype.model <- "
SUBTYPE =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommSui
MDD =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommFatig + CommGuilt + CommConc 

SUBTYPE ~~ 1*SUBTYPE
MDD ~~ 1*MDD
"

subtype.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=subtype.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.739 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = subtype.model):
    ## A difference greater than .025 was observed pre- and post-smoothing in the genetic
    ## covariance matrix. This reflects a large difference and results should be interpreted with
    ## caution!! This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = subtype.model):
    ## A difference greater than .025 was observed pre- and post-smoothing for Z-statistics
    ## in the genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered traits,
    ## and you might consider removing those traits from the model. If you are going to run a
    ## multivariate GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
subtype.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1324.772 151 1.252082e-186 1402.772 0.9511062 0.1739938

``` r
subtype.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('SUBTYPE', 'MDD'), rhs %in% c('SUBTYPE', 'MDD'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##       lhs op rhs STD_Genotype STD_Genotype_SE
    ## 1 SUBTYPE ~~ MDD        -0.86           0.099

## Two-factor models

[Elhai Psychiat
Res 2012](https://www.sciencedirect.com/science/article/pii/S0165178112002685)
compared 3 two-factor models

### Psychological-Somatic (Elhai Model 2a)

[Krause Rehab
Psychol 2008](https://psycnet.apa.org/record/2008-17022-011), [Krause
Arch Psys Med
Rehab 2010](https://www.sciencedirect.com/science/article/pii/S0003999310002443):

> the 2-factor solution with 3 somatic items (sleep disturbance, poor
> energy, appetite change) was a better solution than either a
> unidimensional model or 2-factor model that included psychomotor
> slowing as a fourth somatic item

``` r
psych_soma.model <- "
PSYCH =~ NA*ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
SOMA =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

PSYCH ~~ 1*PSYCH
SOMA ~~ 1*SOMA
GATE ~~ 1*GATE
GATE ~~ 0*PSYCH + 0*SOMA
"

psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    1.65 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
psych_soma.fit$modelfit
```

    ##     chisq  df       p_chisq    AIC       CFI      SRMR
    ## df 1136.5 147 1.809352e-152 1222.5 0.9587821 0.1625164

``` r
psych_soma.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'SOMA'), rhs %in% c('PSYCH', 'SOMA'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ SOMA        -0.98            0.13

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
psych_veg.model <- "
PSYCH =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
VEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

PSYCH ~~ 1*PSYCH
VEG ~~ 1*VEG
GATE ~~ 1*GATE
GATE ~~ 0*PSYCH + 0*VEG
"
psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.807 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
psych_veg.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC      CFI      SRMR
    ## df 1081.984 147 3.556402e-142 1167.984 0.961053 0.1614363

``` r
psych_veg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'VEG'), rhs %in% c('PSYCH', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ VEG         -0.9            0.11

### Affective-Neurovegetative (Elhai Model 2c)

``` r
affect_neuroveg.model <- "
AFFECT =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommSui + UkbDep + UkbAnh
NEUROVEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommConc

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

AFFECT ~~ 1*AFFECT
NEUROVEG ~~ 1*NEUROVEG
GATE ~~ 1*GATE
GATE ~~ 0*AFFECT + 0*NEUROVEG
"
affect_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=affect_neuroveg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.607 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
affect_neuroveg.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 1148.942 147 7.90408e-155 1234.942 0.9582638 0.1608248

``` r
affect_neuroveg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('AFFECT', 'NEUROVEG'), rhs %in% c('AFFECT', 'NEUROVEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##      lhs op      rhs STD_Genotype STD_Genotype_SE
    ## 1 AFFECT ~~ NEUROVEG        -0.89             0.1

## Melancholic and atypical

Account for directionality of symptoms using melancholic and atypical
classifications

``` r
mel_aty.model <- "
MEL =~ NA*CommAnh + UkbAnh + ClinSleDec + CommSleDec + ClinMotoInc + ClinAppDec + CommAppDec + CommGuilt
ATY =~ NA*ClinAppInc + CommAppInc + ClinSleInc + CommSleInc + ClinMotoDec + CommFatig
MDD =~ NA*CommDep + UkbDep + ClinSui + CommConc + CommSui 
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

MEL ~~ 1*MEL
ATY ~~ 1*ATY
MDD ~~ 1*MDD
GATE ~~ 1*GATE
GATE ~~ 0*MDD + 0*MEL + 0*ATY
"

mel_aty.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=mel_aty.model, imp_cov=TRUE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.898 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = mel_aty.model, :
    ## A difference greater than .025 was observed pre- and post-smoothing in the genetic
    ## covariance matrix. This reflects a large difference and results should be interpreted with
    ## caution!! This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = mel_aty.model, :
    ## A difference greater than .025 was observed pre- and post-smoothing for Z-statistics
    ## in the genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered traits,
    ## and you might consider removing those traits from the model. If you are going to run a
    ## multivariate GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
mel_aty.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1044.031 145 6.488099e-136 1134.031 0.9625506 0.1588098

``` r
mel_aty.fit$results[c(1,2,3,6,7)] %>%
     filter(rhs %in% c('MEL', 'ATY', 'MDD'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##   lhs op rhs STD_Genotype STD_Genotype_SE
    ## 1 MEL ~~ ATY         0.84           0.125
    ## 2 MEL ~~ MDD         0.99           0.039
    ## 3 ATY ~~ MDD         0.69           0.129

## Three factor models

### Cognitive-Mood-Neuroveg (Kendler Neale) model

``` r
cog_mood_neuroveg.model <- "
COG =~ NA*ClinSui + CommGuilt + CommConc + CommSui
MOOD =~ NA*CommDep + CommAnh + CommGuilt + UkbDep + UkbAnh
VEG =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

COG ~~ 1*COG
MOOD ~~ 1*MOOD
VEG ~~ 1*VEG
GATE ~~ 1*GATE
GATE ~~ 0*COG + 0*MOOD + 0*VEG
"
cog_mood_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=cog_mood_neuroveg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.866 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_neuroveg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
cog_mood_neuroveg.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1088.314 144 1.104125e-144 1180.314 0.9606643 0.1615335

``` r
cog_mood_neuroveg.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('COG', 'MOOD', 'VEG'), rhs %in% c('COG', 'MOOD', 'VEG'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1  COG ~~ MOOD         0.86            0.12
    ## 2  COG ~~  VEG        -0.87            0.12
    ## 3 MOOD ~~  VEG        -0.79            0.13

### Psych/Somatic + Atypical

Expand on the Psych/Somatic two-factor models by separating out
directional symptoms.

#### Psychological-Somatic-Atypical

``` r
psych_soma_aty.model <- "
PSYCH =~ NA*ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
SOMA =~ NA*ClinAppDec + ClinSleDec + CommAppDec + CommSleDec
ATY =~ NA*ClinAppInc + ClinSleInc + CommAppInc + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

PSYCH ~~ 1*PSYCH
SOMA ~~ 1*SOMA
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*PSYCH + 0*SOMA + 0*ATY
"

psych_soma_aty.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_soma_aty.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.836 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma_aty.model): A difference greater than .025 was observed pre- and post-smoothing
    ## in the genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered traits,
    ## and you might consider removing those traits from the model. If you are going to run a
    ## multivariate GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma_aty.model): A difference greater than .025 was observed pre- and post-smoothing
    ## for Z-statistics in the genetic covariance matrix. This reflects a large difference and
    ## results should be interpreted with caution!! This can often result from including low
    ## powered traits, and you might consider removing those traits from the model. If you are
    ## going to run a multivariate GWAS we strongly recommend setting the smooth_check argument
    ## to true to check smoothing for each SNP.

``` r
psych_soma_aty.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1235.313 145 3.086172e-172 1325.313 0.9545827 0.1548222

``` r
psych_soma_aty.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'SOMA', 'ATY'), rhs %in% c('PSYCH', 'SOMA', 'ATY'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ SOMA         0.46            0.16
    ## 2 PSYCH ~~  ATY         0.82            0.13
    ## 3  SOMA ~~  ATY         0.20            0.15

#### Psychological-Neurovegetative-Atypical

``` r
psych_veg_aty.model <- "
PSYCH =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
VEG =~ NA*ClinAppDec + ClinSleDec + ClinMotoInc +  + CommAppDec + CommSleDec
ATY =~ NA*ClinAppInc + ClinSleInc + ClinMotoDec + CommAppInc + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

PSYCH ~~ 1*PSYCH
VEG ~~ 1*VEG
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*PSYCH + 0*VEG + ATY
"
psych_veg_aty.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_veg_aty.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.692 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg_aty.model): A difference greater than .025 was observed pre- and post-smoothing
    ## in the genetic covariance matrix. This reflects a large difference and results should
    ## be interpreted with caution!! This can often result from including low powered traits,
    ## and you might consider removing those traits from the model. If you are going to run a
    ## multivariate GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_veg_aty.model): A difference greater than .025 was observed pre- and post-smoothing
    ## for Z-statistics in the genetic covariance matrix. This reflects a large difference and
    ## results should be interpreted with caution!! This can often result from including low
    ## powered traits, and you might consider removing those traits from the model. If you are
    ## going to run a multivariate GWAS we strongly recommend setting the smooth_check argument
    ## to true to check smoothing for each SNP.

``` r
psych_veg_aty.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1192.665 144 1.588989e-164 1284.665 0.9563176 0.1511433

``` r
psych_veg_aty.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('PSYCH', 'VEG', 'ATY'), rhs %in% c('PSYCH', 'VEG', 'ATY'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##     lhs op rhs STD_Genotype STD_Genotype_SE
    ## 1 PSYCH ~~ VEG         0.47            0.13
    ## 2 PSYCH ~~ ATY         0.84            0.13
    ## 3   VEG ~~ ATY         0.26            0.16

#### Affective-Neurovegetative-Atypical

``` r
affect_neuroveg_aty.model <- "
AFFECT =~ NA*ClinSui + CommDep + CommAnh + CommGuilt + CommSui + UkbDep + UkbAnh
NEUROVEG =~ NA*ClinAppDec + ClinSleDec + ClinMotoInc + CommAppDec + CommSleDec + CommConc
ATY =~ NA*ClinAppInc + ClinSleInc + ClinMotoDec + CommAppInc + CommSleInc + CommFatig

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

AFFECT ~~ 1*AFFECT
NEUROVEG ~~ 1*NEUROVEG
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*AFFECT + 0*NEUROVEG + 0*ATY
"
affect_neuroveg_aty.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=affect_neuroveg_aty.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.905 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference and
    ## results should be interpreted with caution!! This can often result from including low
    ## powered traits, and you might consider removing those traits from the model. If you are
    ## going to run a multivariate GWAS we strongly recommend setting the smooth_check argument
    ## to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## affect_neuroveg_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from the model.
    ## If you are going to run a multivariate GWAS we strongly recommend setting the smooth_check
    ## argument to true to check smoothing for each SNP.

``` r
affect_neuroveg_aty.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1179.509 145 1.494031e-161 1269.509 0.9569072 0.1570791

``` r
affect_neuroveg_aty.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('AFFECT', 'NEUROVEG', 'ATY'), rhs %in% c('AFFECT', 'NEUROVEG', 'ATY'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##        lhs op      rhs STD_Genotype STD_Genotype_SE
    ## 1   AFFECT ~~ NEUROVEG         0.82            0.13
    ## 2   AFFECT ~~      ATY         0.76            0.13
    ## 3 NEUROVEG ~~      ATY         0.66            0.16

## Four factor models

### Cognitive-Mood-Neuroveg (Kendler Neale) model, plus Melancholic/Atypical split

``` r
cog_mood_mel_aty.model <- "
COG =~ NA*ClinSui + CommGuilt + CommConc + CommSui
MOOD =~ NA*CommDep + CommAnh + CommGuilt + UkbDep + UkbAnh
ATY =~ NA*ClinMotoDec + ClinAppInc +  + ClinSleInc + CommAppInc + CommSleInc + CommFatig
MEL =~ NA*ClinMotoInc + ClinAppDec + CommAppDec + ClinSleDec + CommSleDec

GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

COG ~~ 1*COG
MOOD ~~ 1*MOOD
MEL ~~ 1*MEL
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*COG + 0*MOOD + 0*MEL + 0*ATY
"
cog_mood_mel_aty.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=cog_mood_mel_aty.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.812 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_mel_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from the
    ## model. If you are going to run a multivariate GWAS we strongly recommend setting
    ## the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## cog_mood_mel_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
cog_mood_mel_aty.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1271.655 141 3.898009e-181 1369.655 0.9529022 0.1518656

``` r
cog_mood_mel_aty.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('COG', 'MOOD', 'MEL', 'ATY'), rhs %in% c('COG', 'MOOD', 'MEL', 'ATY'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1  COG ~~ MOOD         0.86            0.12
    ## 2  COG ~~  ATY         0.74            0.14
    ## 3  COG ~~  MEL         0.47            0.14
    ## 4 MOOD ~~  ATY         0.69            0.13
    ## 5 MOOD ~~  MEL         0.43            0.13
    ## 6  ATY ~~  MEL         0.26            0.16

### Model comparisons

``` r
model_list=
list("1a"=list(name="Common", model=commonfactor.fit),
     "1b"=list(name="Directional", model=common_dir.fit),
     "1c"=list(name="Ascertainment", model=clin_comm.fit),
     "1d"=list(name="Gating", model=gate.fit),
     "1e"=list(name="Measurement", model=measure.fit),
     "1f"=list(name="Subtype", model=subtype.fit),
     "2a"=list(name="Psych-Somatic", model=psych_soma.fit),
     "2b"=list(name="Psych-Neuroveg", model=psych_veg.fit),
     "2c"=list(name="Affect-Neuroveg", model=affect_neuroveg.fit),
     "2d"=list(name="Melancholic/Atypical", model=mel_aty.fit),
     "3"=list(name="Cog-Mood-Neuroveg", model=cog_mood_neuroveg.fit),
     "3a"=list(name="Psych-Somatic-Atypical", model=psych_soma_aty.fit),
     "3b"=list(name="Psych-Neuroveg-Atypical", model=psych_veg_aty.fit),
     "3c"=list(name="Affect-Neuroveg-Atypical", model=affect_neuroveg_aty.fit),
     "4"=list(names="Cog-Mood-Mel-Aty", model=cog_mood_mel_aty.fit)
     )

model_fits <- 
data.frame(Model=names(model_list),
           Name=sapply(model_list, function(m) m$name)) %>%
bind_cols(
bind_rows(
lapply(model_list, function(m) m$model$modelfit)
))
rownames(model_fits) <- NULL

knitr::kable(
model_fits %>%
mutate(dAIC=AIC-min(AIC)) %>%
mutate_if(is.numeric, ~signif(., 3))
)
```

| Model | Name                     | chisq |  df | p\_chisq |  AIC |   CFI |  SRMR |  dAIC |
| :---- | :----------------------- | ----: | --: | -------: | ---: | ----: | ----: | ----: |
| 1a    | Common                   |  1340 | 152 |        0 | 1410 | 0.951 | 0.175 | 279.0 |
| 1b    | Directional              |  1220 | 147 |        0 | 1310 | 0.955 | 0.166 | 175.0 |
| 1c    | Ascertainment            |  1330 | 151 |        0 | 1410 | 0.951 | 0.175 | 271.0 |
| 1d    | Gating                   |  1150 | 148 |        0 | 1230 | 0.958 | 0.163 | 101.0 |
| 1e    | Measurement              |  1480 | 149 |        0 | 1560 | 0.945 | 0.160 | 428.0 |
| 1f    | Subtype                  |  1320 | 151 |        0 | 1400 | 0.951 | 0.174 | 269.0 |
| 2a    | Psych-Somatic            |  1140 | 147 |        0 | 1220 | 0.959 | 0.163 |  88.5 |
| 2b    | Psych-Neuroveg           |  1080 | 147 |        0 | 1170 | 0.961 | 0.161 |  34.0 |
| 2c    | Affect-Neuroveg          |  1150 | 147 |        0 | 1230 | 0.958 | 0.161 | 101.0 |
| 2d    | Melancholic/Atypical     |  1040 | 145 |        0 | 1130 | 0.963 | 0.159 |   0.0 |
| 3     | Cog-Mood-Neuroveg        |  1090 | 144 |        0 | 1180 | 0.961 | 0.162 |  46.3 |
| 3a    | Psych-Somatic-Atypical   |  1240 | 145 |        0 | 1330 | 0.955 | 0.155 | 191.0 |
| 3b    | Psych-Neuroveg-Atypical  |  1190 | 144 |        0 | 1280 | 0.956 | 0.151 | 151.0 |
| 3c    | Affect-Neuroveg-Atypical |  1180 | 145 |        0 | 1270 | 0.957 | 0.157 | 135.0 |
| 4     | Cog-Mood-Mel-Aty         |  1270 | 141 |        0 | 1370 | 0.953 | 0.152 | 236.0 |

Check that all symptoms have been included in each model

``` r
sapply(model_list, function(m) all(names(symptoms_S_var) %in% m$model$results$rhs))
```

    ##   1a   1b   1c   1d   1e   1f   2a   2b   2c   2d    3   3a   3b   3c    4 
    ## TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# Exploratory factor analysis

Get the genetic covariance matrix for symptoms with a positive
heritability

``` r
symptoms_cov <- symptoms_covstruct$S
k <- nrow(symptoms_cov)
symptoms_se <- matrix(0, k, k)
symptoms_se[lower.tri(symptoms_se, diag=TRUE)] <- sqrt(diag(symptoms_covstruct$V))

symptoms_se[upper.tri(symptoms_se)] <- t(symptoms_se)[upper.tri(symptoms_se)]

symptoms_cov_keep <- which(diag(symptoms_cov > 0))

symptoms_cov_pos <- symptoms_cov[symptoms_cov_keep,symptoms_cov_keep]
```

Smooth the genetic covariance matrix so that it is positive definite

``` r
# smooth the covariance matrix
symptoms_cov_pd <- as.matrix(Matrix::nearPD(symptoms_covstruct$S, corr=FALSE)$mat)

corrplot(cov2cor(symptoms_cov_pd))
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_pd-1.png)<!-- -->

Check eigen values of the correlation matrix

``` r
symptoms_eigen <- eigen(cov2cor(symptoms_cov_pd)) 

plot(symptoms_eigen$values, ylab='Eigenvalue')
lines(symptoms_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_eigen-1.png)<!-- -->

``` r
symptoms_efa <- factanal(covmat=symptoms_cov_pd, factors=3, rotation='varimax')
print(symptoms_efa, cut=0.4)
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec     ClinSui 
    ##       0.082       0.069       0.005       0.831       0.895       0.968       0.798 
    ##     CommDep     CommAnh  CommAppDec  CommAppInc  CommSleDec  CommSleInc   CommFatig 
    ##       0.112       0.050       0.869       0.290       0.440       0.757       0.646 
    ##   CommGuilt    CommConc     CommSui      UkbDep      UkbAnh 
    ##       0.585       0.583       0.712       0.122       0.168 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec           0.567  -0.753 
    ## ClinAppInc           0.408   0.874 
    ## ClinSleDec           0.990         
    ## ClinSleInc           0.402         
    ## ClinMotoInc                        
    ## ClinMotoDec                        
    ## ClinSui      0.439                 
    ## CommDep      0.910                 
    ## CommAnh      0.946                 
    ## CommAppDec                         
    ## CommAppInc                   0.738 
    ## CommSleDec           0.641         
    ## CommSleInc   0.475                 
    ## CommFatig    0.572                 
    ## CommGuilt    0.577                 
    ## CommConc     0.632                 
    ## CommSui      0.535                 
    ## UkbDep       0.906                 
    ## UkbAnh       0.886                 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      5.423   2.371   2.227
    ## Proportion Var   0.285   0.125   0.117
    ## Cumulative Var   0.285   0.410   0.527
    ## 
    ## The degrees of freedom for the model is 117 and the fit was 76.6301
