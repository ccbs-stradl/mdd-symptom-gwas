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
    ## ── Column specification ───────────────────────────────────────────────────────────
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
    ## ── Column specification ───────────────────────────────────────────────────────────
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
    ##  1 MDD1      Dep          Mood        Depressed mood most of the day, nearly every…
    ##  2 MDD2      Anh          Interest    Markedly diminished interest or pleasure in …
    ##  3 MDD3      App          Weight⇅     Significant change in weight or appetite     
    ##  4 MDD3a     AppDec       Weight⇊     Significant weight loss or decrease in appet…
    ##  5 MDD3b     AppInc       Weight⇈     Significant weight gain or increase in appet…
    ##  6 MDD4      Sle          Sleep⇅      Sleeping too much or not sleeping enough     
    ##  7 MDD4a     SleDec       Sleep⇊      Insomnia nearly every day                    
    ##  8 MDD4b     SleInc       Sleep⇈      Hypersomnia nearly every day                 
    ##  9 MDD5      Moto         Motor⇅      Changes in speed/amount of moving or speaking
    ## 10 MDD5a     MotoInc      Motor⇈      Psychomotor agitation nearly every day       
    ## 11 MDD5b     MotoDec      Motor⇊      Psychomotor slowing nearly every day         
    ## 12 MDD6      Fatig        Fatigue     Fatigue or loss of energy nearly every day   
    ## 13 MDD7      Guilt        Guilt       Feelings of worthlessness or excessive or in…
    ## 14 MDD8      Conc         Concentrate Diminished ability to think or concentrate, …
    ## 15 MDD9      Sui          Suicidality Recurrent thoughts of death or suicide or a …

# GenomicSEM covariance structure

``` r
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 19 Columns: 9
    ## ── Column specification ───────────────────────────────────────────────────────────
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

    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec 
    ## 0.128871019 0.051543416 0.013663704 0.013924961 0.039896861 0.002375007 
    ##     ClinSui     CommDep     CommAnh  CommAppDec  CommAppInc  CommSleDec 
    ## 0.103892914 0.082198484 0.087343514 0.035741853 0.079200630 0.043111949 
    ##  CommSleInc   CommFatig   CommGuilt    CommConc     CommSui      UkbDep 
    ## 0.047553540 0.057539570 0.060012455 0.056887041 0.033166683 0.042463388 
    ##      UkbAnh 
    ## 0.047031487

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

## Correlation among directional symptoms

Directional symptoms that can increase or decrease (appetite/weight
changes, sleep, psychomotor) will likely be negatively correlated, so
add an indepdendent factor to capture this.

First formulation is an independent factor for the directional symptom.

``` r
commonfactor_app.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc

MDD ~~ 1*MDD
APP ~~ 1*APP
MDD ~~ 0*APP
"

commonfactor_app.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.486 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_app.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_app.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
commonfactor_app.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1191.863 148 1.534121e-162 1275.863 0.9565176 0.1578159

``` r
commonfactor_app.fit$results[c(1,2,3,6,7,9)] %>%
     filter(lhs %in% c('MDD', 'APP')) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op         rhs STD_Genotype STD_Genotype_SE  p_value
    ## 1  MDD =~  ClinAppDec       0.0079           0.075  9.2e-01
    ## 2  MDD =~  ClinAppInc      -0.0279           0.104  7.9e-01
    ## 3  MDD =~  ClinSleDec      -0.0957           0.127  4.5e-01
    ## 4  MDD =~  ClinSleInc      -0.1347           0.143  3.5e-01
    ## 5  MDD =~ ClinMotoInc      -0.1785           0.107  9.5e-02
    ## 6  MDD =~ ClinMotoDec      -0.1691           0.165  3.1e-01
    ## 7  MDD =~     ClinSui      -0.4833           0.094  3.1e-07
    ## 8  MDD =~     CommDep      -0.9044           0.047  4.6e-82
    ## 9  MDD =~     CommAnh      -0.9805           0.039 2.0e-142
    ## 10 MDD =~  CommAppDec      -0.1773           0.085  3.8e-02
    ## 11 MDD =~  CommAppInc      -0.3211           0.076  2.7e-05
    ## 12 MDD =~  CommSleDec      -0.5303           0.092  7.2e-09
    ## 13 MDD =~  CommSleInc      -0.4626           0.102  5.8e-06
    ## 14 MDD =~   CommFatig      -0.5744           0.091  3.0e-10
    ## 15 MDD =~   CommGuilt      -0.5820           0.080  2.9e-13
    ## 16 MDD =~    CommConc      -0.6817           0.085  8.3e-16
    ## 17 MDD =~     CommSui      -0.5502           0.090  7.9e-10
    ## 18 MDD =~      UkbDep      -0.8622           0.063  1.3e-42
    ## 19 MDD =~      UkbAnh      -0.8890           0.056  4.4e-57
    ## 20 APP =~  ClinAppDec       0.6207           0.126  8.9e-07
    ## 21 APP =~  ClinAppInc      -0.8232           0.174  2.3e-06
    ## 22 APP =~  CommAppDec       0.3437           0.128  7.2e-03
    ## 23 APP =~  CommAppInc      -0.7992           0.180  9.4e-06
    ## 24 MDD ~~         MDD       1.0000              NA       NA
    ## 25 APP ~~         APP       1.0000              NA       NA
    ## 26 MDD ~~         APP       0.0000              NA       NA

Second formulation is where the appetite symptom is a subfactor that
then loads on the general MDD factor.

``` r
commonfactor_app_sub.model <- "
MDD =~ NA*ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh + APP
APP =~ NA*ClinAppDec + ClinAppInc + CommAppDec + CommAppInc

MDD ~~ 1*MDD
APP ~~ 1*APP
"

commonfactor_app_sub.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor_app_sub.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.828 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_app_sub.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_app_sub.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
commonfactor_app_sub.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1131.777 151 8.338548e-150 1209.777 0.9591455 0.1660853

``` r
commonfactor_app_sub.fit$results[c(1,2,3,6,7,9)] %>%
     filter(lhs %in% c('MDD', 'APP')) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op         rhs STD_Genotype STD_Genotype_SE  p_value
    ## 1  MDD =~  ClinSleDec        0.094           0.127  4.6e-01
    ## 2  MDD =~  ClinSleInc        0.137           0.143  3.3e-01
    ## 3  MDD =~ ClinMotoInc        0.174           0.107  1.0e-01
    ## 4  MDD =~ ClinMotoDec        0.166           0.164  3.1e-01
    ## 5  MDD =~     ClinSui        0.480           0.094  3.5e-07
    ## 6  MDD =~     CommDep        0.904           0.047  6.2e-83
    ## 7  MDD =~     CommAnh        0.982           0.038 7.8e-144
    ## 8  MDD =~  CommSleDec        0.525           0.093  1.4e-08
    ## 9  MDD =~  CommSleInc        0.465           0.102  5.1e-06
    ## 10 MDD =~   CommFatig        0.574           0.091  2.7e-10
    ## 11 MDD =~   CommGuilt        0.583           0.080  2.5e-13
    ## 12 MDD =~    CommConc        0.681           0.085  9.1e-16
    ## 13 MDD =~     CommSui        0.548           0.089  7.4e-10
    ## 14 MDD =~      UkbDep        0.865           0.063  4.9e-43
    ## 15 MDD =~      UkbAnh        0.890           0.055  6.4e-58
    ## 16 MDD =~         APP       -0.143           0.116  2.2e-01
    ## 17 APP =~  ClinAppDec        0.228           0.177  2.0e-01
    ## 18 APP =~  ClinAppInc       -0.332           0.261  2.0e-01
    ## 19 APP =~  CommAppDec        0.050           0.072  4.9e-01
    ## 20 APP =~  CommAppInc       -2.192           1.729  2.0e-01
    ## 21 MDD ~~         MDD        1.000              NA       NA
    ## 22 APP ~~         APP        1.000              NA       NA

``` r
commonfactor_sle.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
SLE =~ NA*ClinSleDec + ClinSleInc + CommSleDec + CommSleInc

MDD ~~ 1*MDD
SLE ~~ 1*SLE
MDD ~~ 0*SLE

c4a > 0
CommSleDec ~~ c4a*CommSleDec
"

commonfactor_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   9.163 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_sle.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_sle.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
commonfactor_sle.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 1301.647 148 1.36484e-183 1385.647 0.9519445 0.1697398

``` r
commonfactor_sle.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('MDD', 'SLE')) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op         rhs STD_Genotype STD_Genotype_SE
    ## 1  MDD =~  ClinAppDec        0.031           0.075
    ## 2  MDD =~  ClinAppInc       -0.056           0.104
    ## 3  MDD =~  ClinSleDec       -0.076           0.128
    ## 4  MDD =~  ClinSleInc       -0.134           0.143
    ## 5  MDD =~ ClinMotoInc       -0.178           0.107
    ## 6  MDD =~ ClinMotoDec       -0.169           0.165
    ## 7  MDD =~     ClinSui       -0.483           0.094
    ## 8  MDD =~     CommDep       -0.903           0.047
    ## 9  MDD =~     CommAnh       -0.981           0.038
    ## 10 MDD =~  CommAppDec       -0.161           0.085
    ## 11 MDD =~  CommAppInc       -0.321           0.076
    ## 12 MDD =~  CommSleDec       -0.535           0.092
    ## 13 MDD =~  CommSleInc       -0.474           0.102
    ## 14 MDD =~   CommFatig       -0.575           0.091
    ## 15 MDD =~   CommGuilt       -0.584           0.080
    ## 16 MDD =~    CommConc       -0.682           0.085
    ## 17 MDD =~     CommSui       -0.550           0.089
    ## 18 MDD =~      UkbDep       -0.861           0.063
    ## 19 MDD =~      UkbAnh       -0.888           0.056
    ## 20 SLE =~  ClinSleDec        0.596           0.637
    ## 21 SLE =~  ClinSleInc        0.044           0.320
    ## 22 SLE =~  CommSleDec        0.869           0.936
    ## 23 SLE =~  CommSleInc       -0.240           0.260
    ## 24 MDD ~~         MDD        1.000              NA
    ## 25 SLE ~~         SLE        1.000              NA
    ## 26 MDD ~~         SLE        0.000              NA

## Ascertainment-specific factors

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

Add the sample factors as a bifactor model to the general MDD factor

``` r
clin_comm_bif.model <- "
MDD =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh
CLIN =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui
COMM =~ NA*CommDep + CommAnh + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommFatig + CommGuilt + CommConc + CommSui + UkbDep + UkbAnh

MDD ~~ 1*MDD
CLIN ~~ 1*CLIN
COMM ~~ 1*COMM

MDD ~~ 0*CLIN + 0*COMM
"

clin_comm_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_comm_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.116 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_comm_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_comm_bif.fit$modelfit
```

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 784.4446 132 2.521846e-93 900.4446 0.9728222 0.1475522

``` r
clin_comm_bif.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('CLIN', 'COMM'), rhs %in% c('CLIN', 'COMM'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op  rhs STD_Genotype STD_Genotype_SE
    ## 1 CLIN ~~ COMM         0.64            0.33

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

Group symptoms from the clinical cohorts (that had positive genetic
variance) together with the same symptoms from the community cohorts,
then a separate factor for all other community cohort symptoms

``` r
divergence.model <- "
CASE =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinMotoDec + ClinSui + CommAppDec + CommAppInc + CommSleDec + CommSleInc + CommSui
CASECON =~ NA*CommDep + CommAnh + UkbDep + UkbAnh + CommFatig + CommGuilt + CommConc 

CASE ~~ 1*CASE
CASECON ~~ 1*CASECON
"

divergence.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=divergence.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.561 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## divergence.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## divergence.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
divergence.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1324.772 151 1.252082e-186 1402.772 0.9511062 0.1739938

``` r
divergence.fit$results[c(1,2,3,6,7)] %>%
     filter(lhs %in% c('CASE', 'CASECON'), rhs %in% c('CASE', 'CASECON'), lhs != rhs) %>%
     mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     print(digits=2)
```

    ##    lhs op     rhs STD_Genotype STD_Genotype_SE
    ## 1 CASE ~~ CASECON        -0.86           0.099

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

mel_aty.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=mel_aty.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.703 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## mel_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## mel_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
mel_aty.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1044.031 145 6.488099e-136 1134.031 0.9625506 0.1588098

``` r
mel_aty.fit$results[c(1,2,3,6,7)]# %>%
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE
    ## 42         MEL =~     CommAnh   0.80192727 0.0666090638138897
    ## 46         MEL =~      UkbAnh   0.67262425 0.0880081685363808
    ## 41         MEL =~  ClinSleDec   0.11420023  0.136430372252711
    ## 45         MEL =~  CommSleDec   0.62015753   0.10611794870195
    ## 40         MEL =~ ClinMotoInc   0.20828075  0.119288111961555
    ## 39         MEL =~  ClinAppDec  -0.06267036 0.0855513086310624
    ## 43         MEL =~  CommAppDec   0.18427682 0.0920449027383444
    ## 44         MEL =~   CommGuilt   0.68227142 0.0907922995390054
    ## 1          ATY =~  ClinAppInc   0.17800744  0.139447015737535
    ## 4          ATY =~  CommAppInc   0.47577795  0.105751413828269
    ## 3          ATY =~  ClinSleInc   0.24260177  0.188899426139744
    ## 6          ATY =~  CommSleInc   0.63192098  0.146325435444057
    ## 2          ATY =~ ClinMotoDec   0.31086862  0.223425351375137
    ## 5          ATY =~   CommFatig   0.85530948  0.137706767813619
    ## 34         MDD =~     CommDep   0.71184016 0.0716735595163633
    ## 36         MDD =~      UkbDep   0.55974425  0.097171910598007
    ## 32         MDD =~     ClinSui   0.59402370  0.109470359890207
    ## 33         MDD =~    CommConc   0.81818514  0.103916366597557
    ## 35         MDD =~     CommSui   0.64020020    0.1016418712075
    ## 28        GATE =~     CommDep   0.61738833 0.0973409455001856
    ## 27        GATE =~     CommAnh   0.53497918 0.0918420235876619
    ## 30        GATE =~      UkbDep   0.84040643  0.117144200431505
    ## 29        GATE =~      UkbAnh   0.60186776   0.11085186006899
    ## 17     CommAnh ~~     CommAnh   0.07070978 0.0504892995511117
    ## 51      UkbAnh ~~      UkbAnh   0.18532916 0.0838793375991046
    ## 14  ClinSleDec ~~  ClinSleDec   0.98695900  0.583175133715307
    ## 24  CommSleDec ~~  CommSleDec   0.61540429  0.295121904660136
    ## 13 ClinMotoInc ~~ ClinMotoInc   0.95662021  0.392550984419721
    ## 10  ClinAppDec ~~  ClinAppDec   0.99607191    0.2187761172164
    ## 18  CommAppDec ~~  CommAppDec   0.96604112  0.238859935745877
    ## 23   CommGuilt ~~   CommGuilt   0.53450617   0.16449030753765
    ## 11  ClinAppInc ~~  ClinAppInc   0.96831518  0.317327938503308
    ## 19  CommAppInc ~~  CommAppInc   0.77363512  0.153481190779067
    ## 15  ClinSleInc ~~  ClinSleInc   0.94114511  0.668312260667035
    ## 25  CommSleInc ~~  CommSleInc   0.60067727  0.253268192601536
    ## 12 ClinMotoDec ~~ ClinMotoDec   0.90335597   1.08571779423162
    ## 22   CommFatig ~~   CommFatig   0.26844441  0.314361145623672
    ## 21     CommDep ~~     CommDep   0.11211379 0.0576587345809116
    ## 52      UkbDep ~~      UkbDep  -0.01959589  0.122324962131421
    ## 16     ClinSui ~~     ClinSui   0.64713559  0.325262920738735
    ## 20    CommConc ~~    CommConc   0.33057229  0.251467046775068
    ## 26     CommSui ~~     CommSui   0.59014411  0.233070482238974
    ## 47         MEL ~~         ATY   0.83555134  0.124836933831232
    ## 49         MEL ~~         MDD   0.99028034 0.0391262518266257
    ## 9          ATY ~~         MDD   0.68620636  0.128574943855314
    ## 50         MEL ~~         MEL   1.00000000                   
    ## 7          ATY ~~         ATY   1.00000000                   
    ## 38         MDD ~~         MDD   1.00000000                   
    ## 31        GATE ~~        GATE   1.00000000                   
    ## 37         MDD ~~        GATE   0.00000000                   
    ## 48         MEL ~~        GATE   0.00000000                   
    ## 8          ATY ~~        GATE   0.00000000

``` r
     #filter(lhs == 'MDD') %>%
     #mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     #print(digits=2)
```

``` r
mel_aty_psych_soma.model <- "
AFFECT =~ NA*CommAnh + UkbAnh + CommGuilt + CommDep + UkbDep + ClinSui + CommSui
SOMA =~ NA*ClinSleDec + CommSleDec + ClinMotoInc + ClinAppDec + CommAppDec + CommConc
ATY =~ NA*ClinAppInc + CommAppInc + ClinSleInc + CommSleInc + ClinMotoDec + CommFatig
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

AFFECT ~~ 1*AFFECT
SOMA ~~ 1*SOMA
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*AFFECT + 0*SOMA + 0*ATY

u1 > 0
UkbDep ~~ u1*UkbDep
"

mel_aty_psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=mel_aty_psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  10.173 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## mel_aty_psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from the
    ## model. If you are going to run a multivariate GWAS we strongly recommend setting
    ## the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## mel_aty_psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
mel_aty_psych_soma.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC       CFI      SRMR
    ## df 1166.902 145 3.791703e-159 1256.902 0.9574324 0.1570993

``` r
mel_aty_psych_soma.fit$results[c(1,2,3,6,7)]# %>%
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE
    ## 2       AFFECT =~     CommAnh  8.305171e-01 0.0761816723933389
    ## 6       AFFECT =~      UkbAnh  6.989384e-01 0.0958558827593718
    ## 4       AFFECT =~   CommGuilt  6.864001e-01 0.0901235253473931
    ## 3       AFFECT =~     CommDep  7.025352e-01 0.0769641504035171
    ## 7       AFFECT =~      UkbDep  5.582458e-01 0.0987148808935767
    ## 1       AFFECT =~     ClinSui  5.808986e-01  0.106935388452116
    ## 5       AFFECT =~     CommSui  6.250716e-01 0.0989807843268335
    ## 44        SOMA =~  ClinSleDec  1.455061e-01  0.154903779263418
    ## 47        SOMA =~  CommSleDec  7.286289e-01   0.12809297202732
    ## 43        SOMA =~ ClinMotoInc  2.408536e-01  0.134144373928615
    ## 42        SOMA =~  ClinAppDec -5.322864e-02 0.0968209145414992
    ## 45        SOMA =~  CommAppDec  2.178311e-01  0.103646906225252
    ## 46        SOMA =~    CommConc  9.441400e-01  0.146483483157185
    ## 12         ATY =~  ClinAppInc  1.678217e-01  0.140987044678836
    ## 15         ATY =~  CommAppInc  4.718197e-01  0.106150784485376
    ## 14         ATY =~  ClinSleInc  2.448854e-01  0.190025202897218
    ## 17         ATY =~  CommSleInc  6.341677e-01  0.147252621331795
    ## 13         ATY =~ ClinMotoDec  3.096935e-01  0.223703483891982
    ## 16         ATY =~   CommFatig  8.561414e-01  0.138828225315499
    ## 38        GATE =~     CommDep  6.235158e-01  0.107279213144974
    ## 37        GATE =~     CommAnh  5.040257e-01  0.111384232822622
    ## 40        GATE =~      UkbDep  8.387049e-01  0.120002646942357
    ## 39        GATE =~      UkbAnh  5.692773e-01  0.125637129620766
    ## 52      UkbDep ~~      UkbDep -7.390970e-07   0.13076256181067
    ## 27     CommAnh ~~     CommAnh  5.619943e-02 0.0449855578858504
    ## 51      UkbAnh ~~      UkbAnh  1.874083e-01  0.084328695176069
    ## 33   CommGuilt ~~   CommGuilt  5.288532e-01  0.166908697358572
    ## 31     CommDep ~~     CommDep  1.176714e-01 0.0538295284441127
    ## 26     ClinSui ~~     ClinSui  6.625437e-01  0.323243991766738
    ## 36     CommSui ~~     CommSui  6.092822e-01  0.232200955414748
    ## 24  ClinSleDec ~~  ClinSleDec  9.788317e-01   0.58137590113302
    ## 34  CommSleDec ~~  CommSleDec  4.691040e-01  0.303530642565606
    ## 23 ClinMotoInc ~~ ClinMotoInc  9.419939e-01  0.390564647439688
    ## 20  ClinAppDec ~~  ClinAppDec  9.971666e-01  0.218873390668526
    ## 28  CommAppDec ~~  CommAppDec  9.525499e-01  0.237112116267575
    ## 30    CommConc ~~    CommConc  1.086020e-01  0.298908961962882
    ## 21  ClinAppInc ~~  ClinAppInc  9.718354e-01  0.317174255592681
    ## 29  CommAppInc ~~  CommAppInc  7.773862e-01  0.153547996425434
    ## 25  ClinSleInc ~~  ClinSleInc  9.400301e-01  0.668262912115568
    ## 35  CommSleInc ~~  CommSleInc  5.978323e-01  0.253859853153597
    ## 22 ClinMotoDec ~~ ClinMotoDec  9.041202e-01   1.08506693959082
    ## 32   CommFatig ~~   CommFatig  2.670191e-01  0.316381516471251
    ## 11      AFFECT ~~        SOMA  8.238282e-01  0.126358121555781
    ## 9       AFFECT ~~         ATY  7.595830e-01  0.125443591232808
    ## 48        SOMA ~~         ATY  6.589604e-01  0.159014402565926
    ## 8       AFFECT ~~      AFFECT  1.000000e+00                   
    ## 50        SOMA ~~        SOMA  1.000000e+00                   
    ## 18         ATY ~~         ATY  1.000000e+00                   
    ## 41        GATE ~~        GATE  1.000000e+00                   
    ## 10      AFFECT ~~        GATE  0.000000e+00                   
    ## 49        SOMA ~~        GATE  0.000000e+00                   
    ## 19         ATY ~~        GATE  0.000000e+00

``` r
     #filter(lhs == 'MDD') %>%
     #mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     #print(digits=2)
```

``` r
psych_soma_aty.model <- "
PSYCH =~ NA*CommAnh + UkbAnh + CommGuilt + CommDep + UkbDep + ClinSui + CommSui + CommConc + CommFatig
SOMA =~ NA*ClinSleDec + CommSleDec + ClinMotoInc + ClinAppDec + CommAppDec 
ATY =~ NA*ClinAppInc + CommAppInc + ClinSleInc + CommSleInc + ClinMotoDec
GATE =~ NA*CommDep + CommAnh + UkbDep + UkbAnh

PSYCH ~~ 1*PSYCH
SOMA ~~ 1*SOMA
ATY ~~ 1*ATY
GATE ~~ 1*GATE
GATE ~~ 0*PSYCH + 0*SOMA + 0*ATY

u1 > 0
UkbDep ~~ u1*UkbDep

co4a > 0
CommSleDec ~~ co4a*CommSleDec
"

psych_soma_aty.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=psych_soma_aty.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  12.436 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0465162451355015 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65918051782736 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from the
    ## model. If you are going to run a multivariate GWAS we strongly recommend setting
    ## the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## psych_soma_aty.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
psych_soma_aty.fit$modelfit
```

    ##       chisq  df       p_chisq      AIC      CFI      SRMR
    ## df 1066.514 145 3.891465e-140 1156.514 0.961614 0.1541632

``` r
psych_soma_aty.fit$results[c(1,2,3,6,7)]# %>%
```

    ##            lhs op         rhs  STD_Genotype    STD_Genotype_SE
    ## 31       PSYCH =~     CommAnh  8.145432e-01 0.0652084481250952
    ## 37       PSYCH =~      UkbAnh  6.854777e-01 0.0888774920587646
    ## 35       PSYCH =~   CommGuilt  6.864946e-01 0.0898784911543086
    ## 33       PSYCH =~     CommDep  6.876657e-01 0.0681474934041511
    ## 38       PSYCH =~      UkbDep  5.449505e-01 0.0940887126295368
    ## 30       PSYCH =~     ClinSui  5.853959e-01  0.106832899990355
    ## 36       PSYCH =~     CommSui  6.236698e-01 0.0992882333021443
    ## 32       PSYCH =~    CommConc  8.082793e-01 0.0992587334133928
    ## 34       PSYCH =~   CommFatig  6.885592e-01  0.101421574149024
    ## 45        SOMA =~  ClinSleDec  2.596966e-01  0.194391868876785
    ## 47        SOMA =~  CommSleDec  1.062567e+00  0.272207138823246
    ## 44        SOMA =~ ClinMotoInc  3.455909e-01  0.170287211196701
    ## 43        SOMA =~  ClinAppDec -1.840620e-03   0.12538156108942
    ## 46        SOMA =~  CommAppDec  3.218455e-01  0.125513345778186
    ## 1          ATY =~  ClinAppInc  2.158547e-01  0.155388808895223
    ## 4          ATY =~  CommAppInc  5.475258e-01  0.127165595142391
    ## 3          ATY =~  ClinSleInc  2.703933e-01  0.212199815511966
    ## 5          ATY =~  CommSleInc  7.253026e-01  0.177103896575548
    ## 2          ATY =~ ClinMotoDec  3.604125e-01  0.245637912916686
    ## 26        GATE =~     CommDep  6.400643e-01 0.0937342780070344
    ## 25        GATE =~     CommAnh  5.258980e-01  0.092548850306844
    ## 28        GATE =~      UkbDep  8.455379e-01  0.114968303695522
    ## 27        GATE =~      UkbAnh  5.855744e-01  0.113736863813743
    ## 52      UkbDep ~~      UkbDep  5.937596e-08  0.123444530989431
    ## 22  CommSleDec ~~  CommSleDec -6.122706e-08  0.582567517165408
    ## 15     CommAnh ~~     CommAnh  5.995085e-02 0.0436612788808545
    ## 51      UkbAnh ~~      UkbAnh  1.872225e-01 0.0843500587628739
    ## 21   CommGuilt ~~   CommGuilt  5.287244e-01  0.166255984771993
    ## 19     CommDep ~~     CommDep  1.174328e-01 0.0530455925311709
    ## 14     ClinSui ~~     ClinSui  6.573087e-01   0.32312760930521
    ## 24     CommSui ~~     CommSui  6.110349e-01   0.23256517535315
    ## 18    CommConc ~~    CommConc  3.466811e-01  0.242531933387732
    ## 20   CommFatig ~~   CommFatig  5.258826e-01  0.296003885048891
    ## 12  ClinSleDec ~~  ClinSleDec  9.326470e-01  0.574777454427649
    ## 11 ClinMotoInc ~~ ClinMotoInc  8.805845e-01  0.391246198104237
    ## 8   ClinAppDec ~~  ClinAppDec  9.999966e-01  0.218584034304953
    ## 16  CommAppDec ~~  CommAppDec  8.964190e-01  0.238835710644493
    ## 9   ClinAppInc ~~  ClinAppInc  9.534064e-01  0.318945770842939
    ## 17  CommAppInc ~~  CommAppInc  7.002154e-01  0.181836055754543
    ## 13  ClinSleInc ~~  ClinSleInc  9.268920e-01  0.672566692598711
    ## 23  CommSleInc ~~  CommSleInc  4.739380e-01   0.29426593072326
    ## 10 ClinMotoDec ~~ ClinMotoDec  8.700948e-01    1.0942972447657
    ## 42       PSYCH ~~        SOMA  5.721856e-01  0.162600280865084
    ## 39       PSYCH ~~         ATY  6.738497e-01  0.139221554021235
    ## 48        SOMA ~~         ATY  3.143211e-01  0.209969123902117
    ## 41       PSYCH ~~       PSYCH  1.000000e+00                   
    ## 50        SOMA ~~        SOMA  1.000000e+00                   
    ## 6          ATY ~~         ATY  1.000000e+00                   
    ## 29        GATE ~~        GATE  1.000000e+00                   
    ## 40       PSYCH ~~        GATE  0.000000e+00                   
    ## 49        SOMA ~~        GATE  0.000000e+00                   
    ## 7          ATY ~~        GATE  0.000000e+00

``` r
     #filter(lhs == 'MDD') %>%
     #mutate(STD_Genotype_SE=as.numeric(STD_Genotype_SE)) %>%
     #print(digits=2)
```

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
     "1b"=list(name="Common (App)", model=commonfactor_app.fit),
     "1c"=list(name="Common (App, sub)", model=commonfactor_app_sub.fit),
     "1d"=list(name="Common (Sle)", model=commonfactor_sle.fit),
     "1e"=list(name="Common (gating)", model=gate.fit),
     "1f"=list(name="Ascertainment", model=clin_comm.fit),
     "1g"=list(name="Measurement", model=measure.fit),
     "1h"=list(name="Divergence", model=divergence.fit),
     "2a"=list(name="Psych-Somatic", model=psych_soma.fit),
     "2b"=list(name="Psych-Neuroveg", model=psych_veg.fit),
     "2c"=list(name="Affect-Neuroveg", model=affect_neuroveg.fit),
     "2d"=list(name="Melancholic/Atypical", model=mel_aty.fit),
     "2e"=list(names="Affect-Somatic-Atypical", model=mel_aty_psych_soma.fit),
     "2f"=list(names="Psych-Somatic-Atypical", model=psych_soma_aty.fit),
     "3"=list(name="Cog-Mood-Neuroveg", model=cog_mood_neuroveg.fit),
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

| Model | Name                    | chisq |  df | p\_chisq |  AIC |   CFI |  SRMR |  dAIC |
| :---- | :---------------------- | ----: | --: | -------: | ---: | ----: | ----: | ----: |
| 1a    | Common                  |  1340 | 152 |        0 | 1410 | 0.951 | 0.175 | 279.0 |
| 1b    | Common (App)            |  1190 | 148 |        0 | 1280 | 0.957 | 0.158 | 142.0 |
| 1c    | Common (App, sub)       |  1130 | 151 |        0 | 1210 | 0.959 | 0.166 |  75.7 |
| 1d    | Common (Sle)            |  1300 | 148 |        0 | 1390 | 0.952 | 0.170 | 252.0 |
| 1e    | Common (gating)         |  1150 | 148 |        0 | 1230 | 0.958 | 0.163 | 101.0 |
| 1f    | Ascertainment           |  1330 | 151 |        0 | 1410 | 0.951 | 0.175 | 271.0 |
| 1g    | Measurement             |  1480 | 149 |        0 | 1560 | 0.945 | 0.160 | 428.0 |
| 1h    | Divergence              |  1320 | 151 |        0 | 1400 | 0.951 | 0.174 | 269.0 |
| 2a    | Psych-Somatic           |  1140 | 147 |        0 | 1220 | 0.959 | 0.163 |  88.5 |
| 2b    | Psych-Neuroveg          |  1080 | 147 |        0 | 1170 | 0.961 | 0.161 |  34.0 |
| 2c    | Affect-Neuroveg         |  1150 | 147 |        0 | 1230 | 0.958 | 0.161 | 101.0 |
| 2d    | Melancholic/Atypical    |  1040 | 145 |        0 | 1130 | 0.963 | 0.159 |   0.0 |
| 2e    | Affect-Somatic-Atypical |  1170 | 145 |        0 | 1260 | 0.957 | 0.157 | 123.0 |
| 2f    | Psych-Somatic-Atypical  |  1070 | 145 |        0 | 1160 | 0.962 | 0.154 |  22.5 |
| 3     | Cog-Mood-Neuroveg       |  1090 | 144 |        0 | 1180 | 0.961 | 0.162 |  46.3 |
| 4     | Cog-Mood-Mel-Aty        |  1270 | 141 |        0 | 1370 | 0.953 | 0.152 | 236.0 |

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
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec 
    ##       0.082       0.069       0.005       0.831       0.895       0.968 
    ##     ClinSui     CommDep     CommAnh  CommAppDec  CommAppInc  CommSleDec 
    ##       0.798       0.112       0.050       0.869       0.290       0.440 
    ##  CommSleInc   CommFatig   CommGuilt    CommConc     CommSui      UkbDep 
    ##       0.757       0.646       0.585       0.583       0.712       0.122 
    ##      UkbAnh 
    ##       0.168 
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
