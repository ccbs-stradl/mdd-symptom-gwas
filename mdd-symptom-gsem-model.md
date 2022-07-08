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
    ## platform       x86_64-generic-linux-gnu    
    ## arch           x86_64                      
    ## os             linux-gnu                   
    ## system         x86_64, linux-gnu           
    ## status                                     
    ## major          4                           
    ## minor          2.1                         
    ## year           2022                        
    ## month          06                          
    ## day            23                          
    ## svn rev        82513                       
    ## language       R                           
    ## version.string R version 4.2.1 (2022-06-23)
    ## nickname       Funny-Looking Kid

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
    ## ── Column specification ────────────────────────────────────────────────
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
    ## ── Column specification ────────────────────────────────────────────────
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
    ##  1 MDD1      Dep          Mood        Depressed mood most of the day, nearly ev…
    ##  2 MDD2      Anh          Interest    Markedly diminished interest or pleasure …
    ##  3 MDD3      App          Weight⇅     Significant change in weight or appetite  
    ##  4 MDD3a     AppDec       Weight⇊     Significant weight loss or decrease in ap…
    ##  5 MDD3b     AppInc       Weight⇈     Significant weight gain or increase in ap…
    ##  6 MDD4      Sle          Sleep⇅      Sleeping too much or not sleeping enough  
    ##  7 MDD4a     SleDec       Sleep⇊      Insomnia nearly every day                 
    ##  8 MDD4b     SleInc       Sleep⇈      Hypersomnia nearly every day              
    ##  9 MDD5      Moto         Motor⇅      Changes in speed/amount of moving or spea…
    ## 10 MDD5a     MotoInc      Motor⇈      Psychomotor agitation nearly every day    
    ## 11 MDD5b     MotoDec      Motor⇊      Psychomotor slowing nearly every day      
    ## 12 MDD6      Fatig        Fatigue     Fatigue or loss of energy nearly every day
    ## 13 MDD7      Guilt        Guilt       Feelings of worthlessness or excessive or…
    ## 14 MDD8      Conc         Concentrate Diminished ability to think or concentrat…
    ## 15 MDD9      Sui          Suicidality Recurrent thoughts of death or suicide or…

# GenomicSEM covariance structure

``` r
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 24 Columns: 9
    ## ── Column specification ────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (5): cohorts, symptom, sumstats, filename, trait_name
    ## dbl (4): Nca, Nco, samp_prev, pop_prev
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

Rename samples: AGDS/PGC is the **Clin**ical sample (`Clin`) and
ALSPAC/UKB is the **Pop**ulation sample (`Pop`); and rename symptoms
numbers (`MDD1`, `MDD2`) to abbreviations (`Dep`, `Anh`)

``` r
cohorts_sample_symptoms <-
sumstats_prevs %>%
left_join(dsm_mdd_symptoms_labels, by=c('symptom'='ref')) %>%
select(cohorts, symptom, trait_name, abbv) %>%
mutate(Sample=case_when(cohorts %in% 'AGDS_PGC' ~ 'Clin',
                        cohorts %in% 'ALSPAC_UKB' ~ 'Pop',
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

    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc     ClinSui 
    ## 0.107310576 0.042492134 0.025429781 0.003850035 0.036308706 0.103892914 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ## 0.082198484 0.087343514 0.033528770 0.054402455 0.032178877 0.024556756 
    ##  PopMotoInc  PopMotoDec    PopFatig    PopGuilt     PopConc      PopSui 
    ## 0.052444420 0.019177797 0.057539570 0.060012455 0.056887041 0.033166683

## Common factor

Common factor across symptoms from both cohorts

``` r
commonfactor.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui + PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
"

commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    2.72 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0313666288930831 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.49543224219422 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model): A difference greater than .025 was observed pre- and post-
    ## smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
commonfactor.fit$modelfit
```

    ##       chisq  df     p_chisq      AIC       CFI      SRMR
    ## df 306.5241 104 7.57386e-22 370.5241 0.9546511 0.1814247

``` r
commonfactor.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE
    ## 1           A1 =~  ClinAppDec   0.08104673  0.088692257809618
    ## 2           A1 =~  ClinAppInc  -0.19252247  0.110239144890949
    ## 4           A1 =~  ClinSleDec  -0.05910385  0.124424048904914
    ## 5           A1 =~  ClinSleInc  -0.14288769  0.165451392570719
    ## 3           A1 =~ ClinMotoInc  -0.18967181  0.115683977419173
    ## 6           A1 =~     ClinSui  -0.58557524 0.0972352466438474
    ## 11          A1 =~      PopDep  -0.83414952 0.0503724828065146
    ## 7           A1 =~      PopAnh  -0.94786091 0.0445806621410705
    ## 8           A1 =~   PopAppDec  -0.17469342 0.0816360764655869
    ## 9           A1 =~   PopAppInc  -0.39226944 0.0737935213630935
    ## 14          A1 =~   PopSleDec  -0.57641902 0.0938473782944009
    ## 15          A1 =~   PopSleInc  -0.45708427 0.0981954915139065
    ## 12          A1 =~    PopFatig  -0.66597064 0.0912349191598555
    ## 13          A1 =~    PopGuilt  -0.64621351 0.0776426242285542
    ## 10          A1 =~     PopConc  -0.70369473 0.0923057396382578
    ## 16          A1 =~      PopSui  -0.58272279 0.0955634257144546
    ## 18  ClinAppDec ~~  ClinAppDec   0.99343314  0.225297847010729
    ## 19  ClinAppInc ~~  ClinAppInc   0.96293253   0.32882932806548
    ## 21  ClinSleDec ~~  ClinSleDec   0.99650374  0.493519574978079
    ## 22  ClinSleInc ~~  ClinSleInc   0.97956091  0.845479581271306
    ## 20 ClinMotoInc ~~ ClinMotoInc   0.96402296  0.432698456246632
    ## 23     ClinSui ~~     ClinSui   0.65710096  0.306087517627822
    ## 28      PopDep ~~      PopDep   0.30419491 0.0746603550465037
    ## 24      PopAnh ~~      PopAnh   0.10155899 0.0623210509107302
    ## 25   PopAppDec ~~   PopAppDec   0.96948351  0.222924875829309
    ## 26   PopAppInc ~~   PopAppInc   0.84612545  0.146051589076692
    ## 31   PopSleDec ~~   PopSleDec   0.66774026  0.287954827295504
    ## 32   PopSleInc ~~   PopSleInc   0.79106882   0.25271446550866
    ## 29    PopFatig ~~    PopFatig   0.55648340  0.308678015647095
    ## 30    PopGuilt ~~    PopGuilt   0.58240650  0.159257041450059
    ## 27     PopConc ~~     PopConc   0.50481509   0.24574581842838
    ## 33      PopSui ~~      PopSui   0.66043312  0.237317281683625
    ## 17          A1 ~~          A1   1.00000000

Correlation among directional symptoms

``` r
commonfactor_dir.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui + PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
ClinAppDec ~~ ClinAppInc
ClinSleDec ~~ ClinSleInc
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"

commonfactor_dir.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor_dir.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   3.238 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0313666288930831 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.49543224219422 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_dir.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## commonfactor_dir.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
commonfactor_dir.fit$modelfit
```

    ##      chisq df      p_chisq     AIC       CFI      SRMR
    ## df 216.872 99 8.449537e-11 290.872 0.9736063 0.1659054

``` r
commonfactor_dir.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppDec   0.09265454 0.0948395328152933 3.285888e-01
    ## 2           A1 =~  ClinAppInc  -0.22677683  0.117690768029817 5.399351e-02
    ## 4           A1 =~  ClinSleDec  -0.06057619  0.130203753541332 6.417609e-01
    ## 5           A1 =~  ClinSleInc  -0.16774807  0.176913842508782 3.430341e-01
    ## 3           A1 =~ ClinMotoInc  -0.20997830  0.124342135808605 9.127476e-02
    ## 6           A1 =~     ClinSui  -0.64306000  0.105472338773873 1.081101e-09
    ## 11          A1 =~      PopDep  -0.67325982 0.0687442541985537 1.198437e-22
    ## 7           A1 =~      PopAnh  -0.79701839 0.0632011142836478 1.840485e-36
    ## 8           A1 =~   PopAppDec  -0.21022580 0.0880189808287135 1.692108e-02
    ## 9           A1 =~   PopAppInc  -0.43634476 0.0795718865718979 4.166142e-08
    ## 14          A1 =~   PopSleDec  -0.65125732   0.10420883118702 4.116695e-10
    ## 15          A1 =~   PopSleInc  -0.50348948  0.105455610498766 1.802301e-06
    ## 12          A1 =~    PopFatig  -0.72860552 0.0970127705283394 5.894356e-14
    ## 13          A1 =~    PopGuilt  -0.70134977  0.083622518456887 4.982971e-17
    ## 10          A1 =~     PopConc  -0.77030805 0.0991520230487558 7.913442e-15
    ## 16          A1 =~      PopSui  -0.61151914  0.101883885413634 1.947601e-09
    ## 19  ClinAppDec ~~  ClinAppInc  -0.38012841  0.212408295464654 7.351603e-02
    ## 23  ClinSleDec ~~  ClinSleInc   0.47376696  0.475346748967055 3.189214e-01
    ## 31      PopDep ~~      PopAnh   0.37004083 0.0903453762454572 4.206240e-05
    ## 28   PopAppDec ~~   PopAppInc  -0.23925516  0.127118166487949 5.981562e-02
    ## 36   PopSleDec ~~   PopSleInc  -0.30564259  0.197904872332296 1.224940e-01
    ## 18  ClinAppDec ~~  ClinAppDec   0.99141662  0.225845329969803 1.134577e-05
    ## 20  ClinAppInc ~~  ClinAppInc   0.94857252  0.330168144707281 4.065983e-03
    ## 22  ClinSleDec ~~  ClinSleDec   0.99633076  0.493412479192737 4.345968e-02
    ## 24  ClinSleInc ~~  ClinSleInc   0.97185771   0.84548331787448 2.503607e-01
    ## 21 ClinMotoInc ~~ ClinMotoInc   0.95590927  0.432530905655498 2.710254e-02
    ## 25     ClinSui ~~     ClinSui   0.58647352   0.31461807063136 6.230947e-02
    ## 32      PopDep ~~      PopDep   0.54672108  0.107892453543588 4.035401e-07
    ## 26      PopAnh ~~      PopAnh   0.36476172 0.0959196267144642 1.430786e-04
    ## 27   PopAppDec ~~   PopAppDec   0.95580294  0.221809603534033 1.639023e-05
    ## 29   PopAppInc ~~   PopAppInc   0.80960186  0.147627628539881 4.156108e-08
    ## 35   PopSleDec ~~   PopSleDec   0.57586408  0.291035207960705 4.785253e-02
    ## 37   PopSleInc ~~   PopSleInc   0.74649969  0.251514376520577 2.997269e-03
    ## 33    PopFatig ~~    PopFatig   0.46913381  0.313326320671658 1.343231e-01
    ## 34    PopGuilt ~~    PopGuilt   0.50811009  0.164432246708526 2.001070e-03
    ## 30     PopConc ~~     PopConc   0.40662644  0.247681698578753 1.006464e-01
    ## 38      PopSui ~~      PopSui   0.62604460  0.240365292942291 9.199385e-03
    ## 17          A1 ~~          A1   1.00000000                              NA

Ascertainment-specific factors

``` r
clin_pop.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui 
A2 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
A2 ~~ 1*A2
a12 > -1
a12 < 1
A1 ~~ a12*A2
ClinAppDec ~~ ClinAppInc
ClinSleDec ~~ ClinSleInc
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"

clin_pop.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_pop.model)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##            lhs op         rhs Unstandardized_Estimate          SE
    ## 1           A1 =~  ClinAppDec              0.03121902 0.034606253
    ## 2           A1 =~  ClinAppInc             -0.05949859 0.042159699
    ## 3           A1 =~  ClinSleDec             -0.01443643 0.031448316
    ## 4           A1 =~  ClinSleInc             -0.02560026 0.029425420
    ## 5           A1 =~ ClinMotoInc             -0.04115131 0.030102032
    ## 6           A1 =~     ClinSui             -0.21414125 0.115454739
    ## 7           A2 =~      PopDep              0.19504058 0.020012711
    ## 8           A2 =~      PopAnh              0.23696295 0.018875802
    ## 9           A2 =~   PopAppDec              0.04030721 0.016876080
    ## 10          A2 =~   PopAppInc              0.10647489 0.019416395
    ## 11          A2 =~   PopSleDec              0.12566958 0.020107594
    ## 12          A2 =~   PopSleInc              0.08547669 0.017903221
    ## 13          A2 =~    PopFatig              0.18364888 0.024453087
    ## 14          A2 =~    PopGuilt              0.17238408 0.020551493
    ## 15          A2 =~     PopConc              0.19381784 0.024946164
    ## 16          A2 =~      PopSui              0.11610901 0.019343881
    ## 19          A1 ~~          A2             -1.00000016 0.555106904
    ## 20  ClinAppDec ~~  ClinAppInc             -0.03360411 0.018816184
    ## 21  ClinSleDec ~~  ClinSleInc              0.01723112 0.017261898
    ## 22      PopDep ~~      PopAnh              0.03187171 0.007797375
    ## 23   PopAppDec ~~   PopAppInc             -0.01119373 0.005947319
    ## 24   PopSleDec ~~   PopSleInc             -0.01001264 0.006483151
    ## 25  ClinAppDec ~~  ClinAppDec              0.11255455 0.025649326
    ## 26  ClinAppInc ~~  ClinAppInc              0.06529586 0.022841906
    ## 27  ClinSleDec ~~  ClinSleDec              0.05658800 0.028034334
    ## 28  ClinSleInc ~~  ClinSleInc              0.02263501 0.019696214
    ## 29 ClinMotoInc ~~ ClinMotoInc              0.03671441 0.016602446
    ## 30     ClinSui ~~     ClinSui              0.06503423 0.055841269
    ## 31      PopDep ~~      PopDep              0.04588295 0.009068575
    ## 32      PopAnh ~~      PopAnh              0.03224294 0.008486012
    ## 33   PopAppDec ~~   PopAppDec              0.03513673 0.008154027
    ## 34   PopAppInc ~~   PopAppInc              0.04820657 0.008790334
    ## 35   PopSleDec ~~   PopSleDec              0.02144251 0.010836463
    ## 36   PopSleInc ~~   PopSleInc              0.02151517 0.007249024
    ## 37    PopFatig ~~    PopFatig              0.02980495 0.019906352
    ## 38    PopGuilt ~~    PopGuilt              0.03069589 0.009931404
    ## 39     PopConc ~~     PopConc              0.02574266 0.015680155
    ## 40      PopSui ~~      PopSui              0.02256925 0.008665397

``` r
clin_pop.fit$modelfit
```

    ## NULL

``` r
clin_pop.fit$results[c(1,2,3,6,7,9)]
```

    ## NULL

``` r
clin_pop_bif.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui 
A2 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui + PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
A2 ~~ 1*A2
A ~~ 1*A
A ~~ 0*A1 + 0*A2
A1 ~~ 0*A2
"

clin_pop_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_pop_bif.model)
```

    ## [1] "Running primary model"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = clin_pop_bif.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

    ##            lhs op         rhs Unstandardized_Estimate
    ## 1           A1 =~  ClinAppDec           -3.907871e-03
    ## 2           A1 =~  ClinAppInc           -2.591313e-03
    ## 3           A1 =~  ClinSleDec           -1.067074e+01
    ## 4           A1 =~  ClinSleInc           -1.579888e-03
    ## 5           A1 =~ ClinMotoInc           -6.898743e-04
    ## 6           A1 =~     ClinSui            1.828071e-03
    ## 7           A2 =~      PopDep            2.521470e-01
    ## 8           A2 =~      PopAnh            2.412534e-01
    ## 9           A2 =~   PopAppDec            2.524625e-02
    ## 10          A2 =~   PopAppInc            2.420650e-02
    ## 11          A2 =~   PopSleDec            4.382233e-02
    ## 12          A2 =~   PopSleInc            5.697106e-02
    ## 13          A2 =~    PopFatig            6.318682e-02
    ## 14          A2 =~    PopGuilt            6.750016e-02
    ## 15          A2 =~     PopConc            7.159646e-02
    ## 16          A2 =~      PopSui            8.203891e-02
    ## 17           A =~  ClinAppDec            5.649715e-02
    ## 18           A =~  ClinAppInc           -9.159902e-02
    ## 19           A =~  ClinSleDec           -2.197105e-02
    ## 20           A =~  ClinSleInc           -3.412306e-02
    ## 21           A =~ ClinMotoInc           -5.198950e-02
    ## 22           A =~     ClinSui           -2.722056e-01
    ## 23           A =~      PopDep           -1.096020e-01
    ## 24           A =~      PopAnh           -1.599323e-01
    ## 25           A =~   PopAppDec           -2.006037e-02
    ## 26           A =~   PopAppInc           -1.151724e-01
    ## 27           A =~   PopSleDec           -1.134208e-01
    ## 28           A =~   PopSleInc           -5.065070e-02
    ## 29           A =~    PopFatig           -1.749063e-01
    ## 30           A =~    PopGuilt           -1.567124e-01
    ## 31           A =~     PopConc           -1.795741e-01
    ## 32           A =~      PopSui           -7.195806e-02
    ## 39  ClinAppDec ~~  ClinAppDec            1.103215e-01
    ## 40  ClinAppInc ~~  ClinAppInc            6.043693e-02
    ## 41  ClinSleDec ~~  ClinSleDec           -1.138085e+02
    ## 42  ClinSleInc ~~  ClinSleInc            2.212343e-02
    ## 43 ClinMotoInc ~~ ClinMotoInc            3.570351e-02
    ## 44     ClinSui ~~     ClinSui            3.679150e-02
    ## 45      PopDep ~~      PopDep            8.333174e-03
    ## 46      PopAnh ~~      PopAnh            4.612999e-03
    ## 47   PopAppDec ~~   PopAppDec            3.572181e-02
    ## 48   PopAppInc ~~   PopAppInc            4.569279e-02
    ## 49   PopSleDec ~~   PopSleDec            2.245043e-02
    ## 50   PopSleInc ~~   PopSleInc            2.301031e-02
    ## 51    PopFatig ~~    PopFatig            2.894753e-02
    ## 52    PopGuilt ~~    PopGuilt            3.129708e-02
    ## 53     PopConc ~~     PopConc            2.593548e-02
    ## 54      PopSui ~~      PopSui            2.414219e-02

``` r
clin_pop_bif.fit$modelfit
```

    ## NULL

``` r
clin_pop_bif.fit$results[c(1,2,3,6,7,9)]
```

    ## NULL

## ADGS-PGC (Clinical)

### Common factor

Common factor model. Allow residual negative correlation between
directional symptoms

``` r
clin_commonfactor.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui
A1 ~~ 1*A1
"
clin_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_commonfactor.model)
```

    ## [1] "Running primary model"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = clin_commonfactor.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

    ##            lhs op         rhs Unstandardized_Estimate
    ## 1           A1 =~  ClinAppDec            2.247074e-03
    ## 2           A1 =~  ClinAppInc            1.868877e-03
    ## 3           A1 =~  ClinSleDec            1.660103e+01
    ## 4           A1 =~  ClinSleInc            1.135322e-03
    ## 5           A1 =~ ClinMotoInc            5.461996e-04
    ## 6           A1 =~     ClinSui           -9.567093e-04
    ## 8   ClinAppDec ~~  ClinAppDec            1.147987e-01
    ## 9   ClinAppInc ~~  ClinAppInc            6.245045e-02
    ## 10  ClinSleDec ~~  ClinSleDec           -2.755377e+02
    ## 11  ClinSleInc ~~  ClinSleInc            1.211893e-02
    ## 12 ClinMotoInc ~~ ClinMotoInc            3.680382e-02
    ## 13     ClinSui ~~     ClinSui            1.069355e-01

``` r
clin_commonfactor.fit$modelfit
```

    ## NULL

``` r
clin_commonfactor.fit$results[c(1,2,3,6,7, 9)]
```

    ## NULL

Add residual covariance between appetite symptoms.

``` r
clin_commonfactor_app.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui
A1 ~~ 1*A1
c3a3b > -1
ClinAppDec ~~ c3a3b*ClinAppInc
"
clin_commonfactor_app.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_commonfactor_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   6.305 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0311960347054173 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.41818584418939 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_commonfactor_app.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_commonfactor_app.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
clin_commonfactor_app.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 2.325305  8 0.9693801 28.32531   1 0.1438518

``` r
clin_commonfactor_app.fit$results[c(1,2,3,6,7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE     p_value
    ## 1           A1 =~  ClinAppDec    0.6634310 0.428690400644662 0.121725888
    ## 2           A1 =~  ClinAppInc    0.7278282 0.484179869694426 0.132784610
    ## 4           A1 =~  ClinSleDec    0.7267902 0.431683423039736 0.092257501
    ## 5           A1 =~  ClinSleInc    0.6208064 0.487198098422708 0.202581095
    ## 3           A1 =~ ClinMotoInc    0.2864365  0.23363640265699 0.220202308
    ## 6           A1 =~     ClinSui    0.0334615  0.17827678796196 0.851114382
    ## 9   ClinAppDec ~~  ClinAppInc   -0.8815747  0.55995562148932 0.115406962
    ## 8   ClinAppDec ~~  ClinAppDec    0.5598592  0.56769604800426 0.324052438
    ## 10  ClinAppInc ~~  ClinAppInc    0.4702662 0.752007458365769 0.531756151
    ## 12  ClinSleDec ~~  ClinSleDec    0.4717756 0.795763911922935 0.553266114
    ## 13  ClinSleInc ~~  ClinSleInc    0.6145997   1.7186638468124 0.720638660
    ## 11 ClinMotoInc ~~ ClinMotoInc    0.9179545 0.460047958933582 0.046005368
    ## 14     ClinSui ~~     ClinSui    0.9988805 0.298668758746582 0.000824483
    ## 7           A1 ~~          A1    1.0000000                            NA

Add residual covariance between sleep symptoms.

``` r
clin_commonfactor_sle.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc + ClinSui
A1 ~~ 1*A1
ClinSleDec ~~ c3a3b*ClinSleInc
"
clin_commonfactor_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_commonfactor_sle.model)
```

    ## [1] "Running primary model"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = clin_commonfactor_sle.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

    ##            lhs op         rhs Unstandardized_Estimate
    ## 1           A1 =~  ClinAppDec            1.985407e+01
    ## 2           A1 =~  ClinAppInc           -1.700226e-03
    ## 3           A1 =~  ClinSleDec            1.878797e-03
    ## 4           A1 =~  ClinSleInc            8.764709e-04
    ## 5           A1 =~ ClinMotoInc            6.654713e-04
    ## 6           A1 =~     ClinSui           -1.213248e-04
    ## 8   ClinSleDec ~~  ClinSleInc            1.884268e-02
    ## 9   ClinAppDec ~~  ClinAppDec           -3.940692e+02
    ## 10  ClinAppInc ~~  ClinAppInc            6.245104e-02
    ## 11  ClinSleDec ~~  ClinSleDec            5.662229e-02
    ## 12  ClinSleInc ~~  ClinSleInc            1.211952e-02
    ## 13 ClinMotoInc ~~ ClinMotoInc            3.680377e-02
    ## 14     ClinSui ~~     ClinSui            1.069364e-01

``` r
clin_commonfactor_sle.fit$modelfit
```

    ## NULL

``` r
clin_commonfactor_sle.fit$results[c(1,2,3,6,7, 9)]
```

    ## NULL

### Somatic factor

``` r
clin_soma.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc
A2 =~ ClinSui
A1 ~~ 1*A1
ClinAppDec ~~ ClinAppInc
ClinSui ~~ 0*ClinSui
"
clin_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.239 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0311960347054173 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.41818584418939 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_soma.model): A difference greater than .025 was observed pre- and post-
    ## smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_soma.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_soma.fit$modelfit
```

    ##       chisq df p_chisq      AIC CFI      SRMR
    ## df 2.325307  8 0.96938 28.32531   1 0.1438516

``` r
clin_soma.fit$results[c(1,2,3,6,7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppDec   0.66343100 0.428690348129595 0.1217242425
    ## 2           A1 =~  ClinAppInc   0.72782809 0.484179818864062 0.1327824813
    ## 4           A1 =~  ClinSleDec   0.72679013 0.431683341319103 0.0922552892
    ## 5           A1 =~  ClinSleInc   0.62080645 0.487198132037538 0.2025786754
    ## 3           A1 =~ ClinMotoInc   0.28643659 0.233636432496519 0.2202002892
    ## 11  ClinAppDec ~~  ClinAppInc  -0.88157467 0.559955481732302 0.1154030513
    ## 10  ClinAppDec ~~  ClinAppDec   0.55985943 0.567695949545366 0.3240382605
    ## 12  ClinAppInc ~~  ClinAppInc   0.47026617 0.752007319664471 0.5317463160
    ## 14  ClinSleDec ~~  ClinSleDec   0.47177613 0.795763739738203 0.5532714767
    ## 15  ClinSleInc ~~  ClinSleInc   0.61459918  1.71866387203324 0.7206414568
    ## 13 ClinMotoInc ~~ ClinMotoInc   0.91795401 0.460047975966117 0.0460056420
    ## 9           A2 ~~          A2   1.00000000 0.298426084011083 0.0008054545
    ## 7           A1 ~~          A2   0.03346141 0.178276802750801 0.8511189551
    ## 8           A2 =~     ClinSui   1.00000000                             NA
    ## 6           A1 ~~          A1   1.00000000                             NA
    ## 16     ClinSui ~~     ClinSui   0.00000000                             NA

``` r
clin_soma_app.model <- "
A1 =~ NA*ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc
A2 =~ ClinAppDec
A3 =~ ClinSui
A1 ~~ 1*A1
ClinAppDec ~~ ClinAppInc
ClinAppDec ~~ 0*ClinAppDec
ClinSui ~~ 0*ClinSui
c4a > 0.001
ClinSleDec ~~ c4a*ClinSleDec
"

clin_soma_app.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_soma_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   6.828 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0311960347054173 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.41818584418939 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_soma_app.model): A difference greater than .025 was observed pre- and post-
    ## smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_soma_app.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_soma_app.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 2.313668  7 0.9404629 30.31367   1 0.1485743

``` r
clin_soma_app.fit$results[c(1,2,3,6,7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppInc   0.85629294 0.619378484978697 1.668175e-01
    ## 3           A1 =~  ClinSleDec   0.61593755 0.416067244728193 1.387722e-01
    ## 4           A1 =~  ClinSleInc   0.56104071 0.479390846897817 2.418726e-01
    ## 2           A1 =~ ClinMotoInc   0.27125530 0.233721508225008 2.458082e-01
    ## 14  ClinAppInc ~~  ClinAppDec  -1.05527303 0.835551112138854 2.066015e-01
    ## 17  ClinSleDec ~~  ClinSleDec   0.62062217 0.699745346532365 3.751196e-01
    ## 15  ClinAppInc ~~  ClinAppInc   0.26676201  1.08509983928275 8.058049e-01
    ## 18  ClinSleInc ~~  ClinSleInc   0.68523356  1.69476885259179 6.859753e-01
    ## 16 ClinMotoInc ~~ ClinMotoInc   0.92642119  0.45911651674238 4.360875e-02
    ## 9           A2 ~~          A2   1.00000000 0.221957242592725 6.625669e-06
    ## 12          A3 ~~          A3   1.00000000 0.298426084011083 8.054545e-04
    ## 6           A1 ~~          A2   0.76674986 0.543575799748169 1.583726e-01
    ## 7           A1 ~~          A3   0.10837347 0.232770979274413 6.415161e-01
    ## 10          A2 ~~          A3  -0.02172102 0.183395936615606 9.057213e-01
    ## 8           A2 =~  ClinAppDec   1.00000000                             NA
    ## 11          A3 =~     ClinSui   1.00000000                             NA
    ## 5           A1 ~~          A1   1.00000000                             NA
    ## 13  ClinAppDec ~~  ClinAppDec   0.00000000                             NA
    ## 19     ClinSui ~~     ClinSui   0.00000000                             NA

## ALSPAC-UKB (Population)

### Common factor

Common factor model

``` r
pop_commonfactor.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
"
pop_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.561 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_commonfactor.fit$modelfit
```

    ##       chisq df     p_chisq      AIC       CFI      SRMR
    ## df 63.07789 35 0.002499044 103.0779 0.9900159 0.1237608

``` r
pop_commonfactor.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep  0.863093526 0.0513537429138102
    ## 1         A1 =~    PopAnh  0.996351120 0.0451433315112963
    ## 2         A1 =~ PopAppDec  0.168551694 0.0875737215245853
    ## 3         A1 =~ PopAppInc  0.403901792 0.0779505460155238
    ## 8         A1 =~ PopSleDec  0.580605068  0.097676101382709
    ## 9         A1 =~ PopSleInc  0.516377956  0.105929022429288
    ## 6         A1 =~  PopFatig  0.664989934 0.0963644775796397
    ## 7         A1 =~  PopGuilt  0.611506224 0.0780285220737901
    ## 4         A1 =~   PopConc  0.707478496 0.0970336138661283
    ## 10        A1 =~    PopSui  0.589830386  0.099771488159136
    ## 16    PopDep ~~    PopDep  0.255069170 0.0792683056178718
    ## 12    PopAnh ~~    PopAnh  0.007284054 0.0659433606867781
    ## 13 PopAppDec ~~ PopAppDec  0.971590824  0.242915351472835
    ## 14 PopAppInc ~~ PopAppInc  0.836863445  0.160437762656559
    ## 19 PopSleDec ~~ PopSleDec  0.662899286  0.314160408130089
    ## 20 PopSleInc ~~ PopSleInc  0.733353922  0.294455468907942
    ## 17  PopFatig ~~  PopFatig  0.557788608  0.337483580540442
    ## 18  PopGuilt ~~  PopGuilt  0.626060735   0.15863027193438
    ## 15   PopConc ~~   PopConc  0.499474088  0.263125445480723
    ## 21    PopSui ~~    PopSui  0.652098387  0.257140069016038
    ## 11        A1 ~~        A1  1.000000000

Remove common variance shared between the gating items (Mood:
`UKB_CIDI1`, Interest: `UKB_CIDI2`) that is uncorrelated with the common
factor variance, to recover the genetic structure among gated items

``` r
pop_commonfactor_gating.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
"
pop_commonfactor_gating.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_gating.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.667 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_gating.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_gating.fit$modelfit
```

    ##      chisq df   p_chisq     AIC       CFI      SRMR
    ## df 40.7478 34 0.1978365 82.7478 0.9976006 0.1166996

``` r
pop_commonfactor_gating.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep    0.7402648 0.0832123815016474
    ## 1         A1 =~    PopAnh    0.8769945 0.0754743520587458
    ## 2         A1 =~ PopAppDec    0.1816943 0.0931950070733182
    ## 3         A1 =~ PopAppInc    0.4348114 0.0855283304440875
    ## 8         A1 =~ PopSleDec    0.6317037  0.107984866769105
    ## 9         A1 =~ PopSleInc    0.5440749  0.114016138469753
    ## 6         A1 =~  PopFatig    0.7183800  0.103746999562885
    ## 7         A1 =~  PopGuilt    0.6597230 0.0852246337165185
    ## 4         A1 =~   PopConc    0.7733652  0.107707393596922
    ## 10        A1 =~    PopSui    0.6237434  0.108177786090839
    ## 16    PopDep ~~    PopAnh    0.2741921  0.117955504536723
    ## 17    PopDep ~~    PopDep    0.4520081  0.135668162784229
    ## 12    PopAnh ~~    PopAnh    0.2308803  0.123780025147359
    ## 13 PopAppDec ~~ PopAppDec    0.9669844  0.242561086380125
    ## 14 PopAppInc ~~ PopAppInc    0.8109383  0.163479374823268
    ## 20 PopSleDec ~~ PopSleDec    0.6009503  0.316643324980521
    ## 21 PopSleInc ~~ PopSleInc    0.7039812  0.295645166009199
    ## 18  PopFatig ~~  PopFatig    0.4839314  0.340028515612521
    ## 19  PopGuilt ~~  PopGuilt    0.5647655  0.162219102690432
    ## 15   PopConc ~~   PopConc    0.4019080  0.264330305642728
    ## 22    PopSui ~~    PopSui    0.6109441  0.261434581322356
    ## 11        A1 ~~        A1    1.0000000

Check if model is improved by allowing residual correlations between the
directional symptoms.

``` r
pop_commonfactor_app.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
"
pop_commonfactor_app.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.619 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_app.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_app.fit$modelfit
```

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 39.46923 33 0.2031463 83.46923 0.9976996 0.1143656

``` r
pop_commonfactor_sle.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
"
pop_commonfactor_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.555 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_sle.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_sle.fit$modelfit
```

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 38.41268 33 0.2375901 82.41268 0.9980753 0.1057011

``` r
pop_commonfactor_app_sle.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
"
pop_commonfactor_app_sle.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_app_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   9.943 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_commonfactor_app_sle.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_commonfactor_app_sle.fit$modelfit
```

    ##       chisq df   p_chisq      AIC      CFI      SRMR
    ## df 36.84266 32 0.2547249 82.84266 0.998278 0.1031013

### Cognitive-Mood-Neuroveg (Kendler Neale) model

``` r
pop_cog_mood_neuroveg.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
"
pop_cog_mood_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg.model)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##          lhs op       rhs Unstandardized_Estimate          SE
    ## 1         A1 =~  PopGuilt             0.134137352 0.089426221
    ## 2         A1 =~   PopConc             0.181602211 0.032239795
    ## 3         A1 =~    PopSui             0.108931581 0.022362524
    ## 4         A2 =~    PopDep             0.252969647 0.014461091
    ## 5         A2 =~    PopAnh             0.309200175 0.014241933
    ## 6         A2 =~  PopGuilt             0.019720341 0.090679998
    ## 7         A3 =~ PopSleDec             0.107608917 0.022370423
    ## 8         A3 =~ PopSleInc             0.078739005 0.019146116
    ## 9         A3 =~  PopFatig             0.157627706 0.030546436
    ## 10        A3 =~ PopAppDec             0.030563065 0.015888144
    ## 11        A3 =~ PopAppInc             0.093801722 0.022703482
    ## 15  PopGuilt ~~  PopGuilt             0.037080005 0.010415243
    ## 16   PopConc ~~   PopConc             0.025481389 0.016826585
    ## 17    PopSui ~~    PopSui             0.021382225 0.008720044
    ## 18    PopDep ~~    PopDep             0.018240816 0.006203009
    ## 19    PopAnh ~~    PopAnh            -0.008236601 0.006941109
    ## 20 PopSleDec ~~ PopSleDec             0.022514470 0.010881372
    ## 21 PopSleInc ~~ PopSleInc             0.018614123 0.007573042
    ## 22  PopFatig ~~  PopFatig             0.032875700 0.020288414
    ## 23 PopAppDec ~~ PopAppDec             0.032784142 0.008188596
    ## 24 PopAppInc ~~ PopAppInc             0.045648645 0.009203128
    ## 25        A1 ~~        A2             0.862046515 0.146615552
    ## 26        A1 ~~        A3             1.234810379 0.248062517
    ## 27        A2 ~~        A3             0.918303623 0.157290398

Add constraints to prevent variances from being negative and
correlations from going out of bounds.

``` r
pop_cog_mood_neuroveg_constr.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
c2 > 0.001
PopAnh ~~ c2*PopAnh
a13 < 0.99
A1 ~~ a13*A3
"
pop_cog_mood_neuroveg_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   20.46 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_cog_mood_neuroveg_constr.model, : A difference greater than .025 was
    ## observed pre- and post-smoothing for Z-statistics in the genetic covariance
    ## matrix. This reflects a large difference and results should be interpreted
    ## with caution!! This can often result from including low powered traits, and you
    ## might consider removing those traits from the model. If you are going to run
    ## a multivariate GWAS we strongly recommend setting the smooth_check argument to
    ## true to check smoothing for each SNP.

``` r
pop_cog_mood_neuroveg_constr.fit$modelfit
```

    ##      chisq df    p_chisq     AIC       CFI      SRMR
    ## df 48.6032 31 0.02301468 96.6032 0.9937406 0.1172712

``` r
pop_cog_mood_neuroveg_constr.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 2         A1 =~  PopGuilt  0.769761002  0.625782101381628
    ## 1         A1 =~   PopConc  0.773881046  0.132766120087573
    ## 3         A1 =~    PopSui  0.629339403  0.124924935934997
    ## 8         A2 =~    PopDep  0.893950083 0.0508618530305653
    ## 7         A2 =~    PopAnh  1.010127746 0.0471688929975766
    ## 9         A2 =~  PopGuilt -0.097070891   0.62359152234629
    ## 15        A3 =~ PopSleDec  0.623336377  0.121020462176279
    ## 16        A3 =~ PopSleInc  0.536774469  0.124184657130501
    ## 14        A3 =~  PopFatig  0.707324529  0.125069397287854
    ## 12        A3 =~ PopAppDec  0.179502580  0.091784801047111
    ## 13        A3 =~ PopAppInc  0.427656396  0.097618335803503
    ## 18    PopAnh ~~    PopAnh  0.001000007 0.0733397835706942
    ## 6         A1 ~~        A3  0.989999989  0.171448122455602
    ## 24  PopGuilt ~~  PopGuilt  0.524625888  0.242235902916805
    ## 21   PopConc ~~   PopConc  0.401107466  0.282143186272646
    ## 27    PopSui ~~    PopSui  0.603930127  0.268746857643935
    ## 22    PopDep ~~    PopDep  0.200853001 0.0757173072891572
    ## 25 PopSleDec ~~ PopSleDec  0.611450384  0.320679500208172
    ## 26 PopSleInc ~~ PopSleInc  0.711871821   0.30697148711204
    ## 23  PopFatig ~~  PopFatig  0.499690639  0.354142841897334
    ## 19 PopAppDec ~~ PopAppDec  0.967778676  0.242807524248992
    ## 20 PopAppInc ~~ PopAppInc  0.817109909  0.170658670818787
    ## 5         A1 ~~        A2  0.847024769  0.140129804655695
    ## 11        A2 ~~        A3  0.869983395  0.135188493668608
    ## 4         A1 ~~        A1  1.000000000                   
    ## 10        A2 ~~        A2  1.000000000                   
    ## 17        A3 ~~        A3  1.000000000

### Two-factor models

[Elhai Psychiat Res
2012](https://www.sciencedirect.com/science/article/pii/S0165178112002685)
compared 3 two-factor models

### Psychological-Somatic (Elhai Model 2a)

[Krause Rehab Psychol
2008](https://psycnet.apa.org/record/2008-17022-011), [Krause Arch Psys
Med Rehab
2010](https://www.sciencedirect.com/science/article/pii/S0003999310002443):

> the 2-factor solution with 3 somatic items (sleep disturbance, poor
> energy, appetite change) was a better solution than either a
> unidimensional model or 2-factor model that included psychomotor
> slowing as a fourth somatic item

``` r
pop_psych_soma.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui 
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
"
pop_psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##          lhs op       rhs Unstandardized_Estimate          SE
    ## 1         A1 =~    PopDep              0.20764171 0.023107667
    ## 2         A1 =~    PopAnh              0.25423902 0.021637051
    ## 3         A1 =~  PopGuilt              0.16196617 0.020788343
    ## 4         A1 =~   PopConc              0.18859556 0.025914409
    ## 5         A1 =~    PopSui              0.11370079 0.019688931
    ## 6         A2 =~ PopAppDec              0.03059595 0.015890477
    ## 7         A2 =~ PopAppInc              0.09361839 0.022671944
    ## 8         A2 =~ PopSleDec              0.10739836 0.022342452
    ## 9         A2 =~ PopSleInc              0.07897112 0.019174033
    ## 10        A2 =~  PopFatig              0.15783688 0.030562436
    ## 13    PopDep ~~    PopAnh              0.02547907 0.009594246
    ## 14    PopDep ~~    PopDep              0.03911940 0.010820913
    ## 15    PopAnh ~~    PopAnh              0.02273067 0.010376362
    ## 16  PopGuilt ~~  PopGuilt              0.03378939 0.009721815
    ## 17   PopConc ~~   PopConc              0.02289253 0.015355123
    ## 18    PopSui ~~    PopSui              0.02032043 0.008683486
    ## 19 PopAppDec ~~ PopAppDec              0.03278214 0.008190261
    ## 20 PopAppInc ~~ PopAppInc              0.04568307 0.009200785
    ## 21 PopSleDec ~~ PopSleDec              0.02255973 0.010882202
    ## 22 PopSleInc ~~ PopSleInc              0.01857750 0.007579743
    ## 23  PopFatig ~~  PopFatig              0.03280988 0.020291600
    ## 24        A1 ~~        A2              1.13026676 0.170834984

``` r
pop_psych_soma_constr.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui 
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
a12 < 1
A1 ~~ a12*A2
"
pop_psych_soma_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma_constr.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  11.208 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_soma_constr.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_psych_soma_constr.fit$modelfit
```

    ##      chisq df   p_chisq     AIC      CFI      SRMR
    ## df 40.7478 33 0.1663865 84.7478 0.997245 0.1166996

``` r
pop_psych_soma_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 3         A1 =~    PopDep    0.7402658 0.0829133067845549
    ## 1         A1 =~    PopAnh    0.8769949 0.0754767590892474
    ## 4         A1 =~  PopGuilt    0.6597223 0.0851875828836091
    ## 2         A1 =~   PopConc    0.7733648  0.107640262975317
    ## 5         A1 =~    PopSui    0.6237430   0.10820101245641
    ## 8         A2 =~ PopAppDec    0.1816942 0.0928542859069906
    ## 9         A2 =~ PopAppInc    0.4348107 0.0979114575875578
    ## 11        A2 =~ PopSleDec    0.6317042  0.121145690669268
    ## 12        A2 =~ PopSleInc    0.5440753  0.124593402942297
    ## 10        A2 =~  PopFatig    0.7183807  0.124508645624776
    ## 18    PopDep ~~    PopAnh    0.2741910  0.118477735352797
    ## 7         A1 ~~        A2    0.9999995  0.133347935537963
    ## 19    PopDep ~~    PopDep    0.4520067  0.136589690351285
    ## 14    PopAnh ~~    PopAnh    0.2308797     0.124506781527
    ## 21  PopGuilt ~~  PopGuilt    0.5647662  0.161926536048416
    ## 17   PopConc ~~   PopConc    0.4019091  0.263084646226254
    ## 24    PopSui ~~    PopSui    0.6109451  0.261229339752704
    ## 15 PopAppDec ~~ PopAppDec    0.9669869  0.242738715219611
    ## 16 PopAppInc ~~ PopAppInc    0.8109385  0.171091199045209
    ## 22 PopSleDec ~~ PopSleDec    0.6009408  0.320757331559754
    ## 23 PopSleInc ~~ PopSleInc    0.7039765  0.307187837495002
    ## 20  PopFatig ~~  PopFatig    0.4839152  0.354353093518351
    ## 6         A1 ~~        A1    1.0000000                   
    ## 13        A2 ~~        A2    1.0000000

Bifactor model

``` r
pop_psych_soma_bif.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A  =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A  ~~ 0*A1
A  ~~ 0*A2
A1 ~~ 0*A2
PopDep ~~ PopAnh
"
pop_psych_soma_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.916 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_soma_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_psych_soma_bif.fit$model
```

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 21.38632 24 0.6158438 83.38632   1 0.0655899

``` r
pop_psych_soma_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 14        A1 =~    PopDep   0.40458812  0.189102175443417
    ## 12        A1 =~    PopAnh   0.15021980  0.155710848556851
    ## 15        A1 =~  PopGuilt   0.29655644  0.181017955005552
    ## 13        A1 =~   PopConc  -0.29818200  0.261424342535502
    ## 16        A1 =~    PopSui   0.74594753  0.347298583758189
    ## 20        A2 =~ PopAppDec   0.11795426  0.202197301070381
    ## 21        A2 =~ PopAppInc   0.17708609   0.16966839831371
    ## 23        A2 =~ PopSleDec   0.70473550  0.409346819925278
    ## 24        A2 =~ PopSleInc  -0.61383273   0.30975673400711
    ## 22        A2 =~  PopFatig  -0.37436989  0.250197845328412
    ## 5          A =~    PopDep   0.65557088 0.0922828488031419
    ## 1          A =~    PopAnh   0.82655874 0.0743470363602854
    ## 7          A =~  PopGuilt   0.60812411 0.0978404196625757
    ## 4          A =~   PopConc   0.92169208  0.125993245037559
    ## 10         A =~    PopSui   0.48553441  0.144288683723201
    ## 2          A =~ PopAppDec   0.18762814 0.0964816247411554
    ## 3          A =~ PopAppInc   0.45855117 0.0901358703872597
    ## 8          A =~ PopSleDec   0.70165833  0.119571986912166
    ## 9          A =~ PopSleInc   0.59030485  0.118916773045259
    ## 6          A =~  PopFatig   0.76663623  0.109345604713364
    ## 31    PopDep ~~    PopAnh   0.32075513  0.120624063387125
    ## 32    PopDep ~~    PopDep   0.40653567  0.181070769658493
    ## 27    PopAnh ~~    PopAnh   0.29423435  0.113034165662565
    ## 34  PopGuilt ~~  PopGuilt   0.54223928  0.173087874467872
    ## 30   PopConc ~~   PopConc   0.06157141    0.3371896498397
    ## 37    PopSui ~~    PopSui   0.20781834  0.528498378580611
    ## 28 PopAppDec ~~ PopAppDec   0.95088163  0.240833945503497
    ## 29 PopAppInc ~~ PopAppInc   0.75837120  0.186748762142194
    ## 35 PopSleDec ~~ PopSleDec   0.01102419  0.695993962945813
    ## 36 PopSleInc ~~ PopSleInc   0.27474890  0.478763119505728
    ## 33  PopFatig ~~  PopFatig   0.27211792  0.363309995042641
    ## 18        A1 ~~        A1   1.00000000                   
    ## 26        A2 ~~        A2   1.00000000                   
    ## 11         A ~~         A   1.00000000                   
    ## 17        A1 ~~         A   0.00000000                   
    ## 25        A2 ~~         A   0.00000000                   
    ## 19        A1 ~~        A2   0.00000000

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
pop_psych_veg.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
"
pop_psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.782 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
pop_psych_veg.fit$modelfit
```

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 40.91469 33 0.1619784 84.91469 0.9971857 0.1136106

``` r
pop_psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 2         A1 =~    PopDep    0.7751677 0.0917873496754169
    ## 1         A1 =~    PopAnh    0.9135854 0.0887692244153966
    ## 3         A1 =~  PopGuilt    0.6528597 0.0856905972816179
    ## 4         A1 =~    PopSui    0.6217293  0.108225022862179
    ## 7         A2 =~ PopAppDec    0.1890142 0.0972511105406584
    ## 8         A2 =~ PopAppInc    0.4576341 0.0930447857780165
    ## 11        A2 =~ PopSleDec    0.6646065  0.115795099333534
    ## 12        A2 =~ PopSleInc    0.5678768  0.122491604131194
    ## 10        A2 =~  PopFatig    0.7566586  0.114516499980711
    ## 9         A2 =~   PopConc    0.8219150  0.115476668155541
    ## 18    PopDep ~~    PopAnh    0.2152187  0.140046836045944
    ## 19    PopDep ~~    PopDep    0.3991155  0.152718249225757
    ## 14    PopAnh ~~    PopAnh    0.1653618  0.152066565141075
    ## 21  PopGuilt ~~  PopGuilt    0.5737749  0.161629174804251
    ## 24    PopSui ~~    PopSui    0.6134539  0.260608791525752
    ## 15 PopAppDec ~~ PopAppDec    0.9642734  0.242529162988306
    ## 16 PopAppInc ~~ PopAppInc    0.7905711  0.168178076991092
    ## 22 PopSleDec ~~ PopSleDec    0.5582988  0.321376248115576
    ## 23 PopSleInc ~~ PopSleInc    0.6775179  0.301244483291873
    ## 20  PopFatig ~~  PopFatig    0.4274689  0.342516786451423
    ## 17   PopConc ~~   PopConc    0.3244553  0.270203969433831
    ## 6         A1 ~~        A2    0.8917916 0.0996620470326232
    ## 5         A1 ~~        A1    1.0000000                   
    ## 13        A2 ~~        A2    1.0000000

``` r
pop_psych_veg_bif.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
"
pop_psych_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.066 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_psych_veg_bif.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI       SRMR
    ## df 15.50072 23 0.8757949 79.50072   1 0.07257442

``` r
pop_psych_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 13        A1 =~    PopDep    0.4753894  0.205240347277867
    ## 12        A1 =~    PopAnh    0.1981545  0.164898102546413
    ## 14        A1 =~  PopGuilt    0.3641352  0.170861245322782
    ## 15        A1 =~    PopSui    0.7583570  0.343933998475301
    ## 19        A2 =~ PopAppDec   -0.1527297  0.192653734593263
    ## 20        A2 =~ PopAppInc    0.5533866  0.313254600965919
    ## 23        A2 =~ PopSleDec    0.2301991  0.243942540946892
    ## 24        A2 =~ PopSleInc   -0.3071509  0.262611885974934
    ## 22        A2 =~  PopFatig   -0.3836750  0.304147338736842
    ## 21        A2 =~   PopConc    0.3122732  0.241961619178714
    ## 5          A =~    PopDep    0.6112869 0.0925636309325131
    ## 1          A =~    PopAnh    0.8181285 0.0846501478513235
    ## 7          A =~  PopGuilt    0.5745891   0.09969235255766
    ## 10         A =~    PopSui    0.4271878  0.139825944531391
    ## 2          A =~ PopAppDec    0.1983106  0.100301824563354
    ## 3          A =~ PopAppInc    0.4723537 0.0996314991769476
    ## 8          A =~ PopSleDec    0.6983005  0.122136151800446
    ## 9          A =~ PopSleInc    0.6401273  0.128484577081147
    ## 6          A =~  PopFatig    0.8145768  0.119423727013029
    ## 4          A =~   PopConc    0.8401019  0.118856740349479
    ## 31    PopDep ~~    PopAnh    0.3290888  0.118125682708932
    ## 36 PopSleDec ~~ PopSleInc   -0.3815444  0.255965714139232
    ## 32    PopDep ~~    PopDep    0.4003335  0.194858503236901
    ## 27    PopAnh ~~    PopAnh    0.2914006  0.115790772772666
    ## 34  PopGuilt ~~  PopGuilt    0.5372526  0.177721763680975
    ## 38    PopSui ~~    PopSui    0.2424054  0.512184978166035
    ## 28 PopAppDec ~~ PopAppDec    0.9373450  0.249524993052742
    ## 29 PopAppInc ~~ PopAppInc    0.4706446  0.367302366759667
    ## 35 PopSleDec ~~ PopSleDec    0.4593841  0.342290294902112
    ## 37 PopSleInc ~~ PopSleInc    0.4958952  0.362130767831745
    ## 33  PopFatig ~~  PopFatig    0.1892588  0.431633576167828
    ## 30   PopConc ~~   PopConc    0.1967126  0.314290743991356
    ## 17        A1 ~~        A1    1.0000000                   
    ## 26        A2 ~~        A2    1.0000000                   
    ## 11         A ~~         A    1.0000000                   
    ## 16        A1 ~~         A    0.0000000                   
    ## 25        A2 ~~         A    0.0000000                   
    ## 18        A1 ~~        A2    0.0000000

### Affective-Neurovegetative (Elhai Model 2c)

``` r
pop_affect_veg.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
"
pop_affect_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.953 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_affect_veg.fit$modelfit
```

    ##      chisq df   p_chisq     AIC CFI      SRMR
    ## df 31.6321 33 0.5351805 75.6321   1 0.1104974

``` r
pop_affect_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7951436 0.0874243873773166
    ## 2         A1 =~  PopGuilt    0.7215656 0.0943260550440023
    ## 3         A1 =~    PopSui    0.6895697  0.119445984241806
    ## 6         A2 =~    PopAnh    0.8847304  0.075704430384117
    ## 8         A2 =~ PopAppInc    0.4527022 0.0883425846843872
    ## 7         A2 =~ PopAppDec    0.1857302 0.0960461566038244
    ## 12        A2 =~ PopSleInc    0.5627518  0.117128000227501
    ## 11        A2 =~ PopSleDec    0.6505323  0.111558728817932
    ## 10        A2 =~  PopFatig    0.7464879   0.10641044494472
    ## 9         A2 =~   PopConc    0.8068758  0.109697986879609
    ## 18    PopDep ~~    PopAnh    0.3484041  0.110936868275062
    ## 19    PopDep ~~    PopDep    0.3677469    0.1494538802891
    ## 21  PopGuilt ~~  PopGuilt    0.4793437  0.177586936002648
    ## 24    PopSui ~~    PopSui    0.5244939  0.272510108423568
    ## 14    PopAnh ~~    PopAnh    0.2172521  0.125548062310234
    ## 16 PopAppInc ~~ PopAppInc    0.7950607   0.16507292638154
    ## 15 PopAppDec ~~ PopAppDec    0.9655045  0.242420488392951
    ## 23 PopSleInc ~~ PopSleInc    0.6833089  0.297427083975384
    ## 22 PopSleDec ~~ PopSleDec    0.5768077  0.317224830556687
    ## 20  PopFatig ~~  PopFatig    0.4427557  0.338290780840395
    ## 17   PopConc ~~   PopConc    0.3489517  0.265110220553829
    ## 5         A1 ~~        A2    0.8173500  0.061561835023154
    ## 4         A1 ~~        A1    1.0000000                   
    ## 13        A2 ~~        A2    1.0000000

Bifactor model

``` r
pop_affect_veg_bif.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A =~ NA*PopDep + PopGuilt + PopSui + PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
PopDep ~~ PopAnh
"
pop_affect_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.864 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pop_affect_veg_bif_constr.model <- "
a2 < 1
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A =~ NA*PopDep + PopGuilt + PopSui + a2*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
PopDep ~~ PopAnh
c4a > 0.001
PopSleDec ~~ c4a*PopSleDec
c4b > 0.001
PopSleInc ~~ c4b*PopSleInc
"
pop_affect_veg_bif_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg_bif_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  17.777 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0019152623787771 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.179489390573754 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif_constr.model, : A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_affect_veg_bif_constr.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI       SRMR
    ## df 17.49647 24 0.8267298 79.49647   1 0.06727487

``` r
pop_affect_veg_bif_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 12        A1 =~    PopDep  0.443619232  0.219857552706333
    ## 13        A1 =~  PopGuilt  0.314149272  0.159461392036603
    ## 14        A1 =~    PopSui  0.667394247  0.324378054550841
    ## 18        A2 =~    PopAnh  0.192786106  0.110445886442298
    ## 20        A2 =~ PopAppInc -0.139102609  0.164993709813016
    ## 19        A2 =~ PopAppDec -0.114468042  0.190024740126896
    ## 24        A2 =~ PopSleInc  0.600862849  0.253115722041161
    ## 23        A2 =~ PopSleDec -0.710346417  0.362934848614404
    ## 22        A2 =~  PopFatig  0.361522455  0.235747763208688
    ## 21        A2 =~   PopConc -0.152272032  0.200916031051972
    ## 5          A =~    PopDep  0.633963226  0.085925852064305
    ## 7          A =~  PopGuilt  0.610798008 0.0899310467895853
    ## 10         A =~    PopSui  0.510708593  0.111216341506241
    ## 1          A =~    PopAnh  0.874900524  0.083873829498763
    ## 3          A =~ PopAppInc  0.465373037 0.0908424959675696
    ## 2          A =~ PopAppDec  0.195414554 0.0973131495699587
    ## 9          A =~ PopSleInc  0.546013814  0.130226372475933
    ## 8          A =~ PopSleDec  0.743176249  0.134937565172349
    ## 6          A =~  PopFatig  0.732223906  0.111853303592385
    ## 4          A =~   PopConc  0.828196089   0.11336880206518
    ## 31    PopDep ~~    PopAnh  0.368745751  0.113564916291489
    ## 35 PopSleDec ~~ PopSleDec  0.000999724  0.650087615461865
    ## 36 PopSleInc ~~ PopSleInc  0.340830481  0.394108166756054
    ## 32    PopDep ~~    PopDep  0.401292662  0.204782758069655
    ## 34  PopGuilt ~~  PopGuilt  0.528236029  0.178156440829026
    ## 37    PopSui ~~    PopSui  0.293762284  0.451186524671523
    ## 27    PopAnh ~~    PopAnh  0.197382834  0.142560163304573
    ## 29 PopAppInc ~~ PopAppInc  0.764079413  0.182938059346263
    ## 28 PopAppDec ~~ PopAppDec  0.948708858  0.240071593852727
    ## 33  PopFatig ~~  PopFatig  0.333146875   0.34932915465057
    ## 30   PopConc ~~   PopConc  0.290903885  0.273873395531704
    ## 16        A1 ~~        A1  1.000000000                   
    ## 26        A2 ~~        A2  1.000000000                   
    ## 11         A ~~         A  1.000000000                   
    ## 15        A1 ~~         A  0.000000000                   
    ## 25        A2 ~~         A  0.000000000                   
    ## 17        A1 ~~        A2  0.000000000

### Model comparisons

``` r
pop_model_list=
list("1a"=list(name="Common", model=pop_commonfactor.fit),
     "1b"=list(name="Common (gating)", model=pop_commonfactor_gating.fit),
     "1c"=list(name="Common (App)", model=pop_commonfactor_app.fit),
     "1d"=list(name="Common (Sle)", model=pop_commonfactor_sle.fit),
     "1er"=list(name="Common (App,Sle)", model=pop_commonfactor_app_sle.fit),
     "2a(i)"=list(name="Psych-Somatic", model=pop_psych_soma_constr.fit),
     "2a(ii)"=list(name="Psych-Somatic (BiF)", model=pop_psych_soma_bif.fit),
     "2b(i)"=list(name="Psych-Neuroveg", model=pop_psych_veg.fit),
     "2b(ii)"=list(name="Psych-Neuroveg (BiF)", model=pop_psych_veg_bif.fit),
     "2c(i)"=list(name="Affect-Neuroveg", model=pop_affect_veg.fit),
     "2c(ii)"=list(name="Affect-Neuroveg (BiF)", model=pop_affect_veg_bif_constr.fit),
     "3"=list(name="Cog-Mood-Neuroveg", model=pop_cog_mood_neuroveg_constr.fit)
     )

pop_model_fits <- 
data.frame(Model=names(pop_model_list),
           Name=sapply(pop_model_list, function(m) m$name)) %>%
bind_cols(
bind_rows(
lapply(pop_model_list, function(m) m$model$modelfit)
))
rownames(pop_model_fits) <- NULL

knitr::kable(
pop_model_fits %>%
mutate(dAIC=AIC-min(AIC)) %>%
mutate_if(is.numeric, ~signif(., 3))
)
```

| Model  | Name                  | chisq |  df | p_chisq |   AIC |   CFI |   SRMR |  dAIC |
|:-------|:----------------------|------:|----:|--------:|------:|------:|-------:|------:|
| 1a     | Common                |  63.1 |  35 |  0.0025 | 103.0 | 0.990 | 0.1240 | 27.40 |
| 1b     | Common (gating)       |  40.7 |  34 |  0.1980 |  82.7 | 0.998 | 0.1170 |  7.12 |
| 1c     | Common (App)          |  39.5 |  33 |  0.2030 |  83.5 | 0.998 | 0.1140 |  7.84 |
| 1d     | Common (Sle)          |  38.4 |  33 |  0.2380 |  82.4 | 0.998 | 0.1060 |  6.78 |
| 1er    | Common (App,Sle)      |  36.8 |  32 |  0.2550 |  82.8 | 0.998 | 0.1030 |  7.21 |
| 2a(i)  | Psych-Somatic         |  40.7 |  33 |  0.1660 |  84.7 | 0.997 | 0.1170 |  9.12 |
| 2a(ii) | Psych-Somatic (BiF)   |  21.4 |  24 |  0.6160 |  83.4 | 1.000 | 0.0656 |  7.75 |
| 2b(i)  | Psych-Neuroveg        |  40.9 |  33 |  0.1620 |  84.9 | 0.997 | 0.1140 |  9.28 |
| 2b(ii) | Psych-Neuroveg (BiF)  |  15.5 |  23 |  0.8760 |  79.5 | 1.000 | 0.0726 |  3.87 |
| 2c(i)  | Affect-Neuroveg       |  31.6 |  33 |  0.5350 |  75.6 | 1.000 | 0.1100 |  0.00 |
| 2c(ii) | Affect-Neuroveg (BiF) |  17.5 |  24 |  0.8270 |  79.5 | 1.000 | 0.0673 |  3.86 |
| 3      | Cog-Mood-Neuroveg     |  48.6 |  31 |  0.0230 |  96.6 | 0.994 | 0.1170 | 21.00 |

# Clinical and population factors

``` r
clin_pop_affect_veg.model <- "
ClinSOMA =~ NA*ClinAppInc + ClinSleDec + ClinSleInc + ClinMotoInc
ClinAPPDEC =~ ClinAppDec
ClinSUI =~ ClinSui
PopAFFECT =~ NA*PopDep + PopGuilt + PopSui
PopNEUROVEG =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
ClinSOMA ~~ 1*ClinSOMA
PopAFFECT ~~ 1*PopAFFECT
PopNEUROVEG ~~ 1*PopNEUROVEG
PopDep ~~ PopAnh
ClinAppDec ~~ ClinAppInc
ClinAppDec ~~ 0*ClinAppDec
ClinSui ~~ 0*ClinSui
"

clin_pop_affect_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_pop_affect_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   3.607 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0313666288930831 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.49543224219422 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_pop_affect_veg.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_pop_affect_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
clin_pop_affect_veg.fit$modelfit
```

    ##      chisq df      p_chisq     AIC       CFI      SRMR
    ## df 176.887 94 5.062223e-07 260.887 0.9814401 0.1535288

``` r
clin_pop_affect_veg.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 12    ClinSOMA =~  ClinAppInc   0.65009389  0.256391862943174 1.122669e-02
    ## 14    ClinSOMA =~  ClinSleDec   0.47030640  0.229438293957574 4.038137e-02
    ## 15    ClinSOMA =~  ClinSleInc   0.45808267   0.28236580566831 1.047363e-01
    ## 13    ClinSOMA =~ ClinMotoInc   0.48655965  0.211191393351506 2.122864e-02
    ## 26   PopAFFECT =~      PopDep   0.72854299 0.0744481791514393 1.294059e-22
    ## 27   PopAFFECT =~    PopGuilt   0.77133437 0.0934025106914945 1.479084e-16
    ## 28   PopAFFECT =~      PopSui   0.68298415  0.111130768017908 7.957539e-10
    ## 39 PopNEUROVEG =~      PopAnh   0.81386024 0.0651826662237478 8.920100e-36
    ## 41 PopNEUROVEG =~   PopAppInc   0.45779861 0.0823074754341747 2.666196e-08
    ## 40 PopNEUROVEG =~   PopAppDec   0.18256193 0.0894094238401614 4.116469e-02
    ## 45 PopNEUROVEG =~   PopSleInc   0.50026401  0.107834742911666 3.498108e-06
    ## 44 PopNEUROVEG =~   PopSleDec   0.65740043  0.106028465078548 5.638342e-10
    ## 43 PopNEUROVEG =~    PopFatig   0.75145946 0.0995915313035136 4.508715e-14
    ## 42 PopNEUROVEG =~     PopConc   0.79775303   0.10110111188002 3.006264e-15
    ## 35      PopDep ~~      PopAnh   0.42945395 0.0867687804647645 7.444194e-07
    ## 7   ClinAppInc ~~  ClinAppDec  -0.81797467  0.357592223949102 2.216901e-02
    ## 8   ClinAppInc ~~  ClinAppInc   0.57737749  0.448240144590205 1.977165e-01
    ## 10  ClinSleDec ~~  ClinSleDec   0.77881358  0.516801796050139 1.318166e-01
    ## 11  ClinSleInc ~~  ClinSleInc   0.79015827  0.870926242428377 3.642731e-01
    ## 9  ClinMotoInc ~~ ClinMotoInc   0.76325836  0.467601352451282 1.026221e-01
    ## 36      PopDep ~~      PopDep   0.46922548  0.121027015834608 1.057570e-04
    ## 38    PopGuilt ~~    PopGuilt   0.40504301  0.184090520004726 2.779054e-02
    ## 49      PopSui ~~      PopSui   0.53353296  0.249924359911278 3.277944e-02
    ## 31      PopAnh ~~      PopAnh   0.33763227  0.100356259744396 7.673314e-04
    ## 33   PopAppInc ~~   PopAppInc   0.79042175  0.150168577925681 1.412990e-07
    ## 32   PopAppDec ~~   PopAppDec   0.96667203  0.222495775859501 1.394851e-05
    ## 48   PopSleInc ~~   PopSleInc   0.74973560  0.255109248451215 3.293963e-03
    ## 47   PopSleDec ~~   PopSleDec   0.56782447  0.292459197230163 5.219107e-02
    ## 37    PopFatig ~~    PopFatig   0.43530859  0.311259513192369 1.619494e-01
    ## 34     PopConc ~~     PopConc   0.36359065  0.247775221476709 1.422634e-01
    ## 3   ClinAPPDEC ~~  ClinAPPDEC   1.00000000   0.22444905317841 8.375486e-06
    ## 23     ClinSUI ~~     ClinSUI   1.00000000   0.28778357857873 5.111707e-04
    ## 16    ClinSOMA ~~  ClinAPPDEC   0.64119216  0.397836969703455 1.070280e-01
    ## 18    ClinSOMA ~~     ClinSUI   0.10313319  0.253792035648623 6.844736e-01
    ## 19    ClinSOMA ~~   PopAFFECT   0.09299871    0.1445324431885 5.199369e-01
    ## 20    ClinSOMA ~~ PopNEUROVEG   0.45264558  0.183230600477497 1.349744e-02
    ## 4   ClinAPPDEC ~~     ClinSUI  -0.02413193  0.181104212348419 8.939697e-01
    ## 5   ClinAPPDEC ~~   PopAFFECT  -0.10829522  0.104739940233647 3.011671e-01
    ## 6   ClinAPPDEC ~~ PopNEUROVEG  -0.10044531  0.108469850131147 3.544365e-01
    ## 24     ClinSUI ~~   PopAFFECT   0.69330977  0.121935397074887 1.301355e-08
    ## 25     ClinSUI ~~ PopNEUROVEG   0.57253622  0.113837601127897 4.919935e-07
    ## 30   PopAFFECT ~~ PopNEUROVEG   0.80479374 0.0626473462877868 9.009547e-38
    ## 1   ClinAPPDEC =~  ClinAppDec   1.00000000                              NA
    ## 21     ClinSUI =~     ClinSui   1.00000000                              NA
    ## 17    ClinSOMA ~~    ClinSOMA   1.00000000                              NA
    ## 29   PopAFFECT ~~   PopAFFECT   1.00000000                              NA
    ## 46 PopNEUROVEG ~~ PopNEUROVEG   1.00000000                              NA
    ## 2   ClinAppDec ~~  ClinAppDec   0.00000000                              NA
    ## 22     ClinSui ~~     ClinSui   0.00000000                              NA

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
symptoms_cov_pd <- as.matrix(Matrix::nearPD(symptoms_cov_pos, corr=FALSE)$mat)

corrplot(cov2cor(symptoms_cov_pd))
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_pd-1.png)<!-- -->

## PGC/AGDS

Check eigen values of the correlation matrix

``` r
symptoms_clin_idx <- which(str_detect(dimnames(symptoms_cov_pd)[[1]], 'Clin'))

symptoms_clin_eigen <- eigen(cov2cor(symptoms_cov_pd[symptoms_clin_idx,symptoms_clin_idx])) 

plot(symptoms_clin_eigen$values, ylab='Eigenvalue')
lines(symptoms_clin_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_clin_efa_eigen-1.png)<!-- -->

``` r
symptoms_clin_efa <- factanal(covmat=symptoms_cov_pd[symptoms_clin_idx,symptoms_clin_idx], factors=3, rotation='varimax')
symptoms_clin_efa
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd[symptoms_clin_idx,     symptoms_clin_idx], rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc     ClinSui 
    ##       0.005       0.005       0.005       0.771       0.648       0.714 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec   0.271   0.951   0.132 
    ## ClinAppInc   0.651  -0.681   0.326 
    ## ClinSleDec   0.982   0.161         
    ## ClinSleInc   0.478                 
    ## ClinMotoInc  0.148   0.103   0.565 
    ## ClinSui     -0.130           0.513 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      1.728   1.412   0.712
    ## Proportion Var   0.288   0.235   0.119
    ## Cumulative Var   0.288   0.523   0.642
    ## 
    ## The degrees of freedom for the model is 0 and the fit was 0.3701

## ALSPAC/UKB

Check eigen values of the correlation matrix

``` r
symptoms_pop_idx <- which(str_detect(dimnames(symptoms_cov_pd)[[1]], 'Pop'))

symptoms_pop_eigen <- eigen(cov2cor(symptoms_cov_pd[symptoms_pop_idx,symptoms_pop_idx])) 

plot(symptoms_pop_eigen$values, ylab='Eigenvalue')
lines(symptoms_pop_eigen$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_pop_efa_eigen-1.png)<!-- -->

``` r
symptoms_pop_efa <- factanal(covmat=symptoms_cov_pd[symptoms_pop_idx,symptoms_pop_idx], factors=3, rotation='varimax')
symptoms_pop_efa
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd[symptoms_pop_idx,     symptoms_pop_idx], rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##     PopDep     PopAnh  PopAppDec  PopAppInc  PopSleDec  PopSleInc PopMotoInc 
    ##      0.005      0.005      0.557      0.702      0.707      0.432      0.789 
    ## PopMotoDec   PopFatig   PopGuilt    PopConc     PopSui 
    ##      0.005      0.520      0.693      0.554      0.584 
    ## 
    ## Loadings:
    ##            Factor1 Factor2 Factor3
    ## PopDep      0.945   0.315         
    ## PopAnh      0.750   0.655         
    ## PopAppDec   0.102   0.136   0.643 
    ## PopAppInc   0.155   0.378  -0.362 
    ## PopSleDec   0.386   0.241  -0.294 
    ## PopSleInc           0.741   0.118 
    ## PopMotoInc                  0.457 
    ## PopMotoDec  0.328  -0.104   0.936 
    ## PopFatig    0.213   0.655         
    ## PopGuilt    0.437   0.322   0.113 
    ## PopConc     0.303   0.593         
    ## PopSui      0.627           0.144 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      2.475   2.195   1.776
    ## Proportion Var   0.206   0.183   0.148
    ## Cumulative Var   0.206   0.389   0.537
    ## 
    ## The degrees of freedom for the model is 33 and the fit was 8.8788

# All symptoms

``` r
all_covstruct_prefix <- 'all.covstruct'
all_covstruct_r <- file.path('ldsc', paste(all_covstruct_prefix, 'deparse.R', sep='.'))
all_covstruct_rds <- file.path('ldsc', paste(all_covstruct_prefix, 'rds', sep='.'))

all_symptoms_covstruct <- dget(all_covstruct_r)

all_sumstats_prevs <- read_tsv(file.path('ldsc', paste(all_covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 12 Columns: 5
    ## ── Column specification ────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (4): filename, sumstats, trait_name, symptom
    ## dbl (1): pop_prev
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
all_cohorts_sample_symptoms <-
all_sumstats_prevs %>%
left_join(dsm_mdd_symptoms_labels, by=c('symptom'='ref')) %>%
select(symptom, trait_name, abbv)

all_sample_symptoms <- all_cohorts_sample_symptoms$abbv
names(all_sample_symptoms) <- all_cohorts_sample_symptoms$trait_name

# rename traits in covstruct
dimnames(all_symptoms_covstruct$S)[[2]] <-
as.vector(all_sample_symptoms[dimnames(all_symptoms_covstruct$S)[[2]]])

# positive variances
dimnames(all_symptoms_covstruct$S)[[2]][which(diag(all_symptoms_covstruct$S) > 0)]
```

    ##  [1] "Dep"     "Anh"     "AppDec"  "AppInc"  "SleDec"  "SleInc"  "MotoInc"
    ##  [8] "Fatig"   "Guilt"   "Conc"    "Sui"

## Common factor

``` r
all_commonfactor.model <- "
A1 =~ NA*Dep + Anh + AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig + Guilt + Conc + Sui
A1 ~~ 1*A1
"
all_commonfactor.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.715 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
all_commonfactor.fit$modelfit
```

    ##       chisq df      p_chisq      AIC     CFI      SRMR
    ## df 166.7673 44 3.508019e-16 210.7673 0.96518 0.1570351

``` r
all_commonfactor.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 5       A1 =~     Dep   0.83968329 0.0477003658466702
    ## 1       A1 =~     Anh   0.98187134 0.0445271644866201
    ## 2       A1 =~  AppDec   0.05529376 0.0812992696381526
    ## 3       A1 =~  AppInc   0.42227844 0.0713889717404505
    ## 9       A1 =~  SleDec   0.38758192 0.0958970348794045
    ## 10      A1 =~  SleInc   0.51023734 0.0905215116163629
    ## 8       A1 =~ MotoInc   0.20096623  0.121063257517133
    ## 6       A1 =~   Fatig   0.63615156 0.0986747020174307
    ## 7       A1 =~   Guilt   0.59632487 0.0730312292400095
    ## 4       A1 =~    Conc   0.73928338  0.102921749946599
    ## 11      A1 =~     Sui   0.64319048 0.0839617489499526
    ## 17     Dep ~~     Dep   0.29493183 0.0719853203584514
    ## 13     Anh ~~     Anh   0.03592706 0.0632420890079731
    ## 14  AppDec ~~  AppDec   0.99694481  0.197985174075566
    ## 15  AppInc ~~  AppInc   0.82168119  0.134232990857521
    ## 21  SleDec ~~  SleDec   0.84977816   0.22262573856092
    ## 22  SleInc ~~  SleInc   0.73965396  0.233763008400871
    ## 20 MotoInc ~~ MotoInc   0.95961349  0.472026748313016
    ## 18   Fatig ~~   Fatig   0.59531572  0.331399688869406
    ## 19   Guilt ~~   Guilt   0.64439664  0.146238977836729
    ## 16    Conc ~~    Conc   0.45345876  0.276572757572967
    ## 23     Sui ~~     Sui   0.58630609  0.177899446798072
    ## 12      A1 ~~      A1   1.00000000

Correlation between cardinal symptoms

``` r
all_commonfactor_cardinal.model <- "
A1 =~ NA*Dep + Anh + AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig + Guilt + Conc + Sui
A1 ~~ 1*A1
Dep ~~ Anh
"
all_commonfactor_cardinal.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_commonfactor_cardinal.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.828 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_commonfactor_cardinal.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
all_commonfactor_cardinal.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 131.9597 43 5.703567e-11 177.9597 0.9747687 0.1525205

``` r
all_commonfactor_cardinal.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 5       A1 =~     Dep   0.70054970  0.073903151298119
    ## 1       A1 =~     Anh   0.84894801 0.0687373545776568
    ## 2       A1 =~  AppDec   0.04385646 0.0869815492144422
    ## 3       A1 =~  AppInc   0.46776516 0.0784550499840306
    ## 9       A1 =~  SleDec   0.41438865   0.10330858855478
    ## 10      A1 =~  SleInc   0.54586163  0.098407633919925
    ## 8       A1 =~ MotoInc   0.21775270  0.130616766201743
    ## 6       A1 =~   Fatig   0.68768252  0.104526923328342
    ## 7       A1 =~   Guilt   0.65211267 0.0804460344805345
    ## 4       A1 =~    Conc   0.81140087  0.112732265046109
    ## 11      A1 =~     Sui   0.68603568 0.0909916620822564
    ## 17     Dep ~~     Anh   0.30932139  0.101438912687904
    ## 18     Dep ~~     Dep   0.50923055  0.114099114636191
    ## 13     Anh ~~     Anh   0.27928535  0.111102486324011
    ## 14  AppDec ~~  AppDec   0.99807408  0.197864459642927
    ## 15  AppInc ~~  AppInc   0.78119557  0.136401605596288
    ## 22  SleDec ~~  SleDec   0.82828188  0.222798891209103
    ## 23  SleInc ~~  SleInc   0.70203555  0.231379394140284
    ## 21 MotoInc ~~ MotoInc   0.95257859  0.470855918746615
    ## 19   Fatig ~~   Fatig   0.52709381  0.327922505562637
    ## 20   Guilt ~~   Guilt   0.57474862  0.147443191004291
    ## 16    Conc ~~    Conc   0.34162960  0.273384928776974
    ## 24     Sui ~~     Sui   0.52935390  0.177912531886564
    ## 12      A1 ~~      A1   1.00000000

Appetite symptom correlation

``` r
all_commonfactor_app.model <- "
A1 =~ NA*Dep + Anh + AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig + Guilt + Conc + Sui
A1 ~~ 1*A1
Dep ~~ Anh
AppDec ~~ AppInc
"
all_commonfactor_app.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_commonfactor_app.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.679 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_commonfactor_app.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
all_commonfactor_app.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 96.46857 42 3.597387e-06 144.4686 0.9845513 0.1411183

``` r
all_commonfactor_app.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 5       A1 =~     Dep    0.7012332  0.073672388999007
    ## 1       A1 =~     Anh    0.8482954 0.0683461485738379
    ## 2       A1 =~  AppDec    0.1050515 0.0890117272275342
    ## 3       A1 =~  AppInc    0.4723937 0.0783442947208167
    ## 9       A1 =~  SleDec    0.4189116  0.103050987465132
    ## 10      A1 =~  SleInc    0.5444870 0.0981726956626097
    ## 8       A1 =~ MotoInc    0.2219199  0.130497346573179
    ## 6       A1 =~   Fatig    0.6872792  0.104826138428038
    ## 7       A1 =~   Guilt    0.6482263 0.0802746992897671
    ## 4       A1 =~    Conc    0.8096262  0.113056133852223
    ## 11      A1 =~     Sui    0.6853122 0.0913054733918014
    ## 18     Dep ~~     Anh    0.3092002  0.100862739604504
    ## 15  AppDec ~~  AppInc   -0.4752880  0.111375435195644
    ## 19     Dep ~~     Dep    0.5082736  0.113721567238765
    ## 13     Anh ~~     Anh    0.2803963  0.110429151293686
    ## 14  AppDec ~~  AppDec    0.9889648  0.197643880686502
    ## 16  AppInc ~~  AppInc    0.7768443  0.136286967399966
    ## 23  SleDec ~~  SleDec    0.8245135  0.222525258122847
    ## 24  SleInc ~~  SleInc    0.7035350   0.23149521071489
    ## 22 MotoInc ~~ MotoInc    0.9507551  0.470511829287279
    ## 20   Fatig ~~   Fatig    0.5276478  0.327827229413937
    ## 21   Guilt ~~   Guilt    0.5798037  0.147484830707026
    ## 17    Conc ~~    Conc    0.3445058  0.272753280160881
    ## 25     Sui ~~     Sui    0.5303480  0.177916024013584
    ## 12      A1 ~~      A1    1.0000000

Sleep symptom correlation

``` r
all_commonfactor_sle.model <- "
A1 =~ NA*Dep + Anh + AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig + Guilt + Conc + Sui
A1 ~~ 1*A1
Dep ~~ Anh
SleDec ~~ SleInc
"
all_commonfactor_sle.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_commonfactor_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.778 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_commonfactor_sle.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
all_commonfactor_sle.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 137.0488 42 5.216867e-12 185.0488 0.9730416 0.1483226

``` r
all_commonfactor_sle.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 5       A1 =~     Dep    0.6942259 0.0727788227089858
    ## 1       A1 =~     Anh    0.8407633 0.0675842313925317
    ## 2       A1 =~  AppDec    0.0447404 0.0870483309075767
    ## 3       A1 =~  AppInc    0.4690943 0.0782922584520168
    ## 9       A1 =~  SleDec    0.4462266  0.103658858713401
    ## 10      A1 =~  SleInc    0.5696780 0.0990868246227844
    ## 8       A1 =~ MotoInc    0.2198810  0.130679570543012
    ## 6       A1 =~   Fatig    0.6884079   0.10430190827867
    ## 7       A1 =~   Guilt    0.6520736 0.0801165221291648
    ## 4       A1 =~    Conc    0.8140762  0.112497624091042
    ## 11      A1 =~     Sui    0.6845436 0.0905582888905505
    ## 17     Dep ~~     Anh    0.3203724 0.0989076334103791
    ## 23  SleDec ~~  SleInc   -0.2963452  0.152708790149111
    ## 18     Dep ~~     Dep    0.5180496  0.112247113608571
    ## 13     Anh ~~     Anh    0.2931168  0.107925697328708
    ## 14  AppDec ~~  AppDec    0.9979990  0.197853870564062
    ## 15  AppInc ~~  AppInc    0.7799499  0.136319635765898
    ## 22  SleDec ~~  SleDec    0.8008817  0.220565402255493
    ## 24  SleInc ~~  SleInc    0.6754675  0.229488421224686
    ## 21 MotoInc ~~ MotoInc    0.9516522  0.470795610115474
    ## 19   Fatig ~~   Fatig    0.5260950  0.327813046345898
    ## 20   Guilt ~~   Guilt    0.5748023  0.147171083001072
    ## 16    Conc ~~    Conc    0.3372808  0.273197364968731
    ## 25     Sui ~~     Sui    0.5313993  0.177850145092662
    ## 12      A1 ~~      A1    1.0000000

``` r
all_commonfactor_app_sle.model <- "
A1 =~ NA*Dep + Anh + AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig + Guilt + Conc + Sui
A1 ~~ 1*A1
Dep ~~ Anh
AppDec ~~ AppInc
SleDec ~~ SleInc
"
all_commonfactor_app_sle.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_commonfactor_app_sle.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    1.71 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_commonfactor_app_sle.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
all_commonfactor_app_sle.fit$modelfit
```

    ##       chisq df      p_chisq     AIC       CFI     SRMR
    ## df 99.40897 41 9.298351e-07 149.409 0.9834337 0.136471

``` r
all_commonfactor_app_sle.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 5       A1 =~     Dep    0.6947461 0.0725249738617444
    ## 1       A1 =~     Anh    0.8399426 0.0671730830918245
    ## 2       A1 =~  AppDec    0.1063360 0.0891340788283055
    ## 3       A1 =~  AppInc    0.4738463 0.0781782266612554
    ## 9       A1 =~  SleDec    0.4509577  0.103421684396248
    ## 10      A1 =~  SleInc    0.5686700 0.0988236790439758
    ## 8       A1 =~ MotoInc    0.2240895  0.130562639937141
    ## 6       A1 =~   Fatig    0.6880225   0.10459930751095
    ## 7       A1 =~   Guilt    0.6481385 0.0799418208003177
    ## 4       A1 =~    Conc    0.8123050  0.112814097835022
    ## 11      A1 =~     Sui    0.6837762 0.0908636010449618
    ## 18     Dep ~~     Anh    0.3205051 0.0982887356918705
    ## 15  AppDec ~~  AppInc   -0.4760506  0.111462567481901
    ## 24  SleDec ~~  SleInc   -0.2985863  0.152733082296259
    ## 19     Dep ~~     Dep    0.5173261  0.111834430937195
    ## 13     Anh ~~     Anh    0.2944975  0.107213001506005
    ## 14  AppDec ~~  AppDec    0.9886917  0.197600892326293
    ## 16  AppInc ~~  AppInc    0.7754696  0.136200632124996
    ## 23  SleDec ~~  SleDec    0.7966420  0.220261765073314
    ## 25  SleInc ~~  SleInc    0.6766122  0.229575780640633
    ## 22 MotoInc ~~ MotoInc    0.9497810  0.470449386574104
    ## 20   Fatig ~~   Fatig    0.5266259  0.327727679517566
    ## 21   Guilt ~~   Guilt    0.5799160   0.14721103124633
    ## 17    Conc ~~    Conc    0.3401598  0.272566685393269
    ## 26     Sui ~~     Sui    0.5324491  0.177846995964415
    ## 12      A1 ~~      A1    1.0000000

## Two-factor models

``` r
all_psych_soma.model <- "
A1 =~ NA*Dep + Anh + Guilt + Conc + Sui 
A2 =~ NA*AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig
A1 ~~ 1*A1
A2 ~~ 1*A2
Dep ~~ Anh
AppDec ~~ AppInc
"
all_psych_soma.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    3.08 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_psych_soma.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
all_psych_soma.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 96.48378 41 2.282649e-06 146.4838 0.9842633 0.1378455

``` r
all_psych_soma.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 3       A1 =~     Dep    0.7227216 0.0770185534588041
    ## 1       A1 =~     Anh    0.8698566 0.0714893884117358
    ## 4       A1 =~   Guilt    0.6452514 0.0811448471114828
    ## 2       A1 =~    Conc    0.8013119  0.113632777907722
    ## 5       A1 =~     Sui    0.6879581 0.0923473638767151
    ## 8       A2 =~  AppDec   -0.1244413  0.099726828875766
    ## 9       A2 =~  AppInc   -0.5254150 0.0894446023795061
    ## 12      A2 =~  SleDec   -0.4539713  0.115462393596611
    ## 13      A2 =~  SleInc   -0.5977304  0.110480070912601
    ## 11      A2 =~ MotoInc   -0.2628048   0.14405005176426
    ## 10      A2 =~   Fatig   -0.7699839   0.12040870415076
    ## 20     Dep ~~     Anh    0.2753883  0.106715807813112
    ## 17  AppDec ~~  AppInc   -0.4910474  0.114088300552753
    ## 21     Dep ~~     Dep    0.4776733   0.11914057475901
    ## 15     Anh ~~     Anh    0.2433491  0.117552252460424
    ## 23   Guilt ~~   Guilt    0.5836491  0.147883057145116
    ## 19    Conc ~~    Conc    0.3578969  0.272528206398725
    ## 27     Sui ~~     Sui    0.5267131  0.177938180635934
    ## 16  AppDec ~~  AppDec    0.9845135  0.197686202282442
    ## 18  AppInc ~~  AppInc    0.7239384   0.14201151228833
    ## 25  SleDec ~~  SleDec    0.7939115  0.228290406829923
    ## 26  SleInc ~~  SleInc    0.6427162  0.233004861706757
    ## 24 MotoInc ~~ MotoInc    0.9309434  0.469858432012614
    ## 22   Fatig ~~   Fatig    0.4071216  0.344573899052707
    ## 7       A1 ~~      A2   -0.8480733  0.095991030111648
    ## 6       A1 ~~      A1    1.0000000                   
    ## 14      A2 ~~      A2    1.0000000

``` r
all_psych_veg.model <- "
A1 =~ NA*Dep + Anh + Guilt + Sui 
A2 =~ NA*AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig + Conc
A1 ~~ 1*A1
A2 ~~ 1*A2
Dep ~~ Anh
AppDec ~~ AppInc
"
all_psych_veg.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_psych_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.226 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_psych_veg.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
all_psych_veg.fit$modelfit
```

    ##       chisq df     p_chisq      AIC       CFI      SRMR
    ## df 95.82477 41 2.78826e-06 145.8248 0.9844503 0.1361692

``` r
all_psych_veg.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 2       A1 =~     Dep    0.7563642 0.0837608922026265
    ## 1       A1 =~     Anh    0.9036945     0.081154687075
    ## 3       A1 =~   Guilt    0.6384848 0.0821096962222822
    ## 4       A1 =~     Sui    0.6860052 0.0937276013552739
    ## 7       A2 =~  AppDec    0.1208109 0.0995535246134471
    ## 8       A2 =~  AppInc    0.5257788 0.0872225069311045
    ## 12      A2 =~  SleDec    0.4564607  0.113624207029404
    ## 13      A2 =~  SleInc    0.5930646  0.107663712600347
    ## 11      A2 =~ MotoInc    0.2619985  0.142276924051984
    ## 10      A2 =~   Fatig    0.7573010  0.115857539097709
    ## 9       A2 =~    Conc    0.9004399   0.12573580780031
    ## 20     Dep ~~     Anh    0.2205304  0.123939326271934
    ## 17  AppDec ~~  AppInc   -0.4891836  0.114140718817417
    ## 21     Dep ~~     Dep    0.4279141  0.131897347379453
    ## 15     Anh ~~     Anh    0.1833370  0.138502714952531
    ## 23   Guilt ~~   Guilt    0.5923363   0.14772406438054
    ## 27     Sui ~~     Sui    0.5293964  0.178124259393937
    ## 16  AppDec ~~  AppDec    0.9854068  0.197413808685382
    ## 18  AppInc ~~  AppInc    0.7235563  0.142884775371217
    ## 25  SleDec ~~  SleDec    0.7916447  0.225671063737058
    ## 26  SleInc ~~  SleInc    0.6482707  0.234805423352101
    ## 24 MotoInc ~~ MotoInc    0.9313557  0.470067946544218
    ## 22   Fatig ~~   Fatig    0.4264971  0.330303344351807
    ## 19    Conc ~~    Conc    0.1892088  0.290671674640801
    ## 6       A1 ~~      A2    0.8079316 0.0883861355804176
    ## 5       A1 ~~      A1    1.0000000                   
    ## 14      A2 ~~      A2    1.0000000

``` r
all_affect_veg.model <- "
A1 =~ NA*Dep + Guilt + Sui 
A2 =~ NA*Anh + AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig + Conc
A1 ~~ 1*A1
A2 ~~ 1*A2
Dep ~~ Anh
AppDec ~~ AppInc
"
all_affect_veg.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_affect_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   2.045 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_affect_veg.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
all_affect_veg.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI     SRMR
    ## df 79.01877 41 0.0003325273 129.0188 0.9892169 0.133106

``` r
all_affect_veg.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 1       A1 =~     Dep    0.7507985 0.0761037224593046
    ## 2       A1 =~   Guilt    0.7308977 0.0870950857627854
    ## 3       A1 =~     Sui    0.7913855  0.102736619277151
    ## 6       A2 =~     Anh    0.8635945 0.0686471311540735
    ## 7       A2 =~  AppDec    0.1116003 0.0939315422013412
    ## 8       A2 =~  AppInc    0.5050069 0.0813264280439135
    ## 12      A2 =~  SleDec    0.4312126  0.106971433089306
    ## 13      A2 =~  SleInc    0.5758865  0.101183554941139
    ## 11      A2 =~ MotoInc    0.2427216  0.135493669785797
    ## 10      A2 =~   Fatig    0.7299895  0.107283529566791
    ## 9       A2 =~    Conc    0.8589718  0.116630052552988
    ## 20     Dep ~~     Anh    0.4179745 0.0901040557119042
    ## 17  AppDec ~~  AppInc   -0.4820232  0.112740383257898
    ## 21     Dep ~~     Dep    0.4363016  0.123042867508004
    ## 23   Guilt ~~   Guilt    0.4657892  0.157142278994902
    ## 27     Sui ~~     Sui    0.3737065   0.18831483766733
    ## 15     Anh ~~     Anh    0.2542039  0.113380737843865
    ## 16  AppDec ~~  AppDec    0.9875454  0.197486626208013
    ## 18  AppInc ~~  AppInc    0.7449682  0.138280568027798
    ## 25  SleDec ~~  SleDec    0.8140570  0.222533669851174
    ## 26  SleInc ~~  SleInc    0.6683521  0.231856224992288
    ## 24 MotoInc ~~ MotoInc    0.9410846  0.470607398435132
    ## 22   Fatig ~~   Fatig    0.4671156  0.323254983910501
    ## 19    Conc ~~    Conc    0.2621648  0.274599196778721
    ## 5       A1 ~~      A2    0.7496738 0.0587062232332387
    ## 4       A1 ~~      A1    1.0000000                   
    ## 14      A2 ~~      A2    1.0000000

### Three factor model

``` r
all_cog_mood_neuroveg.model <- "
A1 =~ NA*Guilt + Conc + Sui 
A2 =~ NA*Dep + Anh + Guilt
A3 =~ NA*AppDec + AppInc + SleDec + SleInc + MotoInc + Fatig
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
AppDec ~~ AppInc
c2 > 0.001
Anh ~~ c2*Anh
a12 < 0.99
A1 ~~ a12*A2
"
all_cog_mood_neuroveg.fit <- usermodel(all_symptoms_covstruct, estimation='DWLS', model=all_cog_mood_neuroveg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  25.012 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936793 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(all_symptoms_covstruct, estimation = "DWLS", model =
    ## all_cog_mood_neuroveg.model): A difference greater than .025 was observed pre-
    ## and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
all_cog_mood_neuroveg.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 103.2905 39 1.011332e-07 157.2905 0.9817655 0.1381106

``` r
all_cog_mood_neuroveg.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs  STD_Genotype    STD_Genotype_SE
    ## 2       A1 =~   Guilt  1.1104558036  0.826930807896118
    ## 1       A1 =~    Conc  0.7775756237  0.123781170701285
    ## 3       A1 =~     Sui  0.6827698338  0.100081897122114
    ## 8       A2 =~     Dep  0.8766076732 0.0472773117098039
    ## 7       A2 =~     Anh  1.0105356493 0.0474178691555137
    ## 9       A2 =~   Guilt -0.4430811203  0.837878731240534
    ## 12      A3 =~  AppDec -0.1246240407 0.0994000988344907
    ## 13      A3 =~  AppInc -0.5245667680 0.0892746611471191
    ## 16      A3 =~  SleDec -0.4542751369  0.115857792635345
    ## 17      A3 =~  SleInc -0.5978144459   0.11034359035845
    ## 15      A3 =~ MotoInc -0.2631201776  0.144032778865397
    ## 14      A3 =~   Fatig -0.7708765585  0.120941368263306
    ## 21  AppDec ~~  AppInc -0.4910377311  0.113983698165299
    ## 19     Anh ~~     Anh  0.0009999897 0.0669219380344484
    ## 5       A1 ~~      A2  0.8726148929   0.11448391980931
    ## 26   Guilt ~~   Guilt  0.4292607181  0.287877477983326
    ## 23    Conc ~~    Conc  0.3953739148  0.273688254551566
    ## 30     Sui ~~     Sui  0.5338233354  0.184324967712537
    ## 24     Dep ~~     Dep  0.2315584325 0.0705861917038074
    ## 20  AppDec ~~  AppDec  0.9844688402  0.197743800779799
    ## 22  AppInc ~~  AppInc  0.7248293246  0.142048241124291
    ## 28  SleDec ~~  SleDec  0.7936329258  0.228751456477724
    ## 29  SleInc ~~  SleInc  0.6426193506  0.233229782470779
    ## 27 MotoInc ~~ MotoInc  0.9307672525  0.469699259347056
    ## 25   Fatig ~~   Fatig  0.4057509932  0.345678020363645
    ## 6       A1 ~~      A3 -0.8019369615  0.103673268788238
    ## 11      A2 ~~      A3 -0.7242639328 0.0982925613996452
    ## 4       A1 ~~      A1  1.0000000000                   
    ## 10      A2 ~~      A2  1.0000000000                   
    ## 18      A3 ~~      A3  1.0000000000

### Model comparisons

``` r
all_model_list=
list("1a"=list(name="Common", model=all_commonfactor.fit),
     "1b"=list(name="Common (cardinal)", model=all_commonfactor_cardinal.fit),
     "1c"=list(name="Common (App)", model=all_commonfactor_app.fit),
     "1d"=list(name="Common (Sle)", model=all_commonfactor_sle.fit),
     "1er"=list(name="Common (App,Sle)", model=all_commonfactor_app_sle.fit),
     "2a(i)"=list(name="Psych-Somatic", model=all_psych_soma.fit),
     "2b(i)"=list(name="Psych-Neuroveg", model=all_psych_veg.fit),
     "2c(i)"=list(name="Affect-Neuroveg", model=all_affect_veg.fit),
     "3"=list(name="Cog-Mood-Neuroveg", model=all_cog_mood_neuroveg.fit)
     )

all_model_fits <- 
data.frame(Model=names(all_model_list),
           Name=sapply(all_model_list, function(m) m$name)) %>%
bind_cols(
bind_rows(
lapply(all_model_list, function(m) m$model$modelfit)
))
rownames(all_model_fits) <- NULL

knitr::kable(
all_model_fits %>%
mutate(dAIC=AIC-min(AIC)) %>%
mutate_if(is.numeric, ~signif(., 3))
)
```

| Model | Name              | chisq |  df |  p_chisq | AIC |   CFI |  SRMR | dAIC |
|:------|:------------------|------:|----:|---------:|----:|------:|------:|-----:|
| 1a    | Common            | 167.0 |  44 | 0.00e+00 | 211 | 0.965 | 0.157 | 81.7 |
| 1b    | Common (cardinal) | 132.0 |  43 | 0.00e+00 | 178 | 0.975 | 0.153 | 48.9 |
| 1c    | Common (App)      |  96.5 |  42 | 3.60e-06 | 144 | 0.985 | 0.141 | 15.4 |
| 1d    | Common (Sle)      | 137.0 |  42 | 0.00e+00 | 185 | 0.973 | 0.148 | 56.0 |
| 1er   | Common (App,Sle)  |  99.4 |  41 | 9.00e-07 | 149 | 0.983 | 0.136 | 20.4 |
| 2a(i) | Psych-Somatic     |  96.5 |  41 | 2.30e-06 | 146 | 0.984 | 0.138 | 17.5 |
| 2b(i) | Psych-Neuroveg    |  95.8 |  41 | 2.80e-06 | 146 | 0.984 | 0.136 | 16.8 |
| 2c(i) | Affect-Neuroveg   |  79.0 |  41 | 3.33e-04 | 129 | 0.989 | 0.133 |  0.0 |
| 3     | Cog-Mood-Neuroveg | 103.0 |  39 | 1.00e-07 | 157 | 0.982 | 0.138 | 28.3 |
