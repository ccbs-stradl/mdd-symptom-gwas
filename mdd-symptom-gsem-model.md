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
```

    ## Loading required package: readr

``` r
require(tidyr)
```

    ## Loading required package: tidyr

``` r
require(stringr)
```

    ## Loading required package: stringr

``` r
require(dplyr)
```

    ## Loading required package: dplyr

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
require(ggplot2)
```

    ## Loading required package: ggplot2

``` r
require(corrplot)
```

    ## Loading required package: corrplot

    ## corrplot 0.92 loaded

``` r
require(GenomicSEM)
```

    ## Loading required package: GenomicSEM

``` r
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
    ## ── Column specification ──────────────────────────────────────────────────
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
    ## ── Column specification ──────────────────────────────────────────────────
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
    ## ── Column specification ──────────────────────────────────────────────────
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

    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec 
    ## 0.128871019 0.051543416 0.013663704 0.013924961 0.039896861 0.002375007 
    ##     ClinSui      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec 
    ## 0.103892914 0.082198484 0.087343514 0.035741853 0.079200630 0.043111949 
    ##   PopSleInc  PopMotoInc  PopMotoDec    PopFatig    PopGuilt     PopConc 
    ## 0.047553540 0.369275531 0.079699383 0.057539570 0.060012455 0.056887041 
    ##      PopSui 
    ## 0.033166683

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
    ##   3.219 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0464961771558631 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65818091523924 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq  df      p_chisq      AIC       CFI      SRMR
    ## df 363.6813 104 1.670625e-30 427.6813 0.9446109 0.1804122

``` r
commonfactor.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE
    ## 1           A1 =~  ClinAppDec   0.08130171 0.0764916664747164
    ## 2           A1 =~  ClinAppInc  -0.14987805   0.10623193331763
    ## 4           A1 =~  ClinSleDec  -0.12195857  0.134137517511587
    ## 5           A1 =~  ClinSleInc  -0.19692788  0.145262383971543
    ## 3           A1 =~ ClinMotoInc  -0.20654637  0.112346761940068
    ## 6           A1 =~     ClinSui  -0.57627045 0.0983562590908127
    ## 11          A1 =~      PopDep  -0.83888334 0.0508263648701571
    ## 7           A1 =~      PopAnh  -0.94442845 0.0448495363572952
    ## 8           A1 =~   PopAppDec  -0.17696093 0.0823865158386624
    ## 9           A1 =~   PopAppInc  -0.39015184  0.074879044752346
    ## 14          A1 =~   PopSleDec  -0.56377208 0.0936828010541065
    ## 15          A1 =~   PopSleInc  -0.45951001 0.0986809312314969
    ## 12          A1 =~    PopFatig  -0.66439698 0.0920345866279589
    ## 13          A1 =~    PopGuilt  -0.63851961 0.0785540781848866
    ## 10          A1 =~     PopConc  -0.70426809 0.0932239630742708
    ## 16          A1 =~      PopSui  -0.59025474 0.0961530135389394
    ## 18  ClinAppDec ~~  ClinAppDec   0.99339007  0.202993435241858
    ## 19  ClinAppInc ~~  ClinAppInc   0.97753717  0.310498997391932
    ## 21  ClinSleDec ~~  ClinSleDec   0.98512633  0.597177281731245
    ## 22  ClinSleInc ~~  ClinSleInc   0.96122099  0.675737153113016
    ## 20 ClinMotoInc ~~ ClinMotoInc   0.95733748  0.433493275444151
    ## 23     ClinSui ~~     ClinSui   0.66791234  0.303367534581501
    ## 28      PopDep ~~      PopDep   0.29627443 0.0767317038170107
    ## 24      PopAnh ~~      PopAnh   0.10805436 0.0626318763425247
    ## 25   PopAppDec ~~   PopAppDec   0.96868191  0.233107614222654
    ## 26   PopAppInc ~~   PopAppInc   0.84778169  0.148689919670437
    ## 31   PopSleDec ~~   PopSleDec   0.68216132  0.277925670222859
    ## 32   PopSleInc ~~   PopSleInc   0.78885096  0.242806137543823
    ## 29    PopFatig ~~    PopFatig   0.55857657  0.306608227771798
    ## 30    PopGuilt ~~    PopGuilt   0.59229268  0.157753428523502
    ## 27     PopConc ~~     PopConc   0.50400584  0.240521098221822
    ## 33      PopSui ~~      PopSui   0.65159930  0.239092393259687
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
    ##    2.41 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0464961771558631 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65818091523924 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 257.097 99 5.327374e-16 331.097 0.9662785 0.1661774

``` r
commonfactor_dir.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppDec   0.09787682 0.0820620181918304 2.329810e-01
    ## 2           A1 =~  ClinAppInc  -0.17948958  0.114057352058428 1.155619e-01
    ## 4           A1 =~  ClinSleDec  -0.12775538  0.142220005109015 3.690286e-01
    ## 5           A1 =~  ClinSleInc  -0.22105198  0.153145939106859 1.489054e-01
    ## 3           A1 =~ ClinMotoInc  -0.22854008  0.120381916178802 5.763662e-02
    ## 6           A1 =~     ClinSui  -0.63224839  0.107031225605526 3.481083e-09
    ## 11          A1 =~      PopDep  -0.67858602 0.0699808160885343 3.112981e-22
    ## 7           A1 =~      PopAnh  -0.79410567 0.0639194877282234 1.947722e-35
    ## 8           A1 =~   PopAppDec  -0.21074650 0.0898762663390702 1.903460e-02
    ## 9           A1 =~   PopAppInc  -0.43407398 0.0804899559093647 6.933192e-08
    ## 14          A1 =~   PopSleDec  -0.63900573  0.102878592304825 5.256464e-10
    ## 15          A1 =~   PopSleInc  -0.50520133   0.10709853956406 2.391584e-06
    ## 12          A1 =~    PopFatig  -0.72805135  0.098187695314417 1.217282e-13
    ## 13          A1 =~    PopGuilt  -0.69354033 0.0858776382636896 6.697275e-16
    ## 10          A1 =~     PopConc  -0.77098447   0.10109517117077 2.415483e-14
    ## 16          A1 =~      PopSui  -0.61887816  0.102387945045138 1.499282e-09
    ## 19  ClinAppDec ~~  ClinAppInc  -0.40823269  0.191176715269992 3.273181e-02
    ## 23  ClinSleDec ~~  ClinSleInc   0.36750371  0.462415873425517 4.267583e-01
    ## 31      PopDep ~~      PopAnh   0.36708113  0.090865523134496 5.349441e-05
    ## 28   PopAppDec ~~   PopAppInc  -0.23286661  0.132991476877965 7.994749e-02
    ## 36   PopSleDec ~~   PopSleInc  -0.34384259  0.190386371048206 7.091651e-02
    ## 18  ClinAppDec ~~  ClinAppDec   0.99041889  0.203295804911015 1.105762e-06
    ## 20  ClinAppInc ~~  ClinAppInc   0.96778215  0.310637043564592 1.836449e-03
    ## 22  ClinSleDec ~~  ClinSleDec   0.98367529  0.596337384711498 9.903839e-02
    ## 24  ClinSleInc ~~  ClinSleInc   0.95113147  0.675042448494817 1.588346e-01
    ## 21 ClinMotoInc ~~ ClinMotoInc   0.94776830  0.432877618596714 2.856325e-02
    ## 25     ClinSui ~~     ClinSui   0.60026198  0.311074718524763 5.365090e-02
    ## 32      PopDep ~~      PopDep   0.53952162  0.109496603550713 8.337941e-07
    ## 26      PopAnh ~~      PopAnh   0.36939641  0.096574162402597 1.307802e-04
    ## 27   PopAppDec ~~   PopAppDec   0.95558856  0.231794445460722 3.746889e-05
    ## 29   PopAppInc ~~   PopAppInc   0.81157811  0.150384849074252 6.788226e-08
    ## 35   PopSleDec ~~   PopSleDec   0.59167344  0.281427762143967 3.551841e-02
    ## 37   PopSleInc ~~   PopSleInc   0.74477335  0.245182923611921 2.384558e-03
    ## 33    PopFatig ~~    PopFatig   0.46994205  0.310919803401855 1.306718e-01
    ## 34    PopGuilt ~~    PopGuilt   0.51900370  0.163489998255946 1.500883e-03
    ## 30     PopConc ~~     PopConc   0.40558319  0.242135429153858 9.392859e-02
    ## 38      PopSui ~~      PopSui   0.61699165  0.242096692576092 1.081786e-02
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
    ## 1           A1 =~  ClinAppDec             -0.03617340 0.033750220
    ## 2           A1 =~  ClinAppInc              0.05193200 0.038488071
    ## 3           A1 =~  ClinSleDec              0.03133530 0.037444607
    ## 4           A1 =~  ClinSleInc              0.03916668 0.030126554
    ## 5           A1 =~ ClinMotoInc              0.04669981 0.031818780
    ## 6           A1 =~     ClinSui              0.21065315 0.099165880
    ## 7           A2 =~      PopDep              0.19694228 0.020281036
    ## 8           A2 =~      PopAnh              0.23672024 0.019065796
    ## 9           A2 =~   PopAppDec              0.04161314 0.017746622
    ## 10          A2 =~   PopAppInc              0.12665671 0.023484907
    ## 11          A2 =~   PopSleDec              0.14472655 0.023300614
    ## 12          A2 =~   PopSleInc              0.11645842 0.024689509
    ## 13          A2 =~    PopFatig              0.18392380 0.024805661
    ## 14          A2 =~    PopGuilt              0.17088600 0.021159727
    ## 15          A2 =~     PopConc              0.19622264 0.025730370
    ## 16          A2 =~      PopSui              0.11708572 0.019369401
    ## 19          A1 ~~          A2              1.00000006 0.467580181
    ## 20  ClinAppDec ~~  ClinAppInc             -0.04365280 0.020467761
    ## 21  ClinSleDec ~~  ClinSleInc              0.01597122 0.019997051
    ## 22      PopDep ~~      PopAnh              0.03175803 0.007853235
    ## 23   PopAppDec ~~   PopAppInc             -0.01341658 0.007662260
    ## 24   PopSleDec ~~   PopSleInc             -0.01795177 0.009939994
    ## 25  ClinAppDec ~~  ClinAppDec              0.13528224 0.027897671
    ## 26  ClinAppInc ~~  ClinAppInc              0.08101577 0.025980404
    ## 27  ClinSleDec ~~  ClinSleDec              0.05917797 0.035901785
    ## 28  ClinSleInc ~~  ClinSleInc              0.02985965 0.021213604
    ## 29 ClinMotoInc ~~ ClinMotoInc              0.03957413 0.018148706
    ## 30     ClinSui ~~     ClinSui              0.06663421 0.048126000
    ## 31      PopDep ~~      PopDep              0.04544420 0.009209213
    ## 32      PopAnh ~~      PopAnh              0.03282517 0.008585491
    ## 33   PopAppDec ~~   PopAppDec              0.03725711 0.009037381
    ## 34   PopAppInc ~~   PopAppInc              0.06909717 0.012803495
    ## 35   PopSleDec ~~   PopSleDec              0.03035077 0.014436268
    ## 36   PopSleInc ~~   PopSleInc              0.03957643 0.013029026
    ## 37    PopFatig ~~    PopFatig              0.02999135 0.019842991
    ## 38    PopGuilt ~~    PopGuilt              0.03150917 0.009924676
    ## 39     PopConc ~~     PopConc              0.02627157 0.015683808
    ## 40      PopSui ~~      PopSui              0.02208390 0.008665682

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
    ## 1           A1 =~  ClinAppDec           -2.020613e+01
    ## 2           A1 =~  ClinAppInc            1.907050e-03
    ## 3           A1 =~  ClinSleDec           -2.712678e-03
    ## 4           A1 =~  ClinSleInc           -1.259496e-03
    ## 5           A1 =~ ClinMotoInc           -6.061003e-04
    ## 6           A1 =~     ClinSui           -6.477481e-04
    ## 7           A2 =~      PopDep            2.546671e-01
    ## 8           A2 =~      PopAnh            2.524779e-01
    ## 9           A2 =~   PopAppDec            2.793639e-02
    ## 10          A2 =~   PopAppInc            3.595431e-02
    ## 11          A2 =~   PopSleDec            6.063403e-02
    ## 12          A2 =~   PopSleInc            8.831097e-02
    ## 13          A2 =~    PopFatig            6.681558e-02
    ## 14          A2 =~    PopGuilt            7.351482e-02
    ## 15          A2 =~     PopConc            8.578286e-02
    ## 16          A2 =~      PopSui            8.517282e-02
    ## 17           A =~  ClinAppDec            8.010776e-02
    ## 18           A =~  ClinAppInc           -8.730696e-02
    ## 19           A =~  ClinSleDec           -4.458396e-02
    ## 20           A =~  ClinSleInc           -5.408140e-02
    ## 21           A =~ ClinMotoInc           -6.233179e-02
    ## 22           A =~     ClinSui           -2.817111e-01
    ## 23           A =~      PopDep           -9.935318e-02
    ## 24           A =~      PopAnh           -1.454918e-01
    ## 25           A =~   PopAppDec           -1.834164e-02
    ## 26           A =~   PopAppInc           -1.357401e-01
    ## 27           A =~   PopSleDec           -1.218829e-01
    ## 28           A =~   PopSleInc           -5.754565e-02
    ## 29           A =~    PopFatig           -1.783349e-01
    ## 30           A =~    PopGuilt           -1.535030e-01
    ## 31           A =~     PopConc           -1.723882e-01
    ## 32           A =~      PopSui           -6.962289e-02
    ## 39  ClinAppDec ~~  ClinAppDec           -4.081575e+02
    ## 40  ClinAppInc ~~  ClinAppInc            7.608654e-02
    ## 41  ClinSleDec ~~  ClinSleDec            5.816470e-02
    ## 42  ClinSleInc ~~  ClinSleInc            2.846736e-02
    ## 43 ClinMotoInc ~~ ClinMotoInc            3.786940e-02
    ## 44     ClinSui ~~     ClinSui            3.164763e-02
    ## 45      PopDep ~~      PopDep            9.504033e-03
    ## 46      PopAnh ~~      PopAnh            3.948709e-03
    ## 47   PopAppDec ~~   PopAppDec            3.787190e-02
    ## 48   PopAppInc ~~   PopAppInc            6.542104e-02
    ## 49   PopSleDec ~~   PopSleDec            3.276460e-02
    ## 50   PopSleInc ~~   PopSleInc            4.202866e-02
    ## 51    PopFatig ~~    PopFatig            2.755169e-02
    ## 52    PopGuilt ~~    PopGuilt            3.174359e-02
    ## 53     PopConc ~~     PopConc            2.769851e-02
    ## 54      PopSui ~~      PopSui            2.369121e-02

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
    ## 1           A1 =~  ClinAppDec            1.724879e+01
    ## 2           A1 =~  ClinAppInc           -2.471221e-03
    ## 3           A1 =~  ClinSleDec            2.765596e-03
    ## 4           A1 =~  ClinSleInc            1.449400e-03
    ## 5           A1 =~ ClinMotoInc            4.904932e-04
    ## 6           A1 =~     ClinSui           -5.325391e-04
    ## 8   ClinAppDec ~~  ClinAppDec           -2.973817e+02
    ## 9   ClinAppInc ~~  ClinAppInc            7.672832e-02
    ## 10  ClinSleDec ~~  ClinSleDec            5.808735e-02
    ## 11  ClinSleInc ~~  ClinSleInc            2.475221e-02
    ## 12 ClinMotoInc ~~ ClinMotoInc            4.009843e-02
    ## 13     ClinSui ~~     ClinSui            1.065372e-01

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
    ##   6.331 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0444310362906712 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.57437590660252 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC CFI     SRMR
    ## df 3.430553  8 0.9045121 29.43055   1 0.133462

``` r
clin_commonfactor_app.fit$results[c(1,2,3,6,7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppDec  0.651426548 0.408626630840785 0.1108950316
    ## 2           A1 =~  ClinAppInc  0.536484445 0.386448228377997 0.1650629155
    ## 4           A1 =~  ClinSleDec  0.865631533 0.510931642461104 0.0902236180
    ## 5           A1 =~  ClinSleInc  0.452433501 0.376628473983664 0.2296459515
    ## 3           A1 =~ ClinMotoInc  0.263772472 0.231860917290513 0.2552741536
    ## 6           A1 =~     ClinSui  0.005897957 0.181448498705048 0.9740698148
    ## 9   ClinAppDec ~~  ClinAppInc -0.762487927 0.424175649580971 0.0722449810
    ## 8   ClinAppDec ~~  ClinAppDec  0.575643339 0.534308319861908 0.2813222365
    ## 10  ClinAppInc ~~  ClinAppInc  0.712184403 0.503362788749754 0.1571138951
    ## 12  ClinSleDec ~~  ClinSleDec  0.250681930  1.10235663672449 0.8201059699
    ## 13  ClinSleInc ~~  ClinSleInc  0.795304209 0.924216842894587 0.3895038141
    ## 11 ClinMotoInc ~~ ClinMotoInc  0.930424290 0.452208203923222 0.0396369184
    ## 14     ClinSui ~~     ClinSui  0.999965284 0.299516373951668 0.0008420187
    ## 7           A1 ~~          A1  1.000000000                             NA

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
    ## 1           A1 =~  ClinAppDec            3.016112e-03
    ## 2           A1 =~  ClinAppInc            1.523429e-03
    ## 3           A1 =~  ClinSleDec            1.740891e+01
    ## 4           A1 =~  ClinSleInc            6.076210e+00
    ## 5           A1 =~ ClinMotoInc            4.538184e-04
    ## 6           A1 =~     ClinSui            7.814385e-05
    ## 8   ClinSleDec ~~  ClinSleInc           -1.057616e+02
    ## 9   ClinAppDec ~~  ClinAppDec            1.388389e-01
    ## 10  ClinAppInc ~~  ClinAppInc            7.673242e-02
    ## 11  ClinSleDec ~~  ClinSleDec           -3.030121e+02
    ## 12  ClinSleInc ~~  ClinSleInc           -3.689557e+01
    ## 13 ClinMotoInc ~~ ClinMotoInc            4.009833e-02
    ## 14     ClinSui ~~     ClinSui            1.065372e-01

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
    ##    1.22 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0444310362906712 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.57437590660252 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df  p_chisq      AIC CFI     SRMR
    ## df 3.430554  8 0.904512 29.43055   1 0.133462

``` r
clin_soma.fit$results[c(1,2,3,6,7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE     p_value
    ## 1           A1 =~  ClinAppDec  0.651426329 0.408626547730766 0.110893947
    ## 2           A1 =~  ClinAppInc  0.536484607 0.386448292216482 0.165061979
    ## 4           A1 =~  ClinSleDec  0.865631382 0.510931570136897 0.090222444
    ## 5           A1 =~  ClinSleInc  0.452433494 0.376628512859596 0.229644673
    ## 3           A1 =~ ClinMotoInc  0.263772572 0.231860950302295 0.255273547
    ## 11  ClinAppDec ~~  ClinAppInc -0.762488046 0.424175631404987 0.072243344
    ## 10  ClinAppDec ~~  ClinAppDec  0.575643753  0.53430805137269 0.281313392
    ## 12  ClinAppInc ~~  ClinAppInc  0.712184310 0.503362933650043 0.157110074
    ## 14  ClinSleDec ~~  ClinSleDec  0.250682361   1.1023564165387 0.820111011
    ## 15  ClinSleInc ~~  ClinSleInc  0.795304060 0.924216843732358 0.389505302
    ## 13 ClinMotoInc ~~ ClinMotoInc  0.930423979 0.452208220848423 0.039637092
    ## 9           A2 ~~          A2  1.000000000 0.299544527725393 0.000842618
    ## 7           A1 ~~          A2  0.005897855 0.181448523336873 0.974068388
    ## 8           A2 =~     ClinSui  1.000000000                            NA
    ## 6           A1 ~~          A1  1.000000000                            NA
    ## 16     ClinSui ~~     ClinSui  0.000000000                            NA

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
    ##   2.626 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0444310362906712 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.57437590660252 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 2.993468  7 0.8856059 30.99347   1 0.1327096

``` r
clin_soma_app.fit$results[c(1,2,3,6,7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppInc   0.61307201 0.439318567812275 1.628627e-01
    ## 3           A1 =~  ClinSleDec   0.74350053 0.460346511906823 1.062905e-01
    ## 4           A1 =~  ClinSleInc   0.45618263 0.378391567245566 2.279770e-01
    ## 2           A1 =~ ClinMotoInc   0.25704956 0.229543995027097 2.627872e-01
    ## 14  ClinAppInc ~~  ClinAppDec  -0.86102374 0.541862200852923 1.120568e-01
    ## 17  ClinSleDec ~~  ClinSleDec   0.44720695    0.934845448277 6.323864e-01
    ## 15  ClinAppInc ~~  ClinAppInc   0.62414270 0.599209284233519 2.975869e-01
    ## 18  ClinSleInc ~~  ClinSleInc   0.79189828 0.920632623649619 3.896982e-01
    ## 16 ClinMotoInc ~~ ClinMotoInc   0.93392679 0.450733542245363 3.826421e-02
    ## 9           A2 ~~          A2   1.00000000 0.199671823846446 5.493632e-07
    ## 12          A3 ~~          A3   1.00000000 0.299544527725393 8.426180e-04
    ## 6           A1 ~~          A2   0.73077209 0.482144497972319 1.296017e-01
    ## 7           A1 ~~          A3   0.14826964 0.245209368938142 5.454023e-01
    ## 10          A2 ~~          A3  -0.07551224 0.165563683620563 6.483244e-01
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
    ##   1.522 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 59.69967 35 0.005743458 99.69967 0.9912037 0.1223695

``` r
pop_commonfactor.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep  0.864944683 0.0510571525602896
    ## 1         A1 =~    PopAnh  0.997506447 0.0455150755405114
    ## 2         A1 =~ PopAppDec  0.160767623 0.0874865109207056
    ## 3         A1 =~ PopAppInc  0.404184077  0.077704024575744
    ## 8         A1 =~ PopSleDec  0.576611981 0.0993927475884935
    ## 9         A1 =~ PopSleInc  0.518572263  0.104036420179739
    ## 6         A1 =~  PopFatig  0.663615266 0.0967892107997069
    ## 7         A1 =~  PopGuilt  0.610657296 0.0787569509442939
    ## 4         A1 =~   PopConc  0.702927927 0.0967391386767126
    ## 10        A1 =~    PopSui  0.588273334 0.0995834868923105
    ## 16    PopDep ~~    PopDep  0.251870964 0.0798577415104395
    ## 12    PopAnh ~~    PopAnh  0.004981294  0.066442148728655
    ## 13 PopAppDec ~~ PopAppDec  0.974154168  0.252205164502798
    ## 14 PopAppInc ~~ PopAppInc  0.836635229  0.160713249832841
    ## 19 PopSleDec ~~ PopSleDec  0.667516825  0.314946273190408
    ## 20 PopSleInc ~~ PopSleInc  0.731083792  0.271128952490912
    ## 17  PopFatig ~~  PopFatig  0.559615003  0.336858183975526
    ## 18  PopGuilt ~~  PopGuilt  0.627097365  0.157475590755255
    ## 15   PopConc ~~   PopConc  0.505892560   0.26054794277742
    ## 21    PopSui ~~    PopSui  0.653935350  0.256498185441794
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
    ##   1.505 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 43.66149 34 0.1239931 85.66149 0.9965593 0.1159836

``` r
pop_commonfactor_gating.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep    0.7463334 0.0839975899254058
    ## 1         A1 =~    PopAnh    0.8818492 0.0759194430382935
    ## 2         A1 =~ PopAppDec    0.1728245 0.0932410194722276
    ## 3         A1 =~ PopAppInc    0.4342164  0.084606416319221
    ## 8         A1 =~ PopSleDec    0.6253953   0.10979510675781
    ## 9         A1 =~ PopSleInc    0.5429831    0.1118825815891
    ## 6         A1 =~  PopFatig    0.7150320  0.104488421438397
    ## 7         A1 =~  PopGuilt    0.6577851 0.0868886411197645
    ## 4         A1 =~   PopConc    0.7662175  0.107488035500716
    ## 10        A1 =~    PopSui    0.6211550  0.107270317973983
    ## 16    PopDep ~~    PopAnh    0.2654812  0.119730914725366
    ## 17    PopDep ~~    PopDep    0.4429866  0.137364143197446
    ## 12    PopAnh ~~    PopAnh    0.2223423  0.126261127882081
    ## 13 PopAppDec ~~ PopAppDec    0.9701328  0.251757679114272
    ## 14 PopAppInc ~~ PopAppInc    0.8114556  0.163711783238023
    ## 20 PopSleDec ~~ PopSleDec    0.6088807  0.317353372002421
    ## 21 PopSleInc ~~ PopSleInc    0.7051707  0.273797348773763
    ## 18  PopFatig ~~  PopFatig    0.4887291  0.338948062281623
    ## 19  PopGuilt ~~  PopGuilt    0.5673188  0.161401068955136
    ## 15   PopConc ~~   PopConc    0.4129112  0.260478549039582
    ## 22    PopSui ~~    PopSui    0.6141676  0.259938341560857
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
    ##   1.583 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 42.73627 33 0.1194377 86.73627 0.9965326 0.1138434

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
    ##   1.625 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##      chisq df   p_chisq     AIC       CFI      SRMR
    ## df 39.0324 33 0.2169376 83.0324 0.9978517 0.1030787

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
    ##  10.248 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 37.49323 32 0.2317148 83.49323 0.9980437 0.1006463

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
    ## 1         A1 =~  PopGuilt             0.139602397 0.094968263
    ## 2         A1 =~   PopConc             0.181956477 0.032196278
    ## 3         A1 =~    PopSui             0.109178167 0.022029014
    ## 4         A2 =~    PopDep             0.253300332 0.014347025
    ## 5         A2 =~    PopAnh             0.308918312 0.014452585
    ## 6         A2 =~  PopGuilt             0.014951296 0.095779083
    ## 7         A3 =~ PopSleDec             0.122285238 0.026506226
    ## 8         A3 =~ PopSleInc             0.108694023 0.026026744
    ## 9         A3 =~  PopFatig             0.156039709 0.030550000
    ## 10        A3 =~ PopAppDec             0.029833953 0.016369917
    ## 11        A3 =~ PopAppInc             0.112260510 0.027145350
    ## 15  PopGuilt ~~  PopGuilt             0.036721331 0.010487134
    ## 16   PopConc ~~   PopConc             0.026014771 0.016989310
    ## 17    PopSui ~~    PopSui             0.021358100 0.008691294
    ## 18    PopDep ~~    PopDep             0.018074624 0.006246406
    ## 19    PopAnh ~~    PopAnh            -0.008079386 0.006910021
    ## 20 PopSleDec ~~ PopSleDec             0.030432149 0.014489360
    ## 21 PopSleInc ~~ PopSleInc             0.036055616 0.013237605
    ## 22  PopFatig ~~  PopFatig             0.033431753 0.020360127
    ## 23 PopAppDec ~~ PopAppDec             0.035139011 0.009092321
    ## 24 PopAppInc ~~ PopAppInc             0.066640761 0.013269944
    ## 25        A1 ~~        A2             0.861152567 0.143350168
    ## 26        A1 ~~        A3             1.228290053 0.246182452
    ## 27        A2 ~~        A3             0.932383684 0.157471572

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
    ##  20.251 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##      chisq df   p_chisq     AIC       CFI      SRMR
    ## df 50.2564 31 0.0157931 98.2564 0.9931422 0.1164101

``` r
pop_cog_mood_neuroveg_constr.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 2         A1 =~  PopGuilt  0.797874584  0.673630969467514
    ## 1         A1 =~   PopConc  0.768102342  0.133052834980531
    ## 3         A1 =~    PopSui  0.628045950  0.123429775631199
    ## 8         A2 =~    PopDep  0.894729284 0.0503963642394777
    ## 7         A2 =~    PopAnh  1.009923222 0.0478850707918392
    ## 9         A2 =~  PopGuilt -0.122992000  0.670051708348395
    ## 15        A3 =~ PopSleDec  0.613564800  0.124372542681528
    ## 16        A3 =~ PopSleInc  0.533467248  0.121116456542721
    ## 14        A3 =~  PopFatig  0.699837589  0.125359057821666
    ## 12        A3 =~ PopAppDec  0.169771097  0.091681784428086
    ## 13        A3 =~ PopAppInc  0.424582116 0.0966693072453703
    ## 18    PopAnh ~~    PopAnh  0.001000109 0.0730501120652715
    ## 6         A1 ~~        A3  0.989999970  0.174045365394957
    ## 24  PopGuilt ~~  PopGuilt  0.514971157  0.252625866099056
    ## 21   PopConc ~~   PopConc  0.410016869   0.28307300636008
    ## 27    PopSui ~~    PopSui  0.605555444  0.268528202868254
    ## 22    PopDep ~~    PopDep  0.199459412 0.0761930410923394
    ## 25 PopSleDec ~~ PopSleDec  0.623533510  0.320166284105861
    ## 26 PopSleInc ~~ PopSleInc  0.715409809  0.278715355213914
    ## 23  PopFatig ~~  PopFatig  0.510220909   0.35495390178922
    ## 19 PopAppDec ~~ PopAppDec  0.971177490  0.252195457494435
    ## 20 PopAppInc ~~ PopAppInc  0.819730091  0.169314812101568
    ## 5         A1 ~~        A2  0.849379795  0.139673482127752
    ## 11        A2 ~~        A3  0.882966311  0.134680002665176
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
    ## 1         A1 =~    PopDep              0.20895926 0.023590234
    ## 2         A1 =~    PopAnh              0.25521142 0.022283575
    ## 3         A1 =~  PopGuilt              0.16157850 0.021242487
    ## 4         A1 =~   PopConc              0.18803340 0.026059182
    ## 5         A1 =~    PopSui              0.11327220 0.019550095
    ## 6         A2 =~ PopAppDec              0.02984344 0.016369248
    ## 7         A2 =~ PopAppInc              0.11206053 0.027097291
    ## 8         A2 =~ PopSleDec              0.12203179 0.026473594
    ## 9         A2 =~ PopSleInc              0.10904579 0.026079580
    ## 10        A2 =~  PopFatig              0.15615479 0.030569745
    ## 13    PopDep ~~    PopAnh              0.02495372 0.009882078
    ## 14    PopDep ~~    PopDep              0.03857167 0.011005221
    ## 15    PopAnh ~~    PopAnh              0.02221835 0.010755374
    ## 16  PopGuilt ~~  PopGuilt              0.03392092 0.009672605
    ## 17   PopConc ~~   PopConc              0.02376639 0.015330675
    ## 18    PopSui ~~    PopSui              0.02044728 0.008648280
    ## 19 PopAppDec ~~ PopAppDec              0.03513845 0.009093980
    ## 20 PopAppInc ~~ PopAppInc              0.06668565 0.013267158
    ## 21 PopSleDec ~~ PopSleDec              0.03049407 0.014490076
    ## 22 PopSleInc ~~ PopSleInc              0.03597900 0.013243305
    ## 23  PopFatig ~~  PopFatig              0.03339586 0.020361574
    ## 24        A1 ~~        A2              1.14047968 0.176097006

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
    ##  11.236 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 43.66147 33 0.1015084 87.66147 0.9962031 0.1159836

``` r
pop_psych_soma_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 3         A1 =~    PopDep    0.7463333 0.0848757053296728
    ## 1         A1 =~    PopAnh    0.8818497 0.0779738306816473
    ## 4         A1 =~  PopGuilt    0.6577839 0.0870802035046724
    ## 2         A1 =~   PopConc    0.7662171  0.107665696341588
    ## 5         A1 =~    PopSui    0.6211545  0.107407206134155
    ## 8         A2 =~ PopAppDec    0.1728245  0.093229556523004
    ## 9         A2 =~ PopAppInc    0.4342165 0.0970796039937196
    ## 11        A2 =~ PopSleDec    0.6253954  0.124565444260839
    ## 12        A2 =~ PopSleInc    0.5429833  0.121738721336936
    ## 10        A2 =~  PopFatig    0.7150330  0.124714248681582
    ## 18    PopDep ~~    PopAnh    0.2654810  0.122650576298681
    ## 7         A1 ~~        A2    1.0000000   0.13507570816249
    ## 19    PopDep ~~    PopDep    0.4429860   0.13950513006424
    ## 14    PopAnh ~~    PopAnh    0.2223412   0.12972885907025
    ## 21  PopGuilt ~~  PopGuilt    0.5673197   0.16114493136982
    ## 17   PopConc ~~   PopConc    0.4129122  0.259678568514862
    ## 24    PopSui ~~    PopSui    0.6141668  0.259996596921635
    ## 15 PopAppDec ~~ PopAppDec    0.9701315  0.252087721068806
    ## 16 PopAppInc ~~ PopAppInc    0.8114557  0.169982328717346
    ## 22 PopSleDec ~~ PopSleDec    0.6088737  0.320192803741142
    ## 23 PopSleInc ~~ PopSleInc    0.7051660  0.279396854294639
    ## 20  PopFatig ~~  PopFatig    0.4887172   0.35534533954665
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
    ##   1.798 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##     chisq df   p_chisq    AIC CFI       SRMR
    ## df 20.944 24 0.6420268 82.944   1 0.06476125

``` r
pop_psych_soma_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 14        A1 =~    PopDep   0.40088364  0.191139920839704
    ## 12        A1 =~    PopAnh   0.14557722  0.156450302595406
    ## 15        A1 =~  PopGuilt   0.29186126  0.181520097361726
    ## 13        A1 =~   PopConc  -0.29170743  0.258306449243955
    ## 16        A1 =~    PopSui   0.74761130  0.355748307490071
    ## 20        A2 =~ PopAppDec   0.14274871  0.193370770888637
    ## 21        A2 =~ PopAppInc   0.15297283  0.157377271667274
    ## 23        A2 =~ PopSleDec   0.81735273  0.432313294296805
    ## 24        A2 =~ PopSleInc  -0.56545376    0.2994146203932
    ## 22        A2 =~  PopFatig  -0.32374135  0.233796783166672
    ## 5          A =~    PopDep   0.66193772 0.0904925926760811
    ## 1          A =~    PopAnh   0.83206978  0.075713342728453
    ## 7          A =~  PopGuilt   0.60829813  0.101794318146327
    ## 4          A =~   PopConc   0.90913573  0.126620005995126
    ## 10         A =~    PopSui   0.48509283   0.14874628217132
    ## 2          A =~ PopAppDec   0.17751833  0.097326410743278
    ## 3          A =~ PopAppInc   0.45535310  0.089043122068905
    ## 8          A =~ PopSleDec   0.69830803   0.11803779405946
    ## 9          A =~ PopSleInc   0.59163294  0.117655038612346
    ## 6          A =~  PopFatig   0.76105375  0.110300568447621
    ## 31    PopDep ~~    PopAnh   0.31449714  0.120492714303575
    ## 32    PopDep ~~    PopDep   0.40113100  0.180758386614144
    ## 27    PopAnh ~~    PopAnh   0.28646723  0.116338278846867
    ## 34  PopGuilt ~~  PopGuilt   0.54479038  0.171855740811904
    ## 30   PopConc ~~   PopConc   0.08838019   0.33700103307499
    ## 37    PopSui ~~    PopSui   0.20576271  0.534655569076136
    ## 28 PopAppDec ~~ PopAppDec   0.94811034   0.25022854736613
    ## 29 PopAppInc ~~ PopAppInc   0.76925292  0.178709405460377
    ## 35 PopSleDec ~~ PopSleDec  -0.15569907  0.819627618074095
    ## 36 PopSleInc ~~ PopSleInc   0.33023272  0.437859249371648
    ## 33  PopFatig ~~  PopFatig   0.31598852  0.367043188096156
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
    ##   1.667 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 43.52263 33 0.1040511 87.52263 0.9962526 0.1132743

``` r
pop_psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 2         A1 =~    PopDep    0.7790812 0.0932775857891792
    ## 1         A1 =~    PopAnh    0.9160790 0.0899196803280113
    ## 3         A1 =~  PopGuilt    0.6511046 0.0875781316689064
    ## 4         A1 =~    PopSui    0.6192688  0.107666654812179
    ## 7         A2 =~ PopAppDec    0.1797426  0.097533205761256
    ## 8         A2 =~ PopAppInc    0.4556646 0.0923795660598046
    ## 11        A2 =~ PopSleDec    0.6556109  0.118053794037893
    ## 12        A2 =~ PopSleInc    0.5636417  0.118386646372817
    ## 10        A2 =~  PopFatig    0.7505886  0.115601361174178
    ## 9         A2 =~   PopConc    0.8111338  0.115091888866875
    ## 18    PopDep ~~    PopAnh    0.2099351  0.142699339061326
    ## 19    PopDep ~~    PopDep    0.3930325  0.155110160873163
    ## 14    PopAnh ~~    PopAnh    0.1607995   0.15472558678029
    ## 21  PopGuilt ~~  PopGuilt    0.5760630   0.16066303148368
    ## 24    PopSui ~~    PopSui    0.6165072  0.259341369025889
    ## 15 PopAppDec ~~ PopAppDec    0.9676934  0.251725588456459
    ## 16 PopAppInc ~~ PopAppInc    0.7923687  0.168625031295301
    ## 22 PopSleDec ~~ PopSleDec    0.5701741  0.321297793303327
    ## 23 PopSleInc ~~ PopSleInc    0.6823082  0.276627522377586
    ## 20  PopFatig ~~  PopFatig    0.4366168  0.342150249172555
    ## 17   PopConc ~~   PopConc    0.3420612  0.265682735349156
    ## 6         A1 ~~        A2    0.8988845  0.101840243117401
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
    ##   1.834 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 17.12126 23 0.8032588 81.12126   1 0.07088151

``` r
pop_psych_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 13        A1 =~    PopDep    0.4627994  0.207253150230514
    ## 12        A1 =~    PopAnh    0.1880397   0.16191562425803
    ## 14        A1 =~  PopGuilt    0.3567082  0.170603199445206
    ## 15        A1 =~    PopSui    0.7699264  0.367395165570965
    ## 19        A2 =~ PopAppDec   -0.1405390  0.197486228291463
    ## 20        A2 =~ PopAppInc    0.5349435  0.306338594673598
    ## 23        A2 =~ PopSleDec    0.2528736  0.251917149427386
    ## 24        A2 =~ PopSleInc   -0.2883060  0.252694944308623
    ## 22        A2 =~  PopFatig   -0.3628558  0.320417342182426
    ## 21        A2 =~   PopConc    0.3365616  0.251176346976691
    ## 5          A =~    PopDep    0.6211088  0.092889201335731
    ## 1          A =~    PopAnh    0.8267491 0.0857122272272616
    ## 7          A =~  PopGuilt    0.5771857  0.102895884270962
    ## 10         A =~    PopSui    0.4279589  0.145603585385493
    ## 2          A =~ PopAppDec    0.1891572  0.100815547361925
    ## 3          A =~ PopAppInc    0.4665192 0.0972241942196086
    ## 8          A =~ PopSleDec    0.6917286  0.120870594645275
    ## 9          A =~ PopSleInc    0.6396647  0.126941343430577
    ## 6          A =~  PopFatig    0.8055356  0.123075166962883
    ## 4          A =~   PopConc    0.8277940  0.117754374753204
    ## 31    PopDep ~~    PopAnh    0.3231089  0.118314545282717
    ## 36 PopSleDec ~~ PopSleInc   -0.4068722   0.24786648429591
    ## 32    PopDep ~~    PopDep    0.4000404  0.192403564542114
    ## 27    PopAnh ~~    PopAnh    0.2811268  0.121832791141398
    ## 34  PopGuilt ~~  PopGuilt    0.5396161  0.176317272204204
    ## 38    PopSui ~~    PopSui    0.2240646  0.549825889239673
    ## 28 PopAppDec ~~ PopAppDec    0.9444682  0.255850608840996
    ## 29 PopAppInc ~~ PopAppInc    0.4961954  0.358073404874782
    ## 35 PopSleDec ~~ PopSleDec    0.4575666  0.349083998983177
    ## 37 PopSleInc ~~ PopSleInc    0.5077083  0.331202448297237
    ## 33  PopFatig ~~  PopFatig    0.2194475  0.432888450597109
    ## 30   PopConc ~~   PopConc    0.2014827  0.311102044454361
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
    ##   1.745 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI      SRMR
    ## df 36.04196 33 0.3281296 80.04196 0.9989167 0.1101096

``` r
pop_affect_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.8007195  0.088870295991572
    ## 2         A1 =~  PopGuilt    0.7178417   0.09487326720562
    ## 3         A1 =~    PopSui    0.6855533  0.117324932542369
    ## 6         A2 =~    PopAnh    0.8901032 0.0762523031891244
    ## 8         A2 =~ PopAppInc    0.4517705 0.0873930320888673
    ## 7         A2 =~ PopAppDec    0.1765200 0.0960674476478214
    ## 12        A2 =~ PopSleInc    0.5605483   0.11481260826087
    ## 11        A2 =~ PopSleDec    0.6432638  0.113427800273179
    ## 10        A2 =~  PopFatig    0.7423364  0.107265176933694
    ## 9         A2 =~   PopConc    0.7986420  0.109564531819925
    ## 18    PopDep ~~    PopAnh    0.3388393  0.112313167701726
    ## 19    PopDep ~~    PopDep    0.3588489  0.151814698101874
    ## 21  PopGuilt ~~  PopGuilt    0.4847044  0.176752687528473
    ## 24    PopSui ~~    PopSui    0.5300172  0.269713449989323
    ## 14    PopAnh ~~    PopAnh    0.2077162  0.128273334772596
    ## 16 PopAppInc ~~ PopAppInc    0.7959035  0.165336854876618
    ## 15 PopAppDec ~~ PopAppDec    0.9688407  0.251641355705667
    ## 23 PopSleInc ~~ PopSleInc    0.6857854  0.274936612124243
    ## 22 PopSleDec ~~ PopSleDec    0.5862115  0.317453950094168
    ## 20  PopFatig ~~  PopFatig    0.4489373  0.337516892973848
    ## 17   PopConc ~~   PopConc    0.3621699  0.261227964707112
    ## 5         A1 ~~        A2    0.8205101 0.0623316128124066
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
    ##   1.778 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ##   17.12 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00227387839846352 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.201559380909594 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 15.26185 24 0.9128626 77.26185   1 0.06597186

``` r
pop_affect_veg_bif_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs  STD_Genotype    STD_Genotype_SE
    ## 12        A1 =~    PopDep  0.4341007279  0.218254488730081
    ## 13        A1 =~  PopGuilt  0.3133817134  0.158905017845656
    ## 14        A1 =~    PopSui  0.6764793245  0.333938348628582
    ## 18        A2 =~    PopAnh  0.1953571384  0.101797743767918
    ## 20        A2 =~ PopAppInc -0.1377905333  0.161774304568031
    ## 19        A2 =~ PopAppDec -0.1160833150  0.182949426091087
    ## 24        A2 =~ PopSleInc  0.6059419816  0.267072456590023
    ## 23        A2 =~ PopSleDec -0.7246674628   0.34714313781337
    ## 22        A2 =~  PopFatig  0.3335137812  0.228337281190834
    ## 21        A2 =~   PopConc -0.1588770912  0.193881255234378
    ## 5          A =~    PopDep  0.6413492489 0.0861325726008565
    ## 7          A =~  PopGuilt  0.6090200594 0.0914806487780227
    ## 10         A =~    PopSui  0.5084470762  0.111541136597094
    ## 1          A =~    PopAnh  0.8822613436 0.0812330761835293
    ## 3          A =~ PopAppInc  0.4646515913 0.0908980029762933
    ## 2          A =~ PopAppDec  0.1864985313 0.0970878095743529
    ## 9          A =~ PopSleInc  0.5480145223  0.129344691152798
    ## 8          A =~ PopSleDec  0.7371693224  0.132657360620863
    ## 6          A =~  PopFatig  0.7256972532  0.111790162982498
    ## 4          A =~   PopConc  0.8203654989  0.113061537595685
    ## 31    PopDep ~~    PopAnh  0.3577971821  0.114273631906121
    ## 35 PopSleDec ~~ PopSleDec  0.0009999992   0.66147076196586
    ## 36 PopSleInc ~~ PopSleInc  0.3325132172  0.404679153224083
    ## 32    PopDep ~~    PopDep  0.4002274892  0.198992216903642
    ## 34  PopGuilt ~~  PopGuilt  0.5308864698  0.177040961096217
    ## 37    PopSui ~~    PopSui  0.2838574754  0.466991339845517
    ## 27    PopAnh ~~    PopAnh  0.1834506232  0.143195094515664
    ## 29 PopAppInc ~~ PopAppInc  0.7651120400  0.180193006262125
    ## 28 PopAppDec ~~ PopAppDec  0.9517428144  0.249066398760246
    ## 33  PopFatig ~~  PopFatig  0.3621324414  0.353188028749244
    ## 30   PopConc ~~   PopConc  0.3017583606  0.267768720447568
    ## 16        A1 ~~        A1  1.0000000000                   
    ## 26        A2 ~~        A2  1.0000000000                   
    ## 11         A ~~         A  1.0000000000                   
    ## 15        A1 ~~         A  0.0000000000                   
    ## 25        A2 ~~         A  0.0000000000                   
    ## 17        A1 ~~        A2  0.0000000000

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

| Model  | Name                  | chisq |  df | p_chisq |  AIC |   CFI |   SRMR |  dAIC |
|:-------|:----------------------|------:|----:|--------:|-----:|------:|-------:|------:|
| 1a     | Common                |  59.7 |  35 | 0.00574 | 99.7 | 0.991 | 0.1220 | 22.40 |
| 1b     | Common (gating)       |  43.7 |  34 | 0.12400 | 85.7 | 0.997 | 0.1160 |  8.40 |
| 1c     | Common (App)          |  42.7 |  33 | 0.11900 | 86.7 | 0.997 | 0.1140 |  9.47 |
| 1d     | Common (Sle)          |  39.0 |  33 | 0.21700 | 83.0 | 0.998 | 0.1030 |  5.77 |
| 1er    | Common (App,Sle)      |  37.5 |  32 | 0.23200 | 83.5 | 0.998 | 0.1010 |  6.23 |
| 2a(i)  | Psych-Somatic         |  43.7 |  33 | 0.10200 | 87.7 | 0.996 | 0.1160 | 10.40 |
| 2a(ii) | Psych-Somatic (BiF)   |  20.9 |  24 | 0.64200 | 82.9 | 1.000 | 0.0648 |  5.68 |
| 2b(i)  | Psych-Neuroveg        |  43.5 |  33 | 0.10400 | 87.5 | 0.996 | 0.1130 | 10.30 |
| 2b(ii) | Psych-Neuroveg (BiF)  |  17.1 |  23 | 0.80300 | 81.1 | 1.000 | 0.0709 |  3.86 |
| 2c(i)  | Affect-Neuroveg       |  36.0 |  33 | 0.32800 | 80.0 | 0.999 | 0.1100 |  2.78 |
| 2c(ii) | Affect-Neuroveg (BiF) |  15.3 |  24 | 0.91300 | 77.3 | 1.000 | 0.0660 |  0.00 |
| 3      | Cog-Mood-Neuroveg     |  50.3 |  31 | 0.01580 | 98.3 | 0.993 | 0.1160 | 21.00 |

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
    ##   3.217 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0464961771558631 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.65818091523924 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df      p_chisq      AIC      CFI     SRMR
    ## df 212.2718 94 3.948104e-11 296.2718 0.974773 0.151407

``` r
clin_pop_affect_veg.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 12    ClinSOMA =~  ClinAppInc   0.46478678  0.213545859937147 2.951540e-02
    ## 14    ClinSOMA =~  ClinSleDec   0.52526183  0.259341130323294 4.282859e-02
    ## 15    ClinSOMA =~  ClinSleInc   0.51696908   0.27111150393746 5.653908e-02
    ## 13    ClinSOMA =~ ClinMotoInc   0.46468518  0.206011522921712 2.409356e-02
    ## 26   PopAFFECT =~      PopDep   0.73098938 0.0752779601712589 2.719247e-22
    ## 27   PopAFFECT =~    PopGuilt   0.76678384 0.0938705690247469 3.122160e-16
    ## 28   PopAFFECT =~      PopSui   0.69122586  0.110017318701051 3.323705e-10
    ## 39 PopNEUROVEG =~      PopAnh   0.80640013 0.0647025883835088 1.185424e-35
    ## 41 PopNEUROVEG =~   PopAppInc   0.45594204 0.0834300222520314 4.630037e-08
    ## 40 PopNEUROVEG =~   PopAppDec   0.18520095 0.0903702428059009 4.042802e-02
    ## 45 PopNEUROVEG =~   PopSleInc   0.49637185  0.107687791560124 4.038969e-06
    ## 44 PopNEUROVEG =~   PopSleDec   0.64157305  0.106087612353727 1.470392e-09
    ## 43 PopNEUROVEG =~    PopFatig   0.75132243  0.100574945638545 8.003035e-14
    ## 42 PopNEUROVEG =~     PopConc   0.79943197  0.103456661641944 1.099200e-14
    ## 35      PopDep ~~      PopAnh   0.42911614   0.08546601159953 5.143050e-07
    ## 7   ClinAppInc ~~  ClinAppDec  -0.71095943  0.274513710838718 9.600926e-03
    ## 8   ClinAppInc ~~  ClinAppInc   0.78397548  0.337469569070914 2.017405e-02
    ## 10  ClinSleDec ~~  ClinSleDec   0.72410111  0.641424534657194 2.589481e-01
    ## 11  ClinSleInc ~~  ClinSleInc   0.73274594  0.724009599774088 3.115075e-01
    ## 9  ClinMotoInc ~~ ClinMotoInc   0.78407039  0.465851077119523 9.235725e-02
    ## 36      PopDep ~~      PopDep   0.46565486  0.121772514041663 1.313199e-04
    ## 38    PopGuilt ~~    PopGuilt   0.41204248   0.18245996735246 2.392855e-02
    ## 49      PopSui ~~      PopSui   0.52220799  0.250658273171886 3.722034e-02
    ## 31      PopAnh ~~      PopAnh   0.34971922  0.098573131372836 3.884779e-04
    ## 33   PopAppInc ~~   PopAppInc   0.79211699  0.152683252124017 2.125847e-07
    ## 32   PopAppDec ~~   PopAppDec   0.96570016  0.232613466279319 3.302605e-05
    ## 48   PopSleInc ~~   PopSleInc   0.75361478  0.246076736456086 2.194791e-03
    ## 47   PopSleDec ~~   PopSleDec   0.58838457   0.28006751181586 3.565290e-02
    ## 37    PopFatig ~~    PopFatig   0.43551541  0.309415260576707 1.592647e-01
    ## 34     PopConc ~~     PopConc   0.36090884  0.243108199821397 1.376618e-01
    ## 3   ClinAPPDEC ~~  ClinAPPDEC   1.00000000  0.202971570106521 8.358793e-07
    ## 23     ClinSUI ~~     ClinSUI   1.00000000  0.287477153866971 5.041583e-04
    ## 16    ClinSOMA ~~  ClinAPPDEC   0.61352787  0.350317817596315 7.988567e-02
    ## 18    ClinSOMA ~~     ClinSUI   0.16653596  0.259032867927785 5.202815e-01
    ## 19    ClinSOMA ~~   PopAFFECT   0.09080505  0.156552205083106 5.618941e-01
    ## 20    ClinSOMA ~~ PopNEUROVEG   0.55119218  0.213934441977209 9.981611e-03
    ## 4   ClinAPPDEC ~~     ClinSUI  -0.07697190  0.163529182953663 6.378749e-01
    ## 5   ClinAPPDEC ~~   PopAFFECT  -0.11933447 0.0919468204812059 1.943323e-01
    ## 6   ClinAPPDEC ~~ PopNEUROVEG  -0.09818471 0.0954491834414475 3.036394e-01
    ## 24     ClinSUI ~~   PopAFFECT   0.68611334   0.12135533510905 1.569908e-08
    ## 25     ClinSUI ~~ PopNEUROVEG   0.56226991  0.115202593970764 1.057086e-06
    ## 30   PopAFFECT ~~ PopNEUROVEG   0.80891869 0.0634340207836651 3.033603e-37
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
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinMotoInc ClinMotoDec 
    ##       0.313       0.005       0.005       0.580       0.971       0.885 
    ##     ClinSui 
    ##       0.008 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec   0.470  -0.679         
    ## ClinAppInc   0.329   0.932   0.132 
    ## ClinSleDec   0.982   0.177         
    ## ClinSleInc   0.383           0.516 
    ## ClinMotoInc          0.144         
    ## ClinMotoDec  0.238           0.236 
    ## ClinSui     -0.153   0.150   0.973 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      1.522   1.414   1.296
    ## Proportion Var   0.217   0.202   0.185
    ## Cumulative Var   0.217   0.420   0.605
    ## 
    ## The degrees of freedom for the model is 3 and the fit was 0.4196

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
    ##      0.121      0.005      0.460      0.649      0.005      0.668      0.398 
    ## PopMotoDec   PopFatig   PopGuilt    PopConc     PopSui 
    ##      0.005      0.641      0.578      0.447      0.654 
    ## 
    ## Loadings:
    ##            Factor1 Factor2 Factor3
    ## PopDep      0.874   0.146   0.307 
    ## PopAnh      0.936           0.339 
    ## PopAppDec   0.105   0.689   0.234 
    ## PopAppInc   0.318  -0.467   0.180 
    ## PopSleDec                   0.993 
    ## PopSleInc   0.571                 
    ## PopMotoInc          0.760   0.142 
    ## PopMotoDec  0.255   0.942  -0.207 
    ## PopFatig    0.510  -0.248   0.194 
    ## PopGuilt    0.431   0.285   0.393 
    ## PopConc     0.374           0.642 
    ## PopSui      0.482   0.279   0.189 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      2.973   2.405   1.991
    ## Proportion Var   0.248   0.200   0.166
    ## Cumulative Var   0.248   0.448   0.614
    ## 
    ## The degrees of freedom for the model is 33 and the fit was 5.0691

# All symptoms

``` r
all_covstruct_prefix <- 'all.covstruct'
all_covstruct_r <- file.path('ldsc', paste(all_covstruct_prefix, 'deparse.R', sep='.'))
all_covstruct_rds <- file.path('ldsc', paste(all_covstruct_prefix, 'rds', sep='.'))

all_symptoms_covstruct <- dget(all_covstruct_r)

all_sumstats_prevs <- read_tsv(file.path('ldsc', paste(all_covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 12 Columns: 5
    ## ── Column specification ──────────────────────────────────────────────────
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
    ##   1.824 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 5       A1 =~     Dep   0.83968329 0.0477003658399802
    ## 1       A1 =~     Anh   0.98187134 0.0445271644844958
    ## 2       A1 =~  AppDec   0.05529376  0.081299269629888
    ## 3       A1 =~  AppInc   0.42227844 0.0713889717331678
    ## 9       A1 =~  SleDec   0.38758192 0.0958970348683241
    ## 10      A1 =~  SleInc   0.51023734 0.0905215116066048
    ## 8       A1 =~ MotoInc   0.20096623  0.121063257504148
    ## 6       A1 =~   Fatig   0.63615156 0.0986747020064327
    ## 7       A1 =~   Guilt   0.59632487 0.0730312292328309
    ## 4       A1 =~    Conc   0.73928338  0.102921749935823
    ## 11      A1 =~     Sui   0.64319048 0.0839617489415137
    ## 17     Dep ~~     Dep   0.29493183 0.0719853203530205
    ## 13     Anh ~~     Anh   0.03592706 0.0632420890196112
    ## 14  AppDec ~~  AppDec   0.99694481  0.197985174074286
    ## 15  AppInc ~~  AppInc   0.82168119  0.134232990857195
    ## 21  SleDec ~~  SleDec   0.84977815  0.222625738559737
    ## 22  SleInc ~~  SleInc   0.73965396  0.233763008396691
    ## 20 MotoInc ~~ MotoInc   0.95961349  0.472026748314184
    ## 18   Fatig ~~   Fatig   0.59531572  0.331399688871215
    ## 19   Guilt ~~   Guilt   0.64439664   0.14623897783607
    ## 16    Conc ~~    Conc   0.45345876  0.276572757575262
    ## 23     Sui ~~     Sui   0.58630609  0.177899446800531
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
    ##   1.572 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 5       A1 =~     Dep   0.70054970 0.0739031513035017
    ## 1       A1 =~     Anh   0.84894801 0.0687373545974629
    ## 2       A1 =~  AppDec   0.04385646 0.0869815492385903
    ## 3       A1 =~  AppInc   0.46776516 0.0784550500054788
    ## 9       A1 =~  SleDec   0.41438865  0.103308588582375
    ## 10      A1 =~  SleInc   0.54586163 0.0984076339452935
    ## 8       A1 =~ MotoInc   0.21775270  0.130616766234914
    ## 6       A1 =~   Fatig   0.68768252  0.104526923354352
    ## 7       A1 =~   Guilt   0.65211267 0.0804460345008033
    ## 4       A1 =~    Conc   0.81140087    0.1127322650755
    ## 11      A1 =~     Sui   0.68603568 0.0909916621050091
    ## 17     Dep ~~     Anh   0.30932139  0.101438912664712
    ## 18     Dep ~~     Dep   0.50923055  0.114099114567104
    ## 13     Anh ~~     Anh   0.27928535  0.111102486357975
    ## 14  AppDec ~~  AppDec   0.99807408  0.197864459642426
    ## 15  AppInc ~~  AppInc   0.78119557   0.13640160559586
    ## 22  SleDec ~~  SleDec   0.82828188  0.222798891210219
    ## 23  SleInc ~~  SleInc   0.70203555  0.231379394135635
    ## 21 MotoInc ~~ MotoInc   0.95257859  0.470855918751975
    ## 19   Fatig ~~   Fatig   0.52709381  0.327922505557557
    ## 20   Guilt ~~   Guilt   0.57474862  0.147443191001018
    ## 16    Conc ~~    Conc   0.34162961  0.273384928772342
    ## 24     Sui ~~     Sui   0.52935390  0.177912531889805
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
    ##   1.723 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 5       A1 =~     Dep    0.7012332 0.0736723891418023
    ## 1       A1 =~     Anh    0.8482954 0.0683461487672083
    ## 2       A1 =~  AppDec    0.1050515 0.0890117271188368
    ## 3       A1 =~  AppInc    0.4723937  0.078344294673896
    ## 9       A1 =~  SleDec    0.4189116  0.103050987366891
    ## 10      A1 =~  SleInc    0.5444870 0.0981726955795622
    ## 8       A1 =~ MotoInc    0.2219199  0.130497346433741
    ## 6       A1 =~   Fatig    0.6872792  0.104826138395791
    ## 7       A1 =~   Guilt    0.6482263 0.0802746992883348
    ## 4       A1 =~    Conc    0.8096262  0.113056133817283
    ## 11      A1 =~     Sui    0.6853122 0.0913054733433305
    ## 18     Dep ~~     Anh    0.3092002  0.100862740111008
    ## 15  AppDec ~~  AppInc   -0.4752880  0.111375435177482
    ## 19     Dep ~~     Dep    0.5082736  0.113721567488788
    ## 13     Anh ~~     Anh    0.2803963  0.110429152074917
    ## 14  AppDec ~~  AppDec    0.9889648  0.197643880683124
    ## 16  AppInc ~~  AppInc    0.7768443  0.136286967365669
    ## 23  SleDec ~~  SleDec    0.8245135  0.222525258100705
    ## 24  SleInc ~~  SleInc    0.7035350    0.2314952106968
    ## 22 MotoInc ~~ MotoInc    0.9507551  0.470511829323942
    ## 20   Fatig ~~   Fatig    0.5276478  0.327827229372466
    ## 21   Guilt ~~   Guilt    0.5798037  0.147484830701272
    ## 17    Conc ~~    Conc    0.3445058  0.272753280180982
    ## 25     Sui ~~     Sui    0.5303480  0.177916023930225
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
    ##   1.637 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 5       A1 =~     Dep    0.6942259 0.0727788226415146
    ## 1       A1 =~     Anh    0.8407633 0.0675842313310787
    ## 2       A1 =~  AppDec    0.0447404 0.0870483309511256
    ## 3       A1 =~  AppInc    0.4690943 0.0782922584830766
    ## 9       A1 =~  SleDec    0.4462266  0.103658858747958
    ## 10      A1 =~  SleInc    0.5696780 0.0990868246447021
    ## 8       A1 =~ MotoInc    0.2198810  0.130679570613235
    ## 6       A1 =~   Fatig    0.6884079  0.104301908278302
    ## 7       A1 =~   Guilt    0.6520736 0.0801165221530283
    ## 4       A1 =~    Conc    0.8140762  0.112497624114562
    ## 11      A1 =~     Sui    0.6845436 0.0905582889232199
    ## 17     Dep ~~     Anh    0.3203724 0.0989076332005006
    ## 23  SleDec ~~  SleInc   -0.2963452  0.152708790112864
    ## 18     Dep ~~     Dep    0.5180496  0.112247113450044
    ## 13     Anh ~~     Anh    0.2931168  0.107925697077128
    ## 14  AppDec ~~  AppDec    0.9979990  0.197853870560069
    ## 15  AppInc ~~  AppInc    0.7799499   0.13631963574815
    ## 22  SleDec ~~  SleDec    0.8008817  0.220565402212937
    ## 24  SleInc ~~  SleInc    0.6754675   0.22948842116874
    ## 21 MotoInc ~~ MotoInc    0.9516522  0.470795610156855
    ## 19   Fatig ~~   Fatig    0.5260950  0.327813046358206
    ## 20   Guilt ~~   Guilt    0.5748023  0.147171083038489
    ## 16    Conc ~~    Conc    0.3372808  0.273197364944619
    ## 25     Sui ~~     Sui    0.5313993  0.177850145132922
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
    ##    1.64 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## 5       A1 =~     Dep    0.6947461 0.0725249739421951
    ## 1       A1 =~     Anh    0.8399426 0.0671730831763142
    ## 2       A1 =~  AppDec    0.1063360 0.0891340788281837
    ## 3       A1 =~  AppInc    0.4738463 0.0781782266681098
    ## 9       A1 =~  SleDec    0.4509577  0.103421684393473
    ## 10      A1 =~  SleInc    0.5686700  0.098823679049497
    ## 8       A1 =~ MotoInc    0.2240895  0.130562639930788
    ## 6       A1 =~   Fatig    0.6880225  0.104599307523278
    ## 7       A1 =~   Guilt    0.6481385 0.0799418208140905
    ## 4       A1 =~    Conc    0.8123050  0.112814097850724
    ## 11      A1 =~     Sui    0.6837762 0.0908636010527709
    ## 18     Dep ~~     Anh    0.3205051 0.0982887358749813
    ## 15  AppDec ~~  AppInc   -0.4760506   0.11146256749066
    ## 24  SleDec ~~  SleInc   -0.2985863  0.152733082301356
    ## 19     Dep ~~     Dep    0.5173261  0.111834431077799
    ## 13     Anh ~~     Anh    0.2944975   0.10721300172793
    ## 14  AppDec ~~  AppDec    0.9886917  0.197600892318499
    ## 16  AppInc ~~  AppInc    0.7754696   0.13620063213386
    ## 23  SleDec ~~  SleDec    0.7966420  0.220261765091348
    ## 25  SleInc ~~  SleInc    0.6766122  0.229575780640618
    ## 22 MotoInc ~~ MotoInc    0.9497810    0.4704493865686
    ## 20   Fatig ~~   Fatig    0.5266259  0.327727679521777
    ## 21   Guilt ~~   Guilt    0.5799160  0.147211031219142
    ## 17    Conc ~~    Conc    0.3401598  0.272566685412859
    ## 26     Sui ~~     Sui    0.5324491  0.177846995948142
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
    ##   1.931 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 96.48374 41 2.282677e-06 146.4837 0.9842634 0.1378455

``` r
all_psych_soma.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 3       A1 =~     Dep    0.7227216 0.0770185534566662
    ## 1       A1 =~     Anh    0.8698566 0.0714893884134782
    ## 4       A1 =~   Guilt    0.6452514 0.0811448471116944
    ## 2       A1 =~    Conc    0.8013119  0.113632777906561
    ## 5       A1 =~     Sui    0.6879581 0.0923473638758412
    ## 8       A2 =~  AppDec   -0.1244413 0.0997268288836779
    ## 9       A2 =~  AppInc   -0.5254150  0.089444602378293
    ## 12      A2 =~  SleDec   -0.4539713  0.115462393600389
    ## 13      A2 =~  SleInc   -0.5977304  0.110480070913635
    ## 11      A2 =~ MotoInc   -0.2628048  0.144050051772288
    ## 10      A2 =~   Fatig   -0.7699839  0.120408704148885
    ## 20     Dep ~~     Anh    0.2753883  0.106715807812682
    ## 17  AppDec ~~  AppInc   -0.4910474  0.114088300556855
    ## 21     Dep ~~     Dep    0.4776733  0.119140574747902
    ## 15     Anh ~~     Anh    0.2433491  0.117552252472401
    ## 23   Guilt ~~   Guilt    0.5836491   0.14788305714655
    ## 19    Conc ~~    Conc    0.3578969  0.272528206397198
    ## 27     Sui ~~     Sui    0.5267131  0.177938180636123
    ## 16  AppDec ~~  AppDec    0.9845135   0.19768620228175
    ## 18  AppInc ~~  AppInc    0.7239384   0.14201151229232
    ## 25  SleDec ~~  SleDec    0.7939115   0.22829040682947
    ## 26  SleInc ~~  SleInc    0.6427162   0.23300486170914
    ## 24 MotoInc ~~ MotoInc    0.9309434  0.469858432012195
    ## 22   Fatig ~~   Fatig    0.4071216  0.344573899058617
    ## 7       A1 ~~      A2   -0.8480733  0.095991030080792
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
    ##   1.905 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##      chisq df      p_chisq      AIC       CFI      SRMR
    ## df 95.8247 41 2.788318e-06 145.8247 0.9844503 0.1361691

``` r
all_psych_veg.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 2       A1 =~     Dep    0.7563642  0.083760892202858
    ## 1       A1 =~     Anh    0.9036945 0.0811546870752059
    ## 3       A1 =~   Guilt    0.6384848 0.0821096962224052
    ## 4       A1 =~     Sui    0.6860052  0.093727601355373
    ## 7       A2 =~  AppDec    0.1208109 0.0995535246136289
    ## 8       A2 =~  AppInc    0.5257788 0.0872225069312676
    ## 12      A2 =~  SleDec    0.4564607  0.113624207029561
    ## 13      A2 =~  SleInc    0.5930646  0.107663712600479
    ## 11      A2 =~ MotoInc    0.2619985  0.142276924052185
    ## 10      A2 =~   Fatig    0.7573010  0.115857539097833
    ## 9       A2 =~    Conc    0.9004399  0.125735807800422
    ## 20     Dep ~~     Anh    0.2205304  0.123939326272302
    ## 17  AppDec ~~  AppInc   -0.4891836  0.114140718817487
    ## 21     Dep ~~     Dep    0.4279141  0.131897347379842
    ## 15     Anh ~~     Anh    0.1833370  0.138502714952848
    ## 23   Guilt ~~   Guilt    0.5923363  0.147724064380569
    ## 27     Sui ~~     Sui    0.5293964  0.178124259393888
    ## 16  AppDec ~~  AppDec    0.9854068  0.197413808685375
    ## 18  AppInc ~~  AppInc    0.7235563  0.142884775371388
    ## 25  SleDec ~~  SleDec    0.7916447  0.225671063737039
    ## 26  SleInc ~~  SleInc    0.6482707  0.234805423352049
    ## 24 MotoInc ~~ MotoInc    0.9313557  0.470067946544241
    ## 22   Fatig ~~   Fatig    0.4264971  0.330303344351861
    ## 19    Conc ~~    Conc    0.1892088   0.29067167464079
    ## 6       A1 ~~      A2    0.8079316  0.088386135580479
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
    ##   1.942 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 79.01883 41 0.0003325217 129.0188 0.9892169 0.133106

``` r
all_affect_veg.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 1       A1 =~     Dep    0.7507985 0.0761037224593543
    ## 2       A1 =~   Guilt    0.7308977 0.0870950857628243
    ## 3       A1 =~     Sui    0.7913855  0.102736619277171
    ## 6       A2 =~     Anh    0.8635945 0.0686471311541394
    ## 7       A2 =~  AppDec    0.1116003 0.0939315422013457
    ## 8       A2 =~  AppInc    0.5050069 0.0813264280439218
    ## 12      A2 =~  SleDec    0.4312126  0.106971433089317
    ## 13      A2 =~  SleInc    0.5758865  0.101183554941174
    ## 11      A2 =~ MotoInc    0.2427216  0.135493669785824
    ## 10      A2 =~   Fatig    0.7299895  0.107283529566813
    ## 9       A2 =~    Conc    0.8589718  0.116630052553003
    ## 20     Dep ~~     Anh    0.4179745 0.0901040557119301
    ## 17  AppDec ~~  AppInc   -0.4820232  0.112740383257897
    ## 21     Dep ~~     Dep    0.4363016  0.123042867508139
    ## 23   Guilt ~~   Guilt    0.4657892  0.157142278994988
    ## 27     Sui ~~     Sui    0.3737065  0.188314837667349
    ## 15     Anh ~~     Anh    0.2542039  0.113380737844086
    ## 16  AppDec ~~  AppDec    0.9875454  0.197486626208001
    ## 18  AppInc ~~  AppInc    0.7449682  0.138280568027754
    ## 25  SleDec ~~  SleDec    0.8140570  0.222533669851152
    ## 26  SleInc ~~  SleInc    0.6683521  0.231856224992316
    ## 24 MotoInc ~~ MotoInc    0.9410846  0.470607398435154
    ## 22   Fatig ~~   Fatig    0.4671156  0.323254983910489
    ## 19    Conc ~~    Conc    0.2621648  0.274599196778721
    ## 5       A1 ~~      A2    0.7496738 0.0587062232331952
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
    ##  24.711 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0108199700205908 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.969657102936796 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 103.2905 39 1.011338e-07 157.2905 0.9817655 0.1381106

``` r
all_cog_mood_neuroveg.fit$results[c(1, 2, 3, 6, 7)]
```

    ##        lhs op     rhs STD_Genotype    STD_Genotype_SE
    ## 2       A1 =~   Guilt  1.110435103  0.826903696664544
    ## 1       A1 =~    Conc  0.777576799  0.123781327940573
    ## 3       A1 =~     Sui  0.682770566  0.100082048953792
    ## 8       A2 =~     Dep  0.876607859 0.0472773177404276
    ## 7       A2 =~     Anh  1.010535600  0.047417868037974
    ## 9       A2 =~   Guilt -0.443060942  0.837851327188041
    ## 12      A3 =~  AppDec -0.124623977 0.0994000728402758
    ## 13      A3 =~  AppInc -0.524566386 0.0892746518902304
    ## 16      A3 =~  SleDec -0.454275104  0.115857768744182
    ## 17      A3 =~  SleInc -0.597814275  0.110343588712852
    ## 15      A3 =~ MotoInc -0.263119972  0.144032752342656
    ## 14      A3 =~   Fatig -0.770876145  0.120941364378673
    ## 21  AppDec ~~  AppInc -0.491037607  0.113983682032011
    ## 19     Anh ~~     Anh  0.001000228 0.0669219346245701
    ## 5       A1 ~~      A2  0.872613264  0.114484010431727
    ## 26   Guilt ~~   Guilt  0.429268277  0.287869466024362
    ## 23    Conc ~~    Conc  0.395371623  0.273688453659517
    ## 30     Sui ~~     Sui  0.533821557  0.184325142445266
    ## 24     Dep ~~     Dep  0.231557975 0.0705862110459131
    ## 20  AppDec ~~  AppDec  0.984468748   0.19774379836734
    ## 22  AppInc ~~  AppInc  0.724830164  0.142048204208288
    ## 28  SleDec ~~  SleDec  0.793632627  0.228751448660134
    ## 29  SleInc ~~  SleInc  0.642619094  0.233229776677567
    ## 27 MotoInc ~~ MotoInc  0.930767041  0.469699262902349
    ## 25   Fatig ~~   Fatig  0.405751197  0.345677997621066
    ## 6       A1 ~~      A3 -0.801938470  0.103673928375522
    ## 11      A2 ~~      A3 -0.724264110 0.0982926598520353
    ## 4       A1 ~~      A1  1.000000000                   
    ## 10      A2 ~~      A2  1.000000000                   
    ## 18      A3 ~~      A3  1.000000000

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
