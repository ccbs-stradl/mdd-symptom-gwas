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
    ## platform       aarch64-apple-darwin20      
    ## arch           aarch64                     
    ## os             darwin20                    
    ## system         aarch64, darwin20           
    ## status                                     
    ## major          4                           
    ## minor          1.0                         
    ## year           2021                        
    ## month          05                          
    ## day            18                          
    ## svn rev        80317                       
    ## language       R                           
    ## version.string R version 4.1.0 (2021-05-18)
    ## nickname       Camp Pontanezen

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

    ## [1] '0.0.3'

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
MDD5;Motor⇅;Motor⇆;Psyc
MDD5a;Motor⇈;Motor⇉;PsycInc
MDD5b;Motor⇊;Motor⇇;PsycDec
MDD6;Fatigue;Fatigue;Fatig
MDD7;Guilt;Guilt;Guilt
MDD8;Concentrate;Concentrate;Conc
MDD9;Suicidality;Suicidality;Sui
", col_names=c('ref', 'h', 'v', 'abbv'), delim=';')
```

    ## Rows: 15 Columns: 4

    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────
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
MDD5b;Psychomotor retardation nearly every day
MDD6;Fatigue or loss of energy nearly every day
MDD7;Feelings of worthlessness or excessive or inappropriate guilt
MDD8;Diminished ability to think or concentrate, or indecisiveness
MDD9;Recurrent thoughts of death or suicide or a suicide attempt or a specific plan for attempting suicide
", col_names=c('Reference', 'Description'), delim=';')
```

    ## Rows: 15 Columns: 2

    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────
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

    ## # A tibble: 15 x 4
    ##    Reference Abbreviation Label     Description                                 
    ##    <chr>     <chr>        <chr>     <chr>                                       
    ##  1 MDD1      Dep          Mood      Depressed mood most of the day, nearly ever…
    ##  2 MDD2      Anh          Interest  Markedly diminished interest or pleasure in…
    ##  3 MDD3      App          Weight⇅   Significant change in weight or appetite    
    ##  4 MDD3a     AppDec       Weight⇊   Significant weight loss or decrease in appe…
    ##  5 MDD3b     AppInc       Weight⇈   Significant weight gain or increase in appe…
    ##  6 MDD4      Sle          Sleep⇅    Sleeping too much or not sleeping enough    
    ##  7 MDD4a     SleDec       Sleep⇊    Insomnia nearly every day                   
    ##  8 MDD4b     SleInc       Sleep⇈    Hypersomnia nearly every day                
    ##  9 MDD5      Psyc         Motor⇅    Changes in speed/amount of moving or speaki…
    ## 10 MDD5a     PsycInc      Motor⇈    Psychomotor agitation nearly every day      
    ## 11 MDD5b     PsycDec      Motor⇊    Psychomotor retardation nearly every day    
    ## 12 MDD6      Fatig        Fatigue   Fatigue or loss of energy nearly every day  
    ## 13 MDD7      Guilt        Guilt     Feelings of worthlessness or excessive or i…
    ## 14 MDD8      Conc         Concentr… Diminished ability to think or concentrate,…
    ## 15 MDD9      Sui          Suicidal… Recurrent thoughts of death or suicide or a…

# GenomicSEM covariance structure

``` r
covstruct_prefix <- 'agds_pgc.alspac_ukb.covstruct'
covstruct_r <- file.path('ldsc', paste(covstruct_prefix, 'deparse.R', sep='.'))
covstruct_rds <- file.path('ldsc', paste(covstruct_prefix, 'rds', sep='.'))

symptoms_covstruct <- dget(covstruct_r)

sumstats_prevs <- read_tsv(file.path('ldsc', paste(covstruct_prefix, 'prevs', 'txt', sep='.')))
```

    ## Rows: 24 Columns: 9

    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────
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

    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ## 0.112544816 0.047192446 0.034299744 0.003892157 0.036548569 0.073599244 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ## 0.080469848 0.086912880 0.034935627 0.078103966 0.039893080 0.044634027 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ## 0.395771366 0.058092210 0.055211615 0.057321248 0.049088126 0.032192455

## Common factor

Common factor across symptoms from both cohorts

``` r
commonfactor.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui + PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
"

commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.677 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0371957759493138 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.36605214821313 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 523.5677 104 3.363284e-57 587.5677 0.9532481 0.1836763

``` r
commonfactor.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE
    ## 1           A1 =~  ClinAppDec   0.06235020 0.0854126196136352
    ## 2           A1 =~  ClinAppInc  -0.18187702  0.101563963131273
    ## 4           A1 =~  ClinSleDec  -0.04736174  0.120812396344169
    ## 5           A1 =~  ClinSleInc  -0.16862764  0.166188706209203
    ## 3           A1 =~ ClinPsycInc  -0.18573946  0.110061754269664
    ## 6           A1 =~     ClinSui  -0.64896425  0.108682004115011
    ## 11          A1 =~      PopDep  -0.83967242 0.0451478511326653
    ## 7           A1 =~      PopAnh  -0.95001597 0.0431080771676742
    ## 8           A1 =~   PopAppDec  -0.18232229 0.0806546223278572
    ## 9           A1 =~   PopAppInc  -0.38908141 0.0756533017474777
    ## 14          A1 =~   PopSleDec  -0.57102353 0.0909327734936484
    ## 15          A1 =~   PopSleInc  -0.45958899 0.0952414160323171
    ## 12          A1 =~    PopFatig  -0.65758373 0.0847854931241009
    ## 13          A1 =~    PopGuilt  -0.65220691 0.0800161046980312
    ## 10          A1 =~     PopConc  -0.69494194 0.0835055998717856
    ## 16          A1 =~      PopSui  -0.57893359 0.0929735308368765
    ## 18  ClinAppDec ~~  ClinAppDec   0.99611395  0.219932605919597
    ## 19  ClinAppInc ~~  ClinAppInc   0.96691962  0.322114277288335
    ## 21  ClinSleDec ~~  ClinSleDec   0.99776135  0.545818607899078
    ## 22  ClinSleInc ~~  ClinSleInc   0.97156967  0.842333895318451
    ## 20 ClinPsycInc ~~ ClinPsycInc   0.96550188  0.441518768528219
    ## 23     ClinSui ~~     ClinSui   0.57884564  0.356780235796075
    ## 28      PopDep ~~      PopDep   0.29494933  0.070022263658254
    ## 24      PopAnh ~~      PopAnh   0.09747049 0.0653102308449891
    ## 25   PopAppDec ~~   PopAppDec   0.96675731  0.231302015717083
    ## 26   PopAppInc ~~   PopAppInc   0.84861527  0.153935778231727
    ## 31   PopSleDec ~~   PopSleDec   0.67393197  0.307886253501401
    ## 32   PopSleInc ~~   PopSleInc   0.78877797  0.251632240687266
    ## 29    PopFatig ~~    PopFatig   0.56758503   0.27583394375877
    ## 30    PopGuilt ~~    PopGuilt   0.57462882   0.16470370090682
    ## 27     PopConc ~~     PopConc   0.51705615  0.244965979408896
    ## 33      PopSui ~~      PopSui   0.66483536  0.230285294604523
    ## 17          A1 ~~          A1   1.00000000

Correlation among directional symptoms

``` r
commonfactor_dir.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui + PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
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
    ##   0.415 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0371957759493138 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.36605214821313 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 366.6889 99 2.174698e-32 440.6889 0.9701718 0.1676507

``` r
commonfactor_dir.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppDec   0.07297861 0.0916897298397981 4.260697e-01
    ## 2           A1 =~  ClinAppInc  -0.21337301  0.109824795810457 5.203476e-02
    ## 4           A1 =~  ClinSleDec  -0.04746294  0.125285098422828 7.048143e-01
    ## 5           A1 =~  ClinSleInc  -0.19520115  0.176058457335191 2.675463e-01
    ## 3           A1 =~ ClinPsycInc  -0.20508141  0.118912061322621 8.459020e-02
    ## 6           A1 =~     ClinSui  -0.70847346  0.117829135736631 1.824429e-09
    ## 11          A1 =~      PopDep  -0.68153483 0.0663157454053075 8.936688e-25
    ## 7           A1 =~      PopAnh  -0.80050755 0.0648225982689927 4.919160e-35
    ## 8           A1 =~   PopAppDec  -0.21688957 0.0863212735070294 1.198515e-02
    ## 9           A1 =~   PopAppInc  -0.42986244 0.0800793383070413 7.963504e-08
    ## 14          A1 =~   PopSleDec  -0.64545160 0.0972331335111101 3.175709e-11
    ## 15          A1 =~   PopSleInc  -0.50717594  0.102415474260813 7.339841e-07
    ## 12          A1 =~    PopFatig  -0.72070774 0.0898143355580135 1.020086e-15
    ## 13          A1 =~    PopGuilt  -0.71117351 0.0881806951614938 7.326341e-16
    ## 10          A1 =~     PopConc  -0.76094430 0.0890427143335892 1.276106e-17
    ## 16          A1 =~      PopSui  -0.60581604 0.0980807452059216 6.545226e-10
    ## 19  ClinAppDec ~~  ClinAppInc  -0.37159847  0.218798973096331 8.944069e-02
    ## 23  ClinSleDec ~~  ClinSleInc   0.48211115  0.505936997584638 3.406388e-01
    ## 31      PopDep ~~      PopAnh   0.35654648 0.0971526112855302 2.425800e-04
    ## 28   PopAppDec ~~   PopAppInc  -0.23096326  0.129835893304005 7.525835e-02
    ## 36   PopSleDec ~~   PopSleInc  -0.32597973  0.192127688043739 8.975693e-02
    ## 18  ClinAppDec ~~  ClinAppDec   0.99467234  0.220417858887529 6.401544e-06
    ## 20  ClinAppInc ~~  ClinAppInc   0.95447259  0.323125287753196 3.138105e-03
    ## 22  ClinSleDec ~~  ClinSleDec   0.99774973  0.545881129322718 6.758381e-02
    ## 24  ClinSleInc ~~  ClinSleInc   0.96187035  0.841547170982877 2.530342e-01
    ## 21 ClinPsycInc ~~ ClinPsycInc   0.95794328  0.442328411634951 3.033572e-02
    ## 25     ClinSui ~~     ClinSui   0.49806444  0.368726148369894 1.767659e-01
    ## 32      PopDep ~~      PopDep   0.53550937  0.108311550045554 7.646646e-07
    ## 26      PopAnh ~~      PopAnh   0.35918740  0.106215223414417 7.204338e-04
    ## 27   PopAppDec ~~   PopAppDec   0.95295900  0.230152753615679 3.464737e-05
    ## 29   PopAppInc ~~   PopAppInc   0.81521533  0.156323211030539 1.838621e-07
    ## 35   PopSleDec ~~   PopSleDec   0.58339189  0.312122909285295 6.160715e-02
    ## 37   PopSleInc ~~   PopSleInc   0.74277217  0.250771347133848 3.056999e-03
    ## 33    PopFatig ~~    PopFatig   0.48058030  0.277412468493067 8.320806e-02
    ## 34    PopGuilt ~~    PopGuilt   0.49422867  0.170875512334294 3.823650e-03
    ## 30     PopConc ~~     PopConc   0.42096384  0.248467632386285 9.021962e-02
    ## 38      PopSui ~~      PopSui   0.63298763   0.23258545623317 6.498180e-03
    ## 17          A1 ~~          A1   1.00000000                              NA

Ascertainment-specific factors

``` r
clin_pop.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui 
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
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   3.411 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0371957759493138 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.36605214821313 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_pop.model): A difference greater than .025 was observed pre- and post-
    ## smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_pop.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_pop.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 366.6886 98 1.117434e-32 442.6886 0.9700604 0.1676507

``` r
clin_pop.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 1           A1 =~  ClinAppDec  -0.07297871 0.0981936484125456 4.573541e-01
    ## 2           A1 =~  ClinAppInc   0.21337282  0.138592824578811 1.236668e-01
    ## 4           A1 =~  ClinSleDec   0.04746212  0.127242216169122 7.091438e-01
    ## 5           A1 =~  ClinSleInc   0.19520394  0.198065455365735 3.243573e-01
    ## 3           A1 =~ ClinPsycInc   0.20508192  0.156204268559555 1.892144e-01
    ## 6           A1 =~     ClinSui   0.70847539  0.364653746989258 5.203238e-02
    ## 13          A2 =~      PopDep   0.68153359 0.0662528712687633 8.077338e-25
    ## 9           A2 =~      PopAnh   0.80050694 0.0647520828293692 4.162259e-35
    ## 10          A2 =~   PopAppDec   0.21688953 0.0863204102029384 1.198421e-02
    ## 11          A2 =~   PopAppInc   0.42986204 0.0800689277781399 7.932841e-08
    ## 16          A2 =~   PopSleDec   0.64545210 0.0972284287849818 3.168701e-11
    ## 17          A2 =~   PopSleInc   0.50717671  0.102415773642554 7.340440e-07
    ## 14          A2 =~    PopFatig   0.72070837 0.0898149280977013 1.020484e-15
    ## 15          A2 =~    PopGuilt   0.71117339 0.0881848300021308 7.348983e-16
    ## 12          A2 =~     PopConc   0.76094507 0.0890510953831263 1.284978e-17
    ## 18          A2 =~      PopSui   0.60581634 0.0980630089546475 6.499070e-10
    ## 8           A1 ~~          A2   0.99999987  0.519042427077257 5.402716e-02
    ## 21  ClinAppDec ~~  ClinAppInc  -0.37159802  0.217870325638212 8.808370e-02
    ## 25  ClinSleDec ~~  ClinSleInc   0.48210827  0.505432684571666 3.401572e-01
    ## 33      PopDep ~~      PopAnh   0.35654859 0.0968213661930228 2.309293e-04
    ## 30   PopAppDec ~~   PopAppInc  -0.23096385   0.12983254752786 7.525083e-02
    ## 38   PopSleDec ~~   PopSleInc  -0.32597953  0.192123961436102 8.975010e-02
    ## 20  ClinAppDec ~~  ClinAppDec   0.99467403  0.220297940051199 6.327757e-06
    ## 22  ClinAppInc ~~  ClinAppInc   0.95447022  0.321812777043377 3.017763e-03
    ## 24  ClinSleDec ~~  ClinSleDec   0.99774725  0.545927684354803 6.760687e-02
    ## 26  ClinSleInc ~~  ClinSleInc   0.96175792   0.84312324416208 2.539230e-01
    ## 23 ClinPsycInc ~~ ClinPsycInc   0.95793867  0.449351871502651 3.302092e-02
    ## 27     ClinSui ~~     ClinSui   0.49804222    0.5536599333717 3.683444e-01
    ## 34      PopDep ~~      PopDep   0.53551207  0.108061605018307 7.210388e-07
    ## 28      PopAnh ~~      PopAnh   0.35918877  0.105959623970949 6.993102e-04
    ## 29   PopAppDec ~~   PopAppDec   0.95295901  0.230152660209782 3.464723e-05
    ## 31   PopAppInc ~~   PopAppInc   0.81521871  0.156328878349443 1.840479e-07
    ## 37   PopSleDec ~~   PopSleDec   0.58339209  0.312121592374173 6.160642e-02
    ## 39   PopSleInc ~~   PopSleInc   0.74277206  0.250771199407349 3.056959e-03
    ## 35    PopFatig ~~    PopFatig   0.48057992  0.277414409226789 8.321065e-02
    ## 36    PopGuilt ~~    PopGuilt   0.49423303   0.17087318732714 3.823215e-03
    ## 32     PopConc ~~     PopConc   0.42096336  0.248465287935383 9.021707e-02
    ## 40      PopSui ~~      PopSui   0.63298699  0.232605144639437 6.502683e-03
    ## 7           A1 ~~          A1   1.00000000                              NA
    ## 19          A2 ~~          A2   1.00000000                              NA

``` r
clin_pop_bif.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui 
A2 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui + PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc +  PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
A2 ~~ 1*A2
A ~~ 1*A
A ~~ 0*A1 + 0*A2
A1 ~~ 0*A2
"

clin_pop_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_pop_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.548 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0371957759493138 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.36605214821313 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_pop_bif.model): A difference greater than .025 was observed pre- and post-
    ## smoothing in the genetic covariance matrix. This reflects a large difference
    ## and results should be interpreted with caution!! This can often result from
    ## including low powered traits, and you might consider removing those traits from
    ## the model. If you are going to run a multivariate GWAS we strongly recommend
    ## setting the smooth_check argument to true to check smoothing for each SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_pop_bif.model): A difference greater than .025 was observed pre- and post-
    ## smoothing for Z-statistics in the genetic covariance matrix. This reflects a
    ## large difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

``` r
clin_pop_bif.fit$modelfit
```

    ##       chisq df     p_chisq      AIC       CFI      SRMR
    ## df 314.5882 88 3.18293e-27 410.5882 0.9747516 0.1550796

``` r
clin_pop_bif.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE      p_value
    ## 18          A1 =~  ClinAppDec -0.189551658 0.335152721396251 5.716492e-01
    ## 19          A1 =~  ClinAppInc -0.168305181 0.330576718163868 6.106325e-01
    ## 21          A1 =~  ClinSleDec -2.674363275  4.97253059249543 5.906603e-01
    ## 22          A1 =~  ClinSleInc -0.203069775 0.418079154074667 6.271358e-01
    ## 20          A1 =~ ClinPsycInc -0.079607371  0.17357378993375 6.464714e-01
    ## 23          A1 =~     ClinSui  0.097638163 0.207017114112648 6.371563e-01
    ## 31          A2 =~      PopDep  0.917708557 0.119949476173568 1.997897e-14
    ## 27          A2 =~      PopAnh  0.734650807 0.146527749587895 5.339590e-07
    ## 28          A2 =~   PopAppDec  0.130920187 0.131943945381653 3.210845e-01
    ## 29          A2 =~   PopAppInc  0.101045511 0.128256488013501 4.307934e-01
    ## 34          A2 =~   PopSleDec  0.229781028 0.183873574067791 2.114253e-01
    ## 35          A2 =~   PopSleInc  0.281202091 0.166908124165233 9.203438e-02
    ## 32          A2 =~    PopFatig  0.170316600 0.190142773272727 3.704011e-01
    ## 33          A2 =~    PopGuilt  0.271620022 0.166920785024478 1.036874e-01
    ## 30          A2 =~     PopConc  0.234557534  0.18736287562333 2.106140e-01
    ## 36          A2 =~      PopSui  0.465276644 0.154903631841632 2.667582e-03
    ## 1            A =~  ClinAppDec  0.127359917 0.109064609100708 2.429103e-01
    ## 2            A =~  ClinAppInc -0.307212184 0.133167650256063 2.105730e-02
    ## 4            A =~  ClinSleDec -0.067888380 0.161153073338744 6.735605e-01
    ## 5            A =~  ClinSleInc -0.255577838 0.203071232948277 2.081899e-01
    ## 3            A =~ ClinPsycInc -0.260583471 0.144385864245294 7.110915e-02
    ## 6            A =~     ClinSui -0.875975938 0.186181280999612 2.539214e-06
    ## 11           A =~      PopDep -0.390037335 0.164609765036795 1.781373e-02
    ## 7            A =~      PopAnh -0.587024170  0.15037956092054 9.476156e-05
    ## 8            A =~   PopAppDec -0.118677239  0.11784579976131 3.139025e-01
    ## 9            A =~   PopAppInc -0.445789000 0.118535977753144 1.693793e-04
    ## 14           A =~   PopSleDec -0.566755664 0.164996695181411 5.926563e-04
    ## 15           A =~   PopSleInc -0.354669257 0.160393918145583 2.701899e-02
    ## 12           A =~    PopFatig -0.757257784 0.128571060730549 3.866648e-09
    ## 13           A =~    PopGuilt -0.640276759  0.12251923478968 1.732883e-07
    ## 10           A =~     PopConc -0.741636602 0.144028176838687 2.615326e-07
    ## 16           A =~      PopSui -0.347665140  0.16324152511154 3.319313e-02
    ## 39  ClinAppDec ~~  ClinAppDec  0.947848931 0.256235226128362 2.163068e-04
    ## 40  ClinAppInc ~~  ClinAppInc  0.877294804 0.363792524048673 1.588618e-02
    ## 42  ClinSleDec ~~  ClinSleDec -6.156826772  26.6176823500605 8.170628e-01
    ## 43  ClinSleInc ~~  ClinSleInc  0.893441476 0.859762832712578 2.987272e-01
    ## 41 ClinPsycInc ~~ ClinPsycInc  0.925753089 0.445116766289596 3.754383e-02
    ## 44     ClinSui ~~     ClinSui  0.223127791 0.478133563610096 6.407259e-01
    ## 49      PopDep ~~      PopDep  0.005681948 0.171375050488871 9.735529e-01
    ## 45      PopAnh ~~      PopAnh  0.115691499 0.088012722177136 1.886798e-01
    ## 46   PopAppDec ~~   PopAppDec  0.968776773 0.231204939513021 2.787997e-05
    ## 47   PopAppInc ~~   PopAppInc  0.791061993 0.158350643810544 5.863811e-07
    ## 52   PopSleDec ~~   PopSleDec  0.625990160 0.322088930744976 5.195294e-02
    ## 53   PopSleInc ~~   PopSleInc  0.795134757 0.252611340156758 1.645834e-03
    ## 50    PopFatig ~~    PopFatig  0.397552702 0.277358606217562 1.517557e-01
    ## 51    PopGuilt ~~    PopGuilt  0.516268912 0.180065134732486 4.142247e-03
    ## 48     PopConc ~~     PopConc  0.394952988 0.268290193057489 1.409867e-01
    ## 54      PopSui ~~      PopSui  0.662646363 0.235311075599635 4.861973e-03
    ## 25          A1 ~~          A1  1.000000000                             NA
    ## 38          A2 ~~          A2  1.000000000                             NA
    ## 17           A ~~           A  1.000000000                             NA
    ## 24          A1 ~~           A  0.000000000                             NA
    ## 37          A2 ~~           A  0.000000000                             NA
    ## 26          A1 ~~          A2  0.000000000                             NA

## ADGS-PGC (Clinical)

### Common factor

Common factor model. Allow residual negative correlation between
directional symptoms

``` r
clin_commonfactor.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
A1 ~~ 1*A1
c3a3b > -1
ClinAppDec ~~ c3a3b*ClinAppInc
"
clin_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.829 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.036707898049133 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.29447064580018 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_commonfactor.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## clin_commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
clin_commonfactor.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 2.604371  8 0.9566862 28.60437   1 0.1435086

``` r
clin_commonfactor.fit$results[c(1,2,3,6,7, 9)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE    p_value
    ## 1           A1 =~  ClinAppDec  0.671726995    0.433105425362 0.12091657
    ## 2           A1 =~  ClinAppInc  0.729670895  0.52150444440806 0.16176786
    ## 4           A1 =~  ClinSleDec  0.746906301 0.474573261514982 0.11552771
    ## 5           A1 =~  ClinSleInc  0.612275763 0.461873619105848 0.18496312
    ## 3           A1 =~ ClinPsycInc  0.279998891  0.23717368047761 0.23777824
    ## 6           A1 =~     ClinSui  0.007294712 0.178671946715864 0.96742911
    ## 9   ClinAppDec ~~  ClinAppInc -0.861144952 0.618264546093302 0.16367832
    ## 8   ClinAppDec ~~  ClinAppDec  0.548782969 0.574059361209168 0.33913056
    ## 10  ClinAppInc ~~  ClinAppInc  0.467581888 0.791246118141881 0.55459493
    ## 12  ClinSleDec ~~  ClinSleDec  0.442128708 0.900082705615737 0.62325288
    ## 13  ClinSleInc ~~  ClinSleInc  0.625191085  1.61291824164948 0.69833329
    ## 11 ClinPsycInc ~~ ClinPsycInc  0.921603214 0.470933722825115 0.05035139
    ## 14     ClinSui ~~     ClinSui  0.999947372 0.363500563648505 0.00594352
    ## 7           A1 ~~          A1  1.000000000                           NA

``` r
fit_graph <- function(results, ...) {

  results_sort <- results %>% arrange(lhs, rhs)

  node_names <- unique(c(results_sort$lhs, results_sort$rhs))
  
  node_idx <- seq_along(node_names)
  names(node_idx) <- node_names
  
  graph <- create_graph(
    nodes_df=create_node_df(n=length(node_names),
                            label=node_labels[node_names], 
                            shape='oval', width=1,
                            fillcolor=node_colors[node_names],
                            fontcolor='black'),
    edges_df=create_edge_df(from=node_idx[results_sort$lhs],
                            to=node_idx[results_sort$rhs],
                            label=round(results_sort$STD_Genotype, 2),
                            penwidth=0.3+abs(2*results_sort$STD_Genotype),
                            dir=edge_dir[results_sort$op]),
    attr_theme="tb")
  
  return(graph)

}

# render_fit <- function(results) render_graph(fit_graph(results))

# render_fit(pgc_commonfactor.fit$results)

# pgc_commonfactor.graph <- fit_graph(pgc_commonfactor.fit$results)

add_rank_same <- function(gv, top, bottom) {

  # add in block to specify node ranks
  gv.list <- str_split(gv, '\n\n')[[1]]
  
  # move the nodes/edges element to the end
  gv.list[6] <- gv.list[5]
  
  # add manually made ranks
  gv.list[5] <- paste("{rank=same",
  paste(paste0("'", top, "'"), collapse=' '),
  "}\n{rank=same",
  paste(paste0("'", bottom, "'"), collapse=' '),
  "}")
  
  rank.gv <- paste(gv.list, collapse='\n\n')

}


# pgc_commonfactor.gv <- add_rank_same(generate_dot(pgc_commonfactor.graph), 1, 2:5)
# grViz(pgc_commonfactor.gv)

# # output as a GraphViz dot file. Replace single quotes with double quotes 
# # as that's what the command line utility expects
# cat(str_replace_all(pgc_commonfactor_rank.gv, "'", '"'), file='mdd-symptom-gsem_files/pgc_commonfactor.gv')
```

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
    ##     0.2 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 59.80875 35 0.005594967 99.80875 0.9928165 0.1257256

``` r
pop_commonfactor.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep   0.86302758 0.0466869423189087
    ## 1         A1 =~    PopAnh   0.99707952 0.0429873302659936
    ## 2         A1 =~ PopAppDec   0.17061418 0.0881755804155387
    ## 3         A1 =~ PopAppInc   0.40619378 0.0809099038334041
    ## 8         A1 =~ PopSleDec   0.57632446 0.0948331691360743
    ## 9         A1 =~ PopSleInc   0.52022007  0.101746593406534
    ## 6         A1 =~  PopFatig   0.67199535  0.093755347882219
    ## 7         A1 =~  PopGuilt   0.62768566  0.080178150443178
    ## 4         A1 =~   PopConc   0.70163978 0.0903611211173149
    ## 10        A1 =~    PopSui   0.58716232 0.0983552072249836
    ## 16    PopDep ~~    PopDep   0.25518298 0.0738079091843754
    ## 12    PopAnh ~~    PopAnh   0.00583183 0.0676894778042077
    ## 13 PopAppDec ~~ PopAppDec   0.97089128  0.251433696769557
    ## 14 PopAppInc ~~ PopAppInc   0.83500651  0.166782794848229
    ## 19 PopSleDec ~~ PopSleDec   0.66785177  0.331077725469588
    ## 20 PopSleInc ~~ PopSleInc   0.72936865  0.284423005895444
    ## 17  PopFatig ~~  PopFatig   0.54842296  0.322517150649473
    ## 18  PopGuilt ~~  PopGuilt   0.60601068   0.16600207289361
    ## 15   PopConc ~~   PopConc   0.50770129  0.268448846854788
    ## 21    PopSui ~~    PopSui   0.65524003  0.253105388878833
    ## 11        A1 ~~        A1   1.00000000

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
    ##   0.214 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df  p_chisq      AIC       CFI     SRMR
    ## df 40.57809 34 0.202952 82.57809 0.9980953 0.117905

``` r
pop_commonfactor_gating.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 5         A1 =~    PopDep    0.7357950 0.0783402055161693
    ## 1         A1 =~    PopAnh    0.8736866   0.07343774915562
    ## 2         A1 =~ PopAppDec    0.1846453 0.0936391312342128
    ## 3         A1 =~ PopAppInc    0.4363186 0.0877965311929842
    ## 8         A1 =~ PopSleDec    0.6278678  0.103964412847252
    ## 9         A1 =~ PopSleInc    0.5482409  0.109301866858449
    ## 6         A1 =~  PopFatig    0.7276567  0.101483918810308
    ## 7         A1 =~  PopGuilt    0.6808367  0.088793518970862
    ## 4         A1 =~   PopConc    0.7679950 0.0985288521213313
    ## 10        A1 =~    PopSui    0.6213229  0.105725469648134
    ## 16    PopDep ~~    PopAnh    0.2809443   0.11916598574347
    ## 17    PopDep ~~    PopDep    0.4586061  0.130383037888472
    ## 12    PopAnh ~~    PopAnh    0.2366716  0.130111227092262
    ## 13 PopAppDec ~~ PopAppDec    0.9659046  0.251348400712642
    ## 14 PopAppInc ~~ PopAppInc    0.8096264  0.171640389390296
    ## 20 PopSleDec ~~ PopSleDec    0.6057816  0.333924086088153
    ## 21 PopSleInc ~~ PopSleInc    0.6994335   0.28533000625151
    ## 18  PopFatig ~~  PopFatig    0.4705159  0.324391546183262
    ## 19  PopGuilt ~~  PopGuilt    0.5364618  0.171168790297976
    ## 15   PopConc ~~   PopConc    0.4101823  0.270813260515178
    ## 22    PopSui ~~    PopSui    0.6139586  0.258378911612517
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
    ##    0.21 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 39.86281 33 0.1912613 83.86281 0.9980128 0.1155375

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
    ##    0.21 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 38.33264 33 0.2403506 82.33264 0.9984559 0.1063218

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
    ##   1.435 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 38.15577 32 0.2097694 84.15577 0.9982176 0.1036837

### Cognitive-Mood-Neuroveg (Kendler Neale) model

``` r
pop_cog_mood_neuroveg.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
"
pop_cog_mood_neuroveg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg.model)
```

    ## [1] "Running primary model"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = pop_cog_mood_neuroveg.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

    ##          lhs op       rhs Unstandardized_Estimate
    ## 1         A1 =~  PopGuilt            1.587337e-01
    ## 2         A1 =~   PopConc            1.687508e-01
    ## 3         A1 =~    PopSui            1.074149e-01
    ## 4         A2 =~    PopDep            3.274805e+00
    ## 5         A2 =~    PopAnh            4.037000e+00
    ## 6         A2 =~  PopGuilt           -2.607563e-04
    ## 7         A3 =~ PopSleDec            1.300047e-01
    ## 8         A3 =~ PopSleInc            1.177593e-01
    ## 9         A3 =~  PopFatig            1.657070e-01
    ## 10        A3 =~ PopAppDec            3.368885e-02
    ## 11        A3 =~ PopAppInc            1.183733e-01
    ## 15    PopDep ~~    PopAnh           -1.314312e+01
    ## 16 PopSleDec ~~ PopSleInc           -1.589248e-02
    ## 17  PopGuilt ~~  PopGuilt            3.214533e-02
    ## 18   PopConc ~~   PopConc            2.302874e-02
    ## 19    PopSui ~~    PopSui            2.078427e-02
    ## 20    PopDep ~~    PopDep           -1.064385e+01
    ## 21    PopAnh ~~    PopAnh           -1.621045e+01
    ## 22 PopSleDec ~~ PopSleDec            2.510722e-02
    ## 23 PopSleInc ~~ PopSleInc            3.099807e-02
    ## 24  PopFatig ~~  PopFatig            2.807104e-02
    ## 25 PopAppDec ~~ PopAppDec            3.404485e-02
    ## 26 PopAppInc ~~ PopAppInc            6.413941e-02
    ## 27        A1 ~~        A2            6.594259e-02
    ## 28        A1 ~~        A3            1.116619e+00
    ## 29        A2 ~~        A3            6.525973e-02

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
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
a13 > -1.0
A1 ~~ a13*A3
"
pop_cog_mood_neuroveg_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_cog_mood_neuroveg_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##          lhs op       rhs Unstandardized_Estimate          SE
    ## 1         A1 =~  PopGuilt             -0.01766959 0.050710091
    ## 2         A1 =~   PopConc              0.11413473 0.057050840
    ## 3         A1 =~    PopSui              0.07335583 0.037405844
    ## 4         A2 =~    PopDep              0.20702603 0.023941247
    ## 5         A2 =~    PopAnh              0.25512942 0.025307876
    ## 6         A2 =~  PopGuilt              0.18864697 0.071049368
    ## 7         A3 =~ PopSleDec              0.12975876 0.025081029
    ## 8         A3 =~ PopSleInc              0.11800587 0.025157860
    ## 9         A3 =~  PopFatig              0.16516547 0.029440994
    ## 10        A3 =~ PopAppDec              0.03355740 0.017206475
    ## 11        A3 =~ PopAppInc              0.11887959 0.027011792
    ## 15    PopDep ~~    PopAnh              0.02445368 0.011416891
    ## 16 PopSleDec ~~ PopSleInc             -0.01589552 0.009423030
    ## 17    PopAnh ~~    PopAnh              0.02182717 0.013375792
    ## 18        A1 ~~        A3              1.74388258 0.802541426
    ## 19  PopGuilt ~~  PopGuilt              0.03168680 0.010475344
    ## 20   PopConc ~~   PopConc              0.03847883 0.016947832
    ## 21    PopSui ~~    PopSui              0.02694116 0.009524750
    ## 22    PopDep ~~    PopDep              0.03763734 0.011401996
    ## 23 PopSleDec ~~ PopSleDec              0.02517111 0.013872450
    ## 24 PopSleInc ~~ PopSleInc              0.03093992 0.012939483
    ## 25  PopFatig ~~  PopFatig              0.02825020 0.017860917
    ## 26 PopAppDec ~~ PopAppDec              0.03405369 0.008836591
    ## 27 PopAppInc ~~ PopAppInc              0.06401928 0.013391865
    ## 28        A1 ~~        A2              1.53756179 0.735616965
    ## 29        A2 ~~        A3              1.03312208 0.167529657

``` r
pop_cog_mood_neuroveg_constr.fit$modelfit
```

    ## NULL

``` r
pop_cog_mood_neuroveg_constr.fit$results[c(1, 2, 3, 6, 7)]
```

    ## NULL

### Two-factor models

[Elhai Psychiat Res
2012](https://www.sciencedirect.com/science/article/pii/S0165178112002685)
compared 3 two-factor models

### Psychological-Somatic (Elhai Model 2a)

[Kruse Rehab Psychol
2008](https://psycnet.apa.org/record/2008-17022-011), [Kruse Arch Psys
Med Rehab
2010](https://www.sciencedirect.com/science/article/pii/S0003999310002443):

> the 2-factor solution with 3 somatic items (sleep disturbance, poor
> energy, appetite change) was a better solution than either a
> unidimensional model or 2-factor model that included psychomotor
> retardation as a fourth somatic item

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
    ## 1         A1 =~    PopDep              0.20414189 0.021805619
    ## 2         A1 =~    PopAnh              0.25222889 0.021822440
    ## 3         A1 =~  PopGuilt              0.16343160 0.021214808
    ## 4         A1 =~   PopConc              0.17576041 0.022291015
    ## 5         A1 =~    PopSui              0.11160297 0.018944333
    ## 6         A2 =~ PopAppDec              0.03170027 0.016348274
    ## 7         A2 =~ PopAppInc              0.11236667 0.027040807
    ## 8         A2 =~ PopSleDec              0.11837466 0.025122228
    ## 9         A2 =~ PopSleInc              0.10665433 0.024434683
    ## 10        A2 =~  PopFatig              0.15586222 0.030523349
    ## 13    PopDep ~~    PopAnh              0.02578163 0.009731854
    ## 14    PopDep ~~    PopDep              0.03882320 0.010197941
    ## 15    PopAnh ~~    PopAnh              0.02329877 0.011290948
    ## 16  PopGuilt ~~  PopGuilt              0.03062643 0.009824874
    ## 17   PopConc ~~   PopConc              0.02061386 0.013897763
    ## 18    PopSui ~~    PopSui              0.01986703 0.008349317
    ## 19 PopAppDec ~~ PopAppDec              0.03417490 0.008844264
    ## 20 PopAppInc ~~ PopAppInc              0.06552539 0.013258849
    ## 21 PopSleDec ~~ PopSleDec              0.02799584 0.013850274
    ## 22 PopSleInc ~~ PopSleInc              0.03349014 0.013092088
    ## 23  PopFatig ~~  PopFatig              0.03123678 0.017848818
    ## 24        A1 ~~        A2              1.13705716 0.181376019

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
    ##   1.492 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI     SRMR
    ## df 40.57811 33 0.1709599 84.57811 0.9978057 0.117905

``` r
pop_psych_soma_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 3         A1 =~    PopDep    0.7357949 0.0792508320063668
    ## 1         A1 =~    PopAnh    0.8736870 0.0765156886919686
    ## 4         A1 =~  PopGuilt    0.6808364 0.0890119915804785
    ## 2         A1 =~   PopConc    0.7679945 0.0986678324379925
    ## 5         A1 =~    PopSui    0.6213223  0.105717653177007
    ## 8         A2 =~ PopAppDec    0.1846452 0.0939232863227464
    ## 9         A2 =~ PopAppInc    0.4363188 0.0977850825190984
    ## 11        A2 =~ PopSleDec    0.6278684  0.122381900618668
    ## 12        A2 =~ PopSleInc    0.5482416  0.117355804363177
    ## 10        A2 =~  PopFatig    0.7276581  0.126881612846273
    ## 18    PopDep ~~    PopAnh    0.2809439  0.121956958964034
    ## 7         A1 ~~        A2    0.9999999  0.141132424247682
    ## 19    PopDep ~~    PopDep    0.4586056  0.131677951429562
    ## 14    PopAnh ~~    PopAnh    0.2366710  0.136399251784379
    ## 21  PopGuilt ~~  PopGuilt    0.5364626  0.171536375923612
    ## 17   PopConc ~~   PopConc    0.4101866  0.270222452822113
    ## 24    PopSui ~~    PopSui    0.6139588  0.258567177531364
    ## 15 PopAppDec ~~ PopAppDec    0.9659055  0.251477286126913
    ## 16 PopAppInc ~~ PopAppInc    0.8096257  0.172906314831802
    ## 22 PopSleDec ~~ PopSleDec    0.6057720  0.330123576753494
    ## 23 PopSleInc ~~ PopSleInc    0.6994251  0.292864388187009
    ## 20  PopFatig ~~  PopFatig    0.4705000  0.323305879957067
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
    ##   0.234 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC CFI       SRMR
    ## df 19.30107 24 0.7357486 81.30107   1 0.06566412

``` r
pop_psych_soma_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 14        A1 =~    PopDep   0.40288594  0.206116952672584
    ## 12        A1 =~    PopAnh   0.14774143  0.171426467168234
    ## 15        A1 =~  PopGuilt   0.26332202  0.180894428119556
    ## 13        A1 =~   PopConc  -0.29735958  0.230311095399884
    ## 16        A1 =~    PopSui   0.75735904  0.370827635635892
    ## 20        A2 =~ PopAppDec   0.11914964  0.200027515899732
    ## 21        A2 =~ PopAppInc   0.17686378  0.175194228919508
    ## 23        A2 =~ PopSleDec   0.69647576  0.358319053772263
    ## 24        A2 =~ PopSleInc  -0.63559168  0.305039925090693
    ## 22        A2 =~  PopFatig  -0.36873796   0.23212546267315
    ## 5          A =~    PopDep   0.65273465 0.0879947284256678
    ## 1          A =~    PopAnh   0.82379014 0.0707052725383307
    ## 7          A =~  PopGuilt   0.63722299  0.101702352971688
    ## 4          A =~   PopConc   0.90658878  0.119618464922695
    ## 10         A =~    PopSui   0.48424380  0.145463123425109
    ## 2          A =~ PopAppDec   0.19156965 0.0967301199876218
    ## 3          A =~ PopAppInc   0.45823790 0.0910152851057044
    ## 8          A =~ PopSleDec   0.70126733  0.110643518245243
    ## 9          A =~ PopSleInc   0.59392098  0.113075791333711
    ## 6          A =~  PopFatig   0.77320673   0.10803272361709
    ## 31    PopDep ~~    PopAnh   0.32655906  0.122413614768284
    ## 32    PopDep ~~    PopDep   0.41162035  0.181202575360938
    ## 27    PopAnh ~~    PopAnh   0.29954193  0.117151140813626
    ## 34  PopGuilt ~~  PopGuilt   0.52460843  0.176553607720371
    ## 30   PopConc ~~   PopConc   0.08967354  0.332437407215755
    ## 37    PopSui ~~    PopSui   0.19191504  0.571443745498145
    ## 28 PopAppDec ~~ PopAppDec   0.94910581  0.251104087294931
    ## 29 PopAppInc ~~ PopAppInc   0.75873714  0.195977289034318
    ## 35 PopSleDec ~~ PopSleDec   0.02314689  0.602652076556542
    ## 36 PopSleInc ~~ PopSleInc   0.24328096  0.467301901728063
    ## 33  PopFatig ~~  PopFatig   0.26618312   0.35990838300219
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
    ##   0.231 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 41.47245 33 0.1478819 85.47245 0.9975468 0.1149995

``` r
pop_psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 2         A1 =~    PopDep    0.7649524 0.0877459944870521
    ## 1         A1 =~    PopAnh    0.9061496 0.0904406172996843
    ## 3         A1 =~  PopGuilt    0.6746060  0.089431525001776
    ## 4         A1 =~    PopSui    0.6202713   0.10612052931221
    ## 7         A2 =~ PopAppDec    0.1909809 0.0970905646167474
    ## 8         A2 =~ PopAppInc    0.4556838 0.0946454958263406
    ## 11        A2 =~ PopSleDec    0.6552295  0.113786093357341
    ## 12        A2 =~ PopSleInc    0.5689441  0.114775073709934
    ## 10        A2 =~  PopFatig    0.7614472   0.11332848545178
    ## 9         A2 =~   PopConc    0.8091796  0.108009766926521
    ## 18    PopDep ~~    PopAnh    0.2306371  0.147433729878873
    ## 19    PopDep ~~    PopDep    0.4148477  0.150292475556504
    ## 14    PopAnh ~~    PopAnh    0.1788934  0.168660935800506
    ## 21  PopGuilt ~~  PopGuilt    0.5449090   0.17213785227647
    ## 24    PopSui ~~    PopSui    0.6152646  0.258901774655008
    ## 15 PopAppDec ~~ PopAppDec    0.9635272  0.251249718530321
    ## 16 PopAppInc ~~ PopAppInc    0.7923520  0.174299658453204
    ## 22 PopSleDec ~~ PopSleDec    0.5706745   0.33269403302518
    ## 23 PopSleInc ~~ PopSleInc    0.6763051  0.287503372559605
    ## 20  PopFatig ~~  PopFatig    0.4201988  0.320052579961974
    ## 17   PopConc ~~   PopConc    0.3452282  0.277230094128007
    ## 6         A1 ~~        A2    0.9053697  0.107973142187666
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
    ##    0.24 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 16.23925 23 0.8448369 80.23925   1 0.07369998

``` r
pop_psych_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 13        A1 =~    PopDep    0.4780922  0.222719038602185
    ## 12        A1 =~    PopAnh    0.2018480  0.182337937522341
    ## 14        A1 =~  PopGuilt    0.3399823  0.175771966132601
    ## 15        A1 =~    PopSui    0.7683686  0.367447068366923
    ## 19        A2 =~ PopAppDec   -0.1397065  0.204427000742679
    ## 20        A2 =~ PopAppInc    0.5471587  0.349724553494874
    ## 23        A2 =~ PopSleDec    0.2431083  0.252899683729667
    ## 24        A2 =~ PopSleInc   -0.3203945  0.231981553591314
    ## 22        A2 =~  PopFatig   -0.4029897  0.318446863654678
    ## 21        A2 =~   PopConc    0.2683295  0.237398361918738
    ## 5          A =~    PopDep    0.6044214 0.0927993800733505
    ## 1          A =~    PopAnh    0.8108809 0.0833271348149793
    ## 7          A =~  PopGuilt    0.6028374  0.105367030973146
    ## 10         A =~    PopSui    0.4222852  0.143066444018353
    ## 2          A =~ PopAppDec    0.1984550 0.0999354785427738
    ## 3          A =~ PopAppInc    0.4775928  0.101548721500163
    ## 8          A =~ PopSleDec    0.6960469  0.115740443391928
    ## 9          A =~ PopSleInc    0.6418012  0.121180946904928
    ## 6          A =~  PopFatig    0.8235212  0.122275976461565
    ## 4          A =~   PopConc    0.8399459   0.10892872022568
    ## 31    PopDep ~~    PopAnh    0.3371827  0.123843029439672
    ## 36 PopSleDec ~~ PopSleInc   -0.3822674  0.245234473289357
    ## 32    PopDep ~~    PopDep    0.4061026  0.201642443541943
    ## 27    PopAnh ~~    PopAnh    0.3017296    0.1210306041663
    ## 34  PopGuilt ~~  PopGuilt    0.5209991  0.180172898040142
    ## 38    PopSui ~~    PopSui    0.2312851  0.565284868312736
    ## 28 PopAppDec ~~ PopAppDec    0.9410980  0.252177717619961
    ## 29 PopAppInc ~~ PopAppInc    0.4725222  0.429665302683225
    ## 35 PopSleDec ~~ PopSleDec    0.4564166  0.343348479097223
    ## 37 PopSleInc ~~ PopSleInc    0.4854381  0.329096313517553
    ## 33  PopFatig ~~  PopFatig    0.1594121  0.438514742019357
    ## 30   PopConc ~~   PopConc    0.2224893  0.319717691575036
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
    ##    0.22 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df  p_chisq      AIC CFI      SRMR
    ## df 30.35261 33 0.599588 74.35261   1 0.1123386

``` r
pop_affect_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7816755   0.08173494743752
    ## 2         A1 =~  PopGuilt    0.7371765 0.0983352075607999
    ## 3         A1 =~    PopSui    0.6837611   0.11602439184229
    ## 6         A2 =~    PopAnh    0.8801197 0.0732212061890407
    ## 8         A2 =~ PopAppInc    0.4522550 0.0902475690044402
    ## 7         A2 =~ PopAppDec    0.1885547  0.096620707897308
    ## 12        A2 =~ PopSleInc    0.5659892  0.111444719555614
    ## 11        A2 =~ PopSleDec    0.6438611   0.10767117031592
    ## 10        A2 =~  PopFatig    0.7562461   0.10420501208481
    ## 9         A2 =~   PopConc    0.7982991  0.100661106758565
    ## 18    PopDep ~~    PopAnh    0.3535100  0.108409798075454
    ## 19    PopDep ~~    PopDep    0.3889835  0.143675486141461
    ## 21  PopGuilt ~~  PopGuilt    0.4565711  0.182906239703438
    ## 24    PopSui ~~    PopSui    0.5324704   0.26533908840485
    ## 14    PopAnh ~~    PopAnh    0.2253898  0.130825469326711
    ## 16 PopAppInc ~~ PopAppInc    0.7954664  0.173231659256943
    ## 15 PopAppDec ~~ PopAppDec    0.9644478   0.25131756737418
    ## 23 PopSleInc ~~ PopSleInc    0.6796565  0.285924221840464
    ## 22 PopSleDec ~~ PopSleDec    0.5854421  0.332410763850476
    ## 20  PopFatig ~~  PopFatig    0.4280928  0.321064642567528
    ## 17   PopConc ~~   PopConc    0.3627187  0.269143534316504
    ## 5         A1 ~~        A2    0.8289468 0.0648763623257505
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
    ##   0.242 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ##   2.497 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150265 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975606 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 17.14818 24 0.8422739 79.14818   1 0.06815124

``` r
pop_affect_veg_bif_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs  STD_Genotype    STD_Genotype_SE
    ## 12        A1 =~    PopDep  0.4511769132  0.265062523011915
    ## 13        A1 =~  PopGuilt  0.2836909544  0.165843436483842
    ## 14        A1 =~    PopSui  0.6637616561  0.367885574707484
    ## 18        A2 =~    PopAnh  0.2100429911   0.11169794603305
    ## 20        A2 =~ PopAppInc -0.1124239673  0.175249222738613
    ## 19        A2 =~ PopAppDec -0.1227738624  0.191515885789678
    ## 24        A2 =~ PopSleInc  0.6422547829  0.256241943995125
    ## 23        A2 =~ PopSleDec -0.6929003958  0.348172290544923
    ## 22        A2 =~  PopFatig  0.3784454165  0.224544821024928
    ## 21        A2 =~   PopConc -0.0838511756  0.190901202768679
    ## 5          A =~    PopDep  0.6273326433 0.0863700198999284
    ## 7          A =~  PopGuilt  0.6412093724 0.0944832941405053
    ## 10         A =~    PopSui  0.5107517510  0.108481042975323
    ## 1          A =~    PopAnh  0.8565593292 0.0813452447541233
    ## 3          A =~ PopAppInc  0.4670050684 0.0914896153042493
    ## 2          A =~ PopAppDec  0.2039859919 0.0994409804933703
    ## 9          A =~ PopSleInc  0.5349462521  0.127271480503417
    ## 8          A =~ PopSleDec  0.7645687054  0.125881565378139
    ## 6          A =~  PopFatig  0.7385208038  0.112053926485566
    ## 4          A =~   PopConc  0.8230037005  0.105970214277626
    ## 31    PopDep ~~    PopAnh  0.3864508634   0.11171851231436
    ## 35 PopSleDec ~~ PopSleDec  0.0009999034  0.597977194848188
    ## 36 PopSleInc ~~ PopSleInc  0.3013390581  0.387422306680823
    ## 32    PopDep ~~    PopDep  0.4028932827  0.227469858410755
    ## 34  PopGuilt ~~  PopGuilt  0.5083699934  0.178850511516754
    ## 37    PopSui ~~    PopSui  0.2985526191  0.514636104805783
    ## 27    PopAnh ~~    PopAnh  0.2221879994  0.140480931368019
    ## 29 PopAppInc ~~ PopAppInc  0.7692667939   0.19003530356453
    ## 28 PopAppDec ~~ PopAppDec  0.9433160170  0.250727912200488
    ## 33  PopFatig ~~  PopFatig  0.3113657657  0.340687373471766
    ## 30   PopConc ~~   PopConc  0.3156335965  0.276353047436225
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
     "2c(ii)"=list(name="Affect-Neuroveg (BiF)", model=pop_affect_veg_bif_constr.fit)
     )

pop_model_fits <- 
data.frame(Model=names(pop_model_list),
           Name=sapply(pop_model_list, function(m) m$name)) %>%
bind_cols(
bind_rows(
lapply(pop_model_list, function(m) m$model$modelfit)
))
rownames(pop_model_fits) <- NULL

pop_model_fits %>%
select(-chisq, -df) %>%
mutate(dAIC=AIC-min(AIC)) %>%
mutate_if(is.numeric, ~round(., 3))
```

    ##     Model                  Name p_chisq    AIC   CFI  SRMR   dAIC
    ## 1      1a                Common   0.006 99.809 0.993 0.126 25.456
    ## 2      1b       Common (gating)   0.203 82.578 0.998 0.118  8.225
    ## 3      1c          Common (App)   0.191 83.863 0.998 0.116  9.510
    ## 4      1d          Common (Sle)   0.240 82.333 0.998 0.106  7.980
    ## 5     1er      Common (App,Sle)   0.210 84.156 0.998 0.104  9.803
    ## 6   2a(i)         Psych-Somatic   0.171 84.578 0.998 0.118 10.225
    ## 7  2a(ii)   Psych-Somatic (BiF)   0.736 81.301 1.000 0.066  6.948
    ## 8   2b(i)        Psych-Neuroveg   0.148 85.472 0.998 0.115 11.120
    ## 9  2b(ii)  Psych-Neuroveg (BiF)   0.845 80.239 1.000 0.074  5.887
    ## 10  2c(i)       Affect-Neuroveg   0.600 74.353 1.000 0.112  0.000
    ## 11 2c(ii) Affect-Neuroveg (BiF)   0.842 79.148 1.000 0.068  4.796

# Clinical and population factors

``` r
clin_pop_affect_veg.model <- "
Soma =~ NA*ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc
Soma ~~ 1*Soma
Affect =~ NA*PopDep + PopGuilt + PopSui
Veg =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
Affect ~~ 1*Affect
Veg ~~ 1*Veg
PopDep ~~ PopAnh
ClinSui ~~ Soma + Affect + Veg
"
clin_pop_affect_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=clin_pop_affect_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   0.476 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0289190810466149 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  1.11962312826327 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 166.5296 84 2.179951e-07 238.5296 0.9895636 0.1492129

``` r
clin_pop_affect_veg.fit$results[c(1,2,3,6,7,9)]
```

    ##            lhs op         rhs STD_Genotype    STD_Genotype_SE      p_value
    ## 23        Soma =~  ClinAppInc   0.63311669  0.249315964519296 1.110364e-02
    ## 25        Soma =~  ClinSleDec   0.39997011  0.285148917126024 1.607154e-01
    ## 26        Soma =~  ClinSleInc   0.52502307  0.323796141121201 1.049170e-01
    ## 24        Soma =~ ClinPsycInc   0.55124682   0.25174645840846 2.854691e-02
    ## 1       Affect =~      PopDep   0.74223676 0.0711883555737211 1.878874e-25
    ## 2       Affect =~    PopGuilt   0.75641804 0.0966276584782944 4.949729e-15
    ## 3       Affect =~      PopSui   0.67514483   0.10836901095603 4.662396e-10
    ## 31         Veg =~      PopAnh   0.81262245 0.0662731036391439 1.453603e-34
    ## 33         Veg =~   PopAppInc   0.42848190 0.0822715940548294 1.907346e-07
    ## 32         Veg =~   PopAppDec   0.20059537 0.0874704946013077 2.183091e-02
    ## 37         Veg =~   PopSleInc   0.49739241  0.103163472766313 1.425547e-06
    ## 36         Veg =~   PopSleDec   0.63042200 0.0996274695454923 2.486952e-10
    ## 35         Veg =~    PopFatig   0.74024906  0.093079382192492 1.822186e-15
    ## 34         Veg =~     PopConc   0.79221142 0.0918017359133403 6.156534e-18
    ## 16      PopDep ~~      PopAnh   0.40446825 0.0925563300545468 1.242495e-05
    ## 28        Soma ~~     ClinSui   0.12275706  0.255142309699262 6.304176e-01
    ## 5       Affect ~~     ClinSui   0.73338255  0.135314986147485 5.966047e-08
    ## 38         Veg ~~     ClinSui   0.65747729  0.132218363326822 6.603947e-07
    ## 7   ClinAppInc ~~  ClinAppInc   0.59916677  0.444388081648588 1.775708e-01
    ## 9   ClinSleDec ~~  ClinSleDec   0.84003033  0.628750873516877 1.815426e-01
    ## 10  ClinSleInc ~~  ClinSleInc   0.72434295  0.882588601524776 4.118152e-01
    ## 8  ClinPsycInc ~~ ClinPsycInc   0.69613007  0.544701812792597 2.012487e-01
    ## 17      PopDep ~~      PopDep   0.44908399  0.124940164766301 3.251602e-04
    ## 19    PopGuilt ~~    PopGuilt   0.42783289  0.180789594345188 1.795970e-02
    ## 22      PopSui ~~      PopSui   0.54418281  0.241437218582983 2.420195e-02
    ## 12      PopAnh ~~      PopAnh   0.33964474  0.111071711534981 2.229033e-03
    ## 14   PopAppInc ~~   PopAppInc   0.81640359  0.154946833406552 1.372301e-07
    ## 13   PopAppDec ~~   PopAppDec   0.95975934  0.231342547177404 3.344220e-05
    ## 21   PopSleInc ~~   PopSleInc   0.75260306  0.252689112396915 2.897887e-03
    ## 20   PopSleDec ~~   PopSleDec   0.60256881  0.295152478043386 4.119623e-02
    ## 18    PopFatig ~~    PopFatig   0.45202976  0.275301216035343 1.006004e-01
    ## 15     PopConc ~~     PopConc   0.37239827  0.251346509715038 1.384421e-01
    ## 11     ClinSui ~~     ClinSui   1.00000043  0.328014352733035 2.298751e-03
    ## 27        Soma ~~      Affect   0.08783379  0.138186258889435 5.250331e-01
    ## 30        Soma ~~         Veg   0.46797492  0.180711881874069 9.608197e-03
    ## 6       Affect ~~         Veg   0.82831046 0.0664911442651041 1.273672e-35
    ## 29        Soma ~~        Soma   1.00000000                              NA
    ## 4       Affect ~~      Affect   1.00000000                              NA
    ## 39         Veg ~~         Veg   1.00000000                              NA

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
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.005       0.173       0.005       0.470       0.918       0.005 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec           0.989  -0.130 
    ## ClinAppInc   0.767  -0.450   0.191 
    ## ClinSleDec   0.940   0.285  -0.174 
    ## ClinSleInc   0.556   0.104   0.459 
    ## ClinPsycInc          0.253   0.130 
    ## ClinSui                      0.997 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      1.782   1.337   1.305
    ## Proportion Var   0.297   0.223   0.217
    ## Cumulative Var   0.297   0.520   0.737
    ## 
    ## The degrees of freedom for the model is 0 and the fit was 0.582

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
    ##     PopDep     PopAnh  PopAppDec  PopAppInc  PopSleDec  PopSleInc PopPsycInc 
    ##      0.005      0.005      0.811      0.490      0.773      0.437      0.005 
    ## PopPsycDec   PopFatig   PopGuilt    PopConc     PopSui 
    ##      0.416      0.507      0.614      0.497      0.571 
    ## 
    ## Loadings:
    ##            Factor1 Factor2 Factor3
    ## PopDep      0.944   0.320         
    ## PopAnh      0.746   0.639  -0.173 
    ## PopAppDec   0.170   0.104   0.386 
    ## PopAppInc   0.194   0.216  -0.652 
    ## PopSleDec   0.297   0.361         
    ## PopSleInc           0.745         
    ## PopPsycInc          0.208   0.975 
    ## PopPsycDec  0.354           0.675 
    ## PopFatig    0.175   0.658  -0.171 
    ## PopGuilt    0.421   0.408   0.204 
    ## PopConc     0.208   0.672         
    ## PopSui      0.629           0.161 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      2.381   2.359   2.128
    ## Proportion Var   0.198   0.197   0.177
    ## Cumulative Var   0.198   0.395   0.572
    ## 
    ## The degrees of freedom for the model is 33 and the fit was 8.1788

## All symptoms

Check eigen values of the correlation matrix

``` r
symptoms_cov_pd.eigen <- eigen(cov2cor(symptoms_cov_pd))

signif(symptoms_cov_pd.eigen$values, 3)
```

    ##  [1] 4.98e+00 3.14e+00 2.60e+00 1.90e+00 1.33e+00 1.21e+00 9.80e-01 4.83e-01
    ##  [9] 4.58e-01 3.89e-01 2.73e-01 2.58e-01 1.27e-07 1.24e-07 8.63e-08 7.59e-08
    ## [17] 6.21e-08 3.37e-08

``` r
plot(eigen(cov2cor(symptoms_cov_pd))$values, ylab='Eigenvalue')
lines(eigen(cov2cor(symptoms_cov_pd))$values)
abline(1, 0, col='red')
```

![](mdd-symptom-gsem-model_files/figure-gfm/mdd_symptom_gsem_efa_eigen-1.png)<!-- -->

Simulate uncertainty in **S** using the **V** matrix

``` r
# replicates
m <- 100

# simulate with S as the mean and V as the variance
S_lowertri_sim = mvtnorm::rmvnorm(m,
    mean=symptoms_covstruct$S[lower.tri(symptoms_covstruct$S, diag=TRUE)],
    sigma=symptoms_covstruct$V)
    
# reshape into a 24x24xm array
S_sim <- plyr::aaply(S_lowertri_sim, 1, function(x, k=nrow(symptoms_covstruct$S)){
    S <- matrix(NA, ncol=k, nrow=k)
    S[lower.tri(S, diag=T)] <- x
    S[upper.tri(S, diag=T)] <- t(S)[upper.tri(S, diag=T)]
    return(S)
})
dimnames(S_sim) <- list(1:m, colnames(symptoms_cov), colnames(symptoms_cov))

# find which symptoms have positive variances across all replicates
sim_cov_keep <- which(colSums(plyr::aaply(S_sim, 1, diag) > 0) == m)

eigen(cov2cor(symptoms_covstruct$S[sim_cov_keep,sim_cov_keep]))$values
```

    ##  [1]  4.341331957  1.815692348  1.250494118  1.036911153  0.676884625
    ##  [6]  0.585888811  0.329218719  0.050928011  0.008914931 -0.096264674

``` r
S_sim_pos <- S_sim[,sim_cov_keep,sim_cov_keep]

S_sim_pos_ev <- plyr::aaply(S_sim_pos, 1, function(x) eigen(cov2cor(x))$values)

summary(S_sim_pos_ev)
```

    ##        1                2               3               4         
    ##  Min.   : 3.930   Min.   :1.514   Min.   :1.082   Min.   :0.7671  
    ##  1st Qu.: 4.324   1st Qu.:1.896   1st Qu.:1.321   1st Qu.:0.9766  
    ##  Median : 4.536   Median :1.988   Median :1.445   Median :1.0579  
    ##  Mean   : 4.674   Mean   :2.063   Mean   :1.435   Mean   :1.0723  
    ##  3rd Qu.: 4.794   3rd Qu.:2.203   3rd Qu.:1.532   3rd Qu.:1.1670  
    ##  Max.   :14.311   Max.   :3.189   Max.   :1.854   Max.   :1.4730  
    ##        5                6                7                 8            
    ##  Min.   :0.5091   Min.   :0.2025   Min.   :0.01294   Min.   :-0.308324  
    ##  1st Qu.:0.6943   1st Qu.:0.4554   1st Qu.:0.15874   1st Qu.:-0.004917  
    ##  Median :0.7410   Median :0.5069   Median :0.26025   Median : 0.036113  
    ##  Mean   :0.7516   Mean   :0.5237   Mean   :0.25195   Mean   : 0.033933  
    ##  3rd Qu.:0.8133   3rd Qu.:0.5973   3rd Qu.:0.33319   3rd Qu.: 0.087804  
    ##  Max.   :1.0174   Max.   :0.7909   Max.   :0.56089   Max.   : 0.294800  
    ##        9                  10         
    ##  Min.   :-0.66013   Min.   :-9.4891  
    ##  1st Qu.:-0.21860   1st Qu.:-0.7115  
    ##  Median :-0.11874   Median :-0.5433  
    ##  Mean   :-0.14676   Mean   :-0.6581  
    ##  3rd Qu.:-0.04559   3rd Qu.:-0.3479  
    ##  Max.   : 0.05182   Max.   :-0.1194

``` r
symptoms_efa3 <- factanal(covmat=symptoms_cov_pd, factors=3, rotation='varimax')

symptoms_efa3
```

    ## 
    ## Call:
    ## factanal(factors = 3, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.895       0.438       0.545       0.271       0.950       0.226 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.114       0.016       0.567       0.494       0.581       0.673 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.005       0.005       0.609       0.600       0.527       0.661 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3
    ## ClinAppDec                  -0.323 
    ## ClinAppInc   0.140   0.152   0.720 
    ## ClinSleDec                   0.671 
    ## ClinSleInc   0.128   0.690   0.486 
    ## ClinPsycInc  0.163          -0.154 
    ## ClinSui      0.491   0.726         
    ## PopDep       0.903          -0.253 
    ## PopAnh       0.977  -0.161         
    ## PopAppDec    0.169   0.390  -0.503 
    ## PopAppInc    0.316  -0.630         
    ## PopSleDec    0.495           0.405 
    ## PopSleInc    0.537   0.165   0.106 
    ## PopPsycInc           0.993         
    ## PopPsycDec   0.178   0.637  -0.747 
    ## PopFatig     0.586           0.196 
    ## PopGuilt     0.596   0.211         
    ## PopConc      0.622   0.157   0.247 
    ## PopSui       0.526   0.107  -0.228 
    ## 
    ##                Factor1 Factor2 Factor3
    ## SS loadings      4.142   3.128   2.553
    ## Proportion Var   0.230   0.174   0.142
    ## Cumulative Var   0.230   0.404   0.546
    ## 
    ## The degrees of freedom for the model is 102 and the fit was 87.385

``` r
symptoms_efa4 <- factanal(covmat=symptoms_cov_pd, factors=4, rotation='varimax')

symptoms_efa4
```

    ## 
    ## Call:
    ## factanal(factors = 4, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.328       0.413       0.005       0.253       0.902       0.157 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.088       0.047       0.451       0.467       0.408       0.712 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.005       0.005       0.494       0.506       0.497       0.629 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3 Factor4
    ## ClinAppDec  -0.233           0.121   0.777 
    ## ClinAppInc   0.144           0.681  -0.310 
    ## ClinSleDec  -0.148           0.958   0.230 
    ## ClinSleInc   0.106   0.634   0.567  -0.108 
    ## ClinPsycInc  0.103                   0.297 
    ## ClinSui      0.586   0.707                 
    ## PopDep       0.829  -0.126           0.456 
    ## PopAnh       0.915  -0.215           0.253 
    ## PopAppDec            0.396  -0.173   0.596 
    ## PopAppInc    0.326  -0.649                 
    ## PopSleDec    0.408           0.621   0.198 
    ## PopSleInc    0.507   0.138                 
    ## PopPsycInc   0.137   0.978   0.138         
    ## PopPsycDec   0.194   0.678  -0.531   0.464 
    ## PopFatig     0.654  -0.126          -0.243 
    ## PopGuilt     0.677   0.182                 
    ## PopConc      0.655   0.106   0.225  -0.110 
    ## PopSui       0.548          -0.107   0.229 
    ## 
    ##                Factor1 Factor2 Factor3 Factor4
    ## SS loadings      4.196   3.055   2.516   1.867
    ## Proportion Var   0.233   0.170   0.140   0.104
    ## Cumulative Var   0.233   0.403   0.543   0.646
    ## 
    ## The degrees of freedom for the model is 87 and the fit was 84.9414

``` r
symptoms_efa6 <- factanal(covmat=symptoms_cov_pd, factors=6, rotation='varimax')

symptoms_efa6
```

    ## 
    ## Call:
    ## factanal(factors = 6, covmat = symptoms_cov_pd, rotation = "varimax")
    ## 
    ## Uniquenesses:
    ##  ClinAppDec  ClinAppInc  ClinSleDec  ClinSleInc ClinPsycInc     ClinSui 
    ##       0.005       0.176       0.005       0.005       0.864       0.005 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.090       0.005       0.388       0.473       0.246       0.177 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.052       0.005       0.163       0.408       0.357       0.574 
    ## 
    ## Loadings:
    ##             Factor1 Factor2 Factor3 Factor4 Factor5 Factor6
    ## ClinAppDec          -0.104                   0.986         
    ## ClinAppInc                   0.752          -0.493         
    ## ClinSleDec          -0.108   0.935           0.262   0.154 
    ## ClinSleInc   0.556  -0.183   0.470   0.213           0.620 
    ## ClinPsycInc          0.149           0.105   0.262  -0.165 
    ## ClinSui      0.779   0.304           0.532                 
    ## PopDep               0.888           0.262   0.226         
    ## PopAnh      -0.141   0.865           0.422   0.111   0.181 
    ## PopAppDec    0.496   0.305          -0.130   0.418  -0.277 
    ## PopAppInc   -0.516   0.340           0.205  -0.251  -0.176 
    ## PopSleDec            0.328   0.683   0.317   0.167  -0.224 
    ## PopSleInc            0.375           0.171           0.799 
    ## PopPsycInc   0.943           0.168                   0.155 
    ## PopPsycDec   0.776   0.313  -0.458           0.273         
    ## PopFatig    -0.101   0.201           0.872           0.159 
    ## PopGuilt     0.250   0.518           0.425  -0.176  -0.210 
    ## PopConc      0.112   0.290   0.185   0.710                 
    ## PopSui       0.225   0.585                  -0.129   0.109 
    ## 
    ##                Factor1 Factor2 Factor3 Factor4 Factor5 Factor6
    ## SS loadings      3.104   3.001   2.435   2.253   1.824   1.385
    ## Proportion Var   0.172   0.167   0.135   0.125   0.101   0.077
    ## Cumulative Var   0.172   0.339   0.474   0.600   0.701   0.778
    ## 
    ## The degrees of freedom for the model is 60 and the fit was 79.9357
