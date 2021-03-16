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
    ## platform       x86_64-conda-linux-gnu      
    ## arch           x86_64                      
    ## os             linux-gnu                   
    ## system         x86_64, linux-gnu           
    ## status                                     
    ## major          4                           
    ## minor          0.3                         
    ## year           2020                        
    ## month          10                          
    ## day            10                          
    ## svn rev        79318                       
    ## language       R                           
    ## version.string R version 4.0.3 (2020-10-10)
    ## nickname       Bunny-Wunnies Freak Out

Package installation

``` r
required_packages <- c('devtools', 'readr', 'tidyr', 'dplyr', 'ggplot2', 'stringr', 'corrplot')
for(pack in required_packages) if(!require(pack, character.only=TRUE)) install.packages(pack)

if(!require(GenomicSEM)) remotes::install_github("MichelNivard/GenomicSEM")

if(!require(tidySEM)) remotes::install_github("cjvanlissa/tidySEM")
```

GenomicSEM version

``` r
require(readr)
require(tidyr)
require(stringr)
require(dplyr)
require(ggplot2)
require(corrplot)
require(tidySEM)
require(GenomicSEM)

packageVersion("GenomicSEM")
```

    ## [1] '0.0.2'

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

    ## 
    ## ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────
    ## cols(
    ##   cohorts = col_character(),
    ##   symptom = col_character(),
    ##   sumstats = col_character(),
    ##   Nca = col_double(),
    ##   Nco = col_double(),
    ##   samp_prev = col_double(),
    ##   filename = col_character(),
    ##   pop_prev = col_double(),
    ##   trait_name = col_character()
    ## )

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

## ADGS-PGC

### Common factor

Common factor model. Allow residual negative correlation between
directional symptoms

``` r
pgc_commonfactor.model <- "
A1 =~ NA*ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
A1 ~~ 1*A1
c3a3b > -1
ClinAppDec ~~ c3a3b*ClinAppInc
"
pgc_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pgc_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   3.111 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0367078980491331 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  Inf . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pgc_commonfactor.model): A difference greater than .025 was observed pre-
    ## and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pgc_commonfactor.model): A difference greater than .025 was observed pre- and
    ## post-smoothing for Z-statistics in the genetic covariance matrix. This reflects
    ## a large difference and results should be interpreted with caution!! This can
    ## often result from including low powered traits, and you might consider removing
    ## those traits from the model. If you are going to run a multivariate GWAS we
    ## strongly recommend setting the smooth_check argument to true to check smoothing
    ## for each SNP.

``` r
pgc_commonfactor.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 2.604373  8 0.9566861 28.60437   1 0.1435089

``` r
pgc_commonfactor.fit$results[c(1,2,3,6,7)]
```

    ##            lhs op         rhs STD_Genotype   STD_Genotype_SE
    ## 1           A1 =~  ClinAppDec  0.671742966 0.433122026645965
    ## 2           A1 =~  ClinAppInc  0.729689563 0.521524441752345
    ## 3           A1 =~  ClinSleDec  0.746888149 0.474570264504223
    ## 4           A1 =~  ClinSleInc  0.612271642 0.461874360406521
    ## 5           A1 =~ ClinPsycInc  0.279997292 0.237173888205013
    ## 6           A1 =~     ClinSui  0.007296054 0.178669809982275
    ## 7           A1 ~~          A1  1.000000000                  
    ## 8   ClinAppDec ~~  ClinAppDec  0.548761368 0.574094487130573
    ## 9   ClinAppDec ~~  ClinAppInc -0.861169631 0.618301627013953
    ## 10  ClinAppInc ~~  ClinAppInc  0.467553163  0.79129051639139
    ## 11  ClinSleDec ~~  ClinSleDec  0.442157966 0.900065138452166
    ## 12  ClinSleInc ~~  ClinSleInc  0.625124988  1.61291792720481
    ## 13 ClinPsycInc ~~ ClinPsycInc  0.921601537 0.470933828537936
    ## 14     ClinSui ~~     ClinSui  0.999946802 0.363500519041763

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
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_commonfactor.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.083 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 61.30386 33 0.001975283 105.3039 0.9918045 0.1147905

``` r
pop_commonfactor.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.8607640  0.046690610179535
    ## 2         A1 =~    PopSui    0.5869717  0.098351794371004
    ## 3         A1 =~    PopAnh    0.9920183 0.0430337978866502
    ## 4         A1 =~ PopAppDec    0.1833142 0.0894139625344264
    ## 5         A1 =~ PopAppInc    0.4102110 0.0805622852724592
    ## 6         A1 =~ PopSleDec    0.5959608 0.0937668177610682
    ## 7         A1 =~ PopSleInc    0.5412317  0.102203925548959
    ## 8         A1 =~  PopFatig    0.6726953 0.0937812042741406
    ## 9         A1 =~  PopGuilt    0.6276465 0.0800933702538832
    ## 10        A1 =~   PopConc    0.7028475 0.0903103532806777
    ## 11        A1 ~~        A1    1.0000000                   
    ## 12    PopDep ~~    PopDep    0.2590854 0.0734986748136565
    ## 13    PopSui ~~    PopSui    0.6554639  0.253091999036569
    ## 14    PopAnh ~~    PopAnh    0.0158997 0.0672034527221335
    ## 15 PopAppDec ~~ PopAppDec    0.9663962  0.250347667091379
    ## 16 PopAppDec ~~ PopAppInc   -0.1514390  0.141684492913413
    ## 17 PopAppInc ~~ PopAppInc    0.8317261  0.166123963867933
    ## 18 PopSleDec ~~ PopSleDec    0.6448284  0.331827080112782
    ## 19 PopSleDec ~~ PopSleInc   -0.3359880  0.212587725438769
    ## 20 PopSleInc ~~ PopSleInc    0.7070673  0.282532187235941
    ## 21  PopFatig ~~  PopFatig    0.5474803  0.322349579364684
    ## 22  PopGuilt ~~  PopGuilt    0.6060609  0.165904819996586
    ## 23   PopConc ~~   PopConc    0.5060052  0.268317184474319

Remove common variance shared between the gating items (Mood:
`UKB_CIDI1`, Interest: `UKB_CIDI2`) that is uncorrelated with the common
factor variance, to recover the genetic structure among gated items

``` r
pop_commonfactor_gating.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_commonfactor_gating.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_commonfactor_gating.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.054 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 38.15576 32 0.2097698 84.15576 0.9982176 0.1036837

``` r
pop_commonfactor_gating.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7246509 0.0763743523874054
    ## 2         A1 =~    PopSui    0.6214609  0.105625167050975
    ## 3         A1 =~    PopAnh    0.8601181 0.0715887366883272
    ## 4         A1 =~ PopAppDec    0.2023905 0.0956710047479845
    ## 5         A1 =~ PopAppInc    0.4427121 0.0872178868578012
    ## 6         A1 =~ PopSleDec    0.6589985  0.103074320429907
    ## 7         A1 =~ PopSleInc    0.5802037   0.11056641098677
    ## 8         A1 =~  PopFatig    0.7307021  0.101382678954278
    ## 9         A1 =~  PopGuilt    0.6821160  0.088600662742439
    ## 10        A1 =~   PopConc    0.7723562 0.0982782672731918
    ## 11        A1 ~~        A1    1.0000000                   
    ## 12    PopDep ~~    PopDep    0.4748809  0.126069489691032
    ## 13    PopDep ~~    PopAnh    0.3005128  0.114419113125973
    ## 14    PopSui ~~    PopSui    0.6137866  0.258385211982711
    ## 15    PopAnh ~~    PopAnh    0.2601970  0.124998871722654
    ## 16 PopAppDec ~~ PopAppDec    0.9590382  0.249868147482632
    ## 17 PopAppDec ~~ PopAppInc   -0.1658429  0.141371760096243
    ## 18 PopAppInc ~~ PopAppInc    0.8040073  0.170690977987863
    ## 19 PopSleDec ~~ PopSleDec    0.5657224  0.335291834375094
    ## 20 PopSleDec ~~ PopSleInc   -0.3957891  0.211691692224068
    ## 21 PopSleInc ~~ PopSleInc    0.6633639  0.282985160089131
    ## 22  PopFatig ~~  PopFatig    0.4660745  0.324157816895692
    ## 23  PopGuilt ~~  PopGuilt    0.5347171  0.171268212578693
    ## 24   PopConc ~~   PopConc    0.4034663  0.270476790692872

### Kendler Neale model

``` r
pop_kendler_neale.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
PopAppDec ~~ PopAppInc
"
pop_kendler_neale.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_kendler_neale.model)
```

    ## [1] "Running primary model"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = pop_kendler_neale.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

    ##           lhs op       rhs Unstandardized_Estimate
    ## 1          A1 =~  PopGuilt           -1.587317e-01
    ## 2          A1 =~   PopConc           -1.687143e-01
    ## 3          A1 =~    PopSui           -1.074497e-01
    ## 4          A2 =~    PopDep            3.257335e+00
    ## 5          A2 =~    PopAnh            4.014837e+00
    ## 6          A2 =~  PopGuilt           -2.624842e-04
    ## 7          A3 =~ PopSleDec            1.320644e-01
    ## 8          A3 =~ PopSleInc            1.194367e-01
    ## 9          A3 =~  PopFatig            1.680700e-01
    ## 10         A3 =~ PopAppDec            3.698864e-02
    ## 11         A3 =~ PopAppInc            1.212150e-01
    ## 15     PopDep ~~    PopAnh           -1.300040e+01
    ## 16  PopSleDec ~~ PopSleInc           -1.635657e-02
    ## 17  PopAppDec ~~ PopAppInc           -8.481291e-03
    ## 133    PopDep ~~    PopDep           -1.052974e+01
    ## 134    PopAnh ~~    PopAnh           -1.603200e+01
    ## 135 PopAppDec ~~ PopAppDec            3.381163e-02
    ## 136 PopAppInc ~~ PopAppInc            6.345856e-02
    ## 137 PopSleDec ~~ PopSleDec            2.456743e-02
    ## 138 PopSleInc ~~ PopSleInc            3.060018e-02
    ## 139  PopFatig ~~  PopFatig            2.728230e-02
    ## 140  PopGuilt ~~  PopGuilt            3.214604e-02
    ## 141   PopConc ~~   PopConc            2.304104e-02
    ## 142    PopSui ~~    PopSui            2.077681e-02
    ## 195        A1 ~~        A2           -6.630526e-02
    ## 196        A1 ~~        A3           -1.096331e+00
    ## 197        A2 ~~        A3            6.439752e-02

Add constraints to prevent variances from being negative and
correlations from going out of bounds.

``` r
pop_kendler_neale_constr.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
PopDep ~~ PopAnh
PopSleDec ~~ PopSleInc
PopAppDec ~~ PopAppInc
c2 > 0.001
PopAnh ~~ c2*PopAnh
a13 > -1.0
A1 ~~ a13*A3
"
pop_kendler_neale_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_kendler_neale_constr.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    8.07 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_kendler_neale_constr.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_kendler_neale_constr.fit$modelfit
```

    ##     chisq df    p_chisq    AIC      CFI      SRMR
    ## df 38.485 28 0.08956883 92.485 0.996964 0.1034545

``` r
pop_kendler_neale_constr.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopSui   0.40812631  0.208103708831721
    ## 2         A1 =~  PopGuilt  -0.07345475  0.211342376306899
    ## 3         A1 =~   PopConc   0.50277825  0.251333320462741
    ## 4         A1 ~~        A1   1.00000000                   
    ## 5         A1 ~~        A2   1.53797968  0.735737883976811
    ## 6         A1 ~~        A3   1.71222908  0.784418260902913
    ## 7         A2 =~    PopDep   0.72962243 0.0844352429616578
    ## 8         A2 =~    PopAnh   0.86516791 0.0858289439411759
    ## 9         A2 =~  PopGuilt   0.78746635  0.296277503826445
    ## 10        A2 ~~        A2   1.00000000                   
    ## 11        A2 ~~        A3   1.01406796  0.160411495645851
    ## 12        A3 =~ PopAppDec   0.19657622 0.0946055813127982
    ## 13        A3 =~ PopAppInc   0.43542653 0.0956821673874205
    ## 14        A3 =~ PopSleDec   0.64316133  0.121428998511701
    ## 15        A3 =~ PopSleInc   0.56507097  0.117485895533465
    ## 16        A3 =~  PopFatig   0.71093275  0.123634929781074
    ## 17        A3 ~~        A3   1.00000000                   
    ## 18    PopDep ~~    PopDep   0.46765113  0.141645139230262
    ## 19    PopDep ~~    PopAnh   0.29255256  0.136434315852559
    ## 20    PopSui ~~    PopSui   0.83343046  0.294706962559891
    ## 21    PopAnh ~~    PopAnh   0.25148441  0.153844883182964
    ## 22 PopAppDec ~~ PopAppDec   0.96135682  0.249622147199745
    ## 23 PopAppDec ~~ PopAppInc  -0.16183684  0.140945795040004
    ## 24 PopAppInc ~~ PopAppInc   0.81040419  0.170094672488897
    ## 25 PopSleDec ~~ PopSleDec   0.58634096  0.329259854076718
    ## 26 PopSleDec ~~ PopSleInc  -0.37686382  0.217024607044167
    ## 27 PopSleInc ~~ PopSleInc   0.68069864  0.287671473332995
    ## 28  PopFatig ~~  PopFatig   0.49457474   0.32138896180926
    ## 29  PopGuilt ~~  PopGuilt   0.55242426  0.182793159275404
    ## 30   PopConc ~~   PopConc   0.74721633  0.328998385787683

Check that an orthogonal model has poor fit.

``` r
pop_kendler_neale_orth.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
PopSleDec ~~ PopSleInc
PopAppDec ~~ PopAppInc
A1 ~~ 0*A2
A1 ~~ 0*A3
A2 ~~ 0*A3
a7 > 0.001
PopGuilt ~~ a7*PopGuilt
"
pop_kendler_neale_orth.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_kendler_neale_orth.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"

    ## Warning in sqrt(1/diag(V)): NaNs produced

    ## Warning in cov2cor(Sigma.hat): diag(.) had 0 or NA entries; non-finite result is
    ## doubtful

    ## Warning in sqrt(1/diag(V)): NaNs produced

    ## Warning in cov2cor(Sigma.hat): diag(.) had 0 or NA entries; non-finite result is
    ## doubtful

    ## [1] "Calculating SRMR"

    ## Warning in sqrt(1/diag(V)): NaNs produced

    ## Warning in sqrt(1/diag(V)): diag(.) had 0 or NA entries; non-finite result is
    ## doubtful

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_kendler_neale_orth.model): CFI estimates below 0 should not be trusted, and
    ## indicate that the other model fit estimates should be interpreted with caution.
    ## A negative CFI estimates typically appears due to negative residual variances.

    ## elapsed 
    ##   3.288 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_kendler_neale_orth.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_kendler_neale_orth.fit$modelfit
```

    ##     chisq df p_chisq    AIC       CFI     SRMR
    ## df 593087 32       0 593133 -170.7224 2.717298

``` r
pop_kendler_neale_orth.fit$results[c(1, 2, 3, 6, 7)]
```

    ##          lhs op       rhs  STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopSui  0.5490753478  0.219388894665628
    ## 2         A1 =~  PopGuilt  0.8515294724   0.33849587231824
    ## 3         A1 =~   PopConc  0.5752659870  0.225884005549328
    ## 4         A1 ~~        A1  1.0000000000                   
    ## 5         A2 =~    PopDep  0.9454733422 0.0634537881417689
    ## 6         A2 =~    PopAnh  0.9787822366 0.0691237124410743
    ## 7         A2 =~  PopGuilt  0.5556669821 0.0862445443303973
    ## 8         A2 ~~        A2  1.0000000000                   
    ## 9         A3 =~ PopAppDec  0.2916993554  0.193645914657204
    ## 10        A3 =~ PopAppInc  0.3410010371  0.173828358667647
    ## 11        A3 =~ PopSleDec  0.7871353633   0.36316786621252
    ## 12        A3 =~ PopSleInc  0.6901433198  0.328511409149719
    ## 13        A3 =~  PopFatig  0.6225957497  0.199025839974674
    ## 14        A3 ~~        A3  1.0000000000                   
    ## 15    PopDep ~~    PopDep  0.1060804340   0.11953858305352
    ## 16    PopSui ~~    PopSui  0.6985137520  0.330027700708191
    ## 17    PopAnh ~~    PopAnh  0.0419846771   0.12247123763652
    ## 18 PopAppDec ~~ PopAppDec  0.9149116198  0.266631203398247
    ## 19 PopAppDec ~~ PopAppInc -0.1757122303  0.155248197613913
    ## 20 PopAppInc ~~ PopAppInc  0.8837180365  0.179986298570089
    ## 21 PopSleDec ~~ PopSleDec  0.3804160596  0.626904612353322
    ## 22 PopSleDec ~~ PopSleInc -0.5566712621  0.428187375740196
    ## 23 PopSleInc ~~ PopSleInc  0.5237013893  0.480650914327075
    ## 24  PopFatig ~~  PopFatig  0.6123748767  0.366599657088121
    ## 25  PopGuilt ~~  PopGuilt  0.0009991597  0.602676440728349
    ## 26   PopConc ~~   PopConc  0.6690423736  0.351840815177709

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
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_psych_soma.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma.model)
```

    ## [1] "Running primary model"
    ## [1] "Error: The primary model produced correlations among your latent variables that are either greater than 1 or less than -1, or the latent variables have negative variances. \n              Consequently, model fit estimates could not be computed and results should likely not be interpreted. Results are provided below \n              to enable troubleshooting. A model constraint that constrains the latent correlations to be above -1, less than 1, or to have positive variances is suggested."
    ##           lhs op       rhs Unstandardized_Estimate          SE
    ## 1          A1 =~    PopDep             0.204328193 0.021831039
    ## 2          A1 =~    PopAnh             0.252141713 0.021811659
    ## 3          A1 =~  PopGuilt             0.163417009 0.021216029
    ## 4          A1 =~   PopConc             0.175667925 0.022290435
    ## 5          A1 =~    PopSui             0.111678516 0.018967663
    ## 6          A2 =~ PopAppDec             0.036939510 0.017677739
    ## 7          A2 =~ PopAppInc             0.121231693 0.026664591
    ## 8          A2 =~ PopSleDec             0.131767596 0.024798301
    ## 9          A2 =~ PopSleInc             0.119751996 0.024859277
    ## 10         A2 =~  PopFatig             0.168092440 0.029266122
    ## 13     PopDep ~~    PopAnh             0.025752494 0.009733175
    ## 14  PopAppDec ~~ PopAppInc            -0.008475956 0.007388728
    ## 15  PopSleDec ~~ PopSleInc            -0.016362657 0.009419571
    ## 111    PopDep ~~    PopDep             0.038747143 0.010209979
    ## 112    PopAnh ~~    PopAnh             0.023342741 0.011283331
    ## 113 PopAppDec ~~ PopAppDec             0.033815261 0.008783860
    ## 114 PopAppInc ~~ PopAppInc             0.063454492 0.013283958
    ## 115 PopSleDec ~~ PopSleDec             0.024645708 0.013856102
    ## 116 PopSleInc ~~ PopSleInc             0.030524774 0.012911799
    ## 117  PopFatig ~~  PopFatig             0.027274744 0.017846624
    ## 118  PopGuilt ~~  PopGuilt             0.030631241 0.009826456
    ## 119   PopConc ~~   PopConc             0.020646342 0.013894571
    ## 120    PopSui ~~    PopSui             0.019850151 0.008351031
    ## 173        A1 ~~        A2             1.035710365 0.147300968

``` r
pop_psych_soma_constr.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui 
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
a1 > 0.001
PopDep ~~ a1*PopDep
a2 > 0.001
PopAnh ~~ a2*PopAnh
a12 < 1
A1 ~~ a12*A2
"
pop_psych_soma_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma_constr.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   7.945 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC      CFI      SRMR
    ## df 38.15578 31 0.1761158 86.15578 0.997928 0.1036837

``` r
pop_psych_soma_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7246508 0.0775829496084442
    ## 2         A1 =~    PopSui    0.6214608  0.105608298282961
    ## 3         A1 =~    PopAnh    0.8601181 0.0746335392077814
    ## 4         A1 =~  PopGuilt    0.6821157 0.0887167291747381
    ## 5         A1 =~   PopConc    0.7723560 0.0983376230720235
    ## 6         A1 ~~        A1    1.0000000                   
    ## 7         A1 ~~        A2    0.9999999  0.138107401798997
    ## 8         A2 =~ PopAppDec    0.2023907  0.096242804229106
    ## 9         A2 =~ PopAppInc    0.4427122 0.0957164237154349
    ## 10        A2 =~ PopSleDec    0.6589995  0.121362037193989
    ## 11        A2 =~ PopSleInc    0.5802043   0.11839857969779
    ## 12        A2 =~  PopFatig    0.7307031  0.123576947330043
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopDep ~~    PopDep    0.4748811  0.128159469175355
    ## 15    PopDep ~~    PopAnh    0.3005131  0.117830103992584
    ## 16    PopSui ~~    PopSui    0.6137865  0.258441632016061
    ## 17    PopAnh ~~    PopAnh    0.2601968  0.131502440495488
    ## 18 PopAppDec ~~ PopAppDec    0.9590375  0.249638493891708
    ## 19 PopAppDec ~~ PopAppInc   -0.1658432  0.141243260593445
    ## 20 PopAppInc ~~ PopAppInc    0.8040054  0.170817214322806
    ## 21 PopSleDec ~~ PopSleDec    0.5657113   0.33024282168413
    ## 22 PopSleDec ~~ PopSleInc   -0.3957915   0.21733657576381
    ## 23 PopSleInc ~~ PopSleInc    0.6633571  0.288309709940671
    ## 24  PopFatig ~~  PopFatig    0.4660634  0.321972813038065
    ## 25  PopGuilt ~~  PopGuilt    0.5347181  0.171437353757826
    ## 26   PopConc ~~   PopConc    0.4034674  0.269866821744325

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
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_psych_soma_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_soma_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.148 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 17.38451 22 0.7417938 83.38451   1 0.05908159

``` r
pop_psych_soma_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1          A =~    PopDep   0.64647220 0.0885699318417481
    ## 2          A =~    PopSui   0.47828533  0.143138439990277
    ## 3          A =~    PopAnh   0.81989287  0.071452753220714
    ## 4          A =~ PopAppDec   0.21798394  0.100595273573109
    ## 5          A =~ PopAppInc   0.46941788 0.0909783955451097
    ## 6          A =~ PopSleDec   0.69695990  0.112597907130706
    ## 7          A =~ PopSleInc   0.59917069  0.114396722042244
    ## 8          A =~  PopFatig   0.77865776  0.108844820844914
    ## 9          A =~  PopGuilt   0.63402541  0.102884635489159
    ## 10         A =~   PopConc   0.90418226  0.118287997656812
    ## 11         A ~~         A   1.00000000                   
    ## 12        A1 =~    PopDep   0.40850580  0.205112931359441
    ## 13        A1 =~    PopSui   0.76535147  0.371260928089008
    ## 14        A1 =~    PopAnh   0.15280906  0.170476204366291
    ## 15        A1 =~  PopGuilt   0.27230502  0.182028242166464
    ## 16        A1 =~   PopConc  -0.27794994  0.222349961666451
    ## 17        A1 ~~        A1   1.00000000                   
    ## 18        A2 =~ PopAppDec   0.19911726  0.232546545338366
    ## 19        A2 =~ PopAppInc   0.23424750  0.252266754176024
    ## 20        A2 =~ PopSleDec   0.62313528  0.648635392187037
    ## 21        A2 =~ PopSleInc  -0.58716132  0.503226587953096
    ## 22        A2 =~  PopFatig  -0.38693898  0.337161364171656
    ## 23        A2 ~~        A2   1.00000000                   
    ## 24    PopDep ~~    PopDep   0.41519659  0.181886005226791
    ## 25    PopDep ~~    PopAnh   0.33133678   0.12249527326515
    ## 26    PopSui ~~    PopSui   0.18548010  0.577268307841389
    ## 27    PopAnh ~~    PopAnh   0.30442505  0.117310956838274
    ## 28 PopAppDec ~~ PopAppDec   0.91283565  0.250996162075526
    ## 29 PopAppDec ~~ PopAppInc  -0.22521076  0.160902574249922
    ## 30 PopAppInc ~~ PopAppInc   0.72477466  0.224953301363296
    ## 31 PopSleDec ~~ PopSleDec   0.12594949  0.835339406614409
    ## 32 PopSleDec ~~ PopSleInc  -0.06515137  0.658895615847185
    ## 33 PopSleInc ~~ PopSleInc   0.29623589  0.664725748894597
    ## 34  PopFatig ~~  PopFatig   0.24397109  0.445917359074471
    ## 35  PopGuilt ~~  PopGuilt   0.52386195  0.177047610220656
    ## 36   PopConc ~~   PopConc   0.10519833  0.323401932682589

### Psychological-Neurovegetative (Elhai Model 2b)

``` r
pop_psych_veg.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_psych_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.563 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##       chisq df   p_chisq      AIC       CFI       SRMR
    ## df 39.37797 31 0.1436769 87.37797 0.9975741 0.09756749

``` r
pop_psych_veg.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7657804 0.0876349999792341
    ## 2         A1 =~    PopSui    0.6209478  0.106202737847909
    ## 3         A1 =~    PopAnh    0.9055401 0.0901606452800706
    ## 4         A1 =~  PopGuilt    0.6742687 0.0893563522759828
    ## 5         A1 ~~        A1    1.0000000                   
    ## 6         A1 ~~        A2    0.8709567  0.101798824132556
    ## 7         A2 =~ PopAppDec    0.2145114  0.101006509988333
    ## 8         A2 =~ PopAppInc    0.4682453  0.094106504567065
    ## 9         A2 =~ PopSleDec    0.7039474   0.11389932088019
    ## 10        A2 =~ PopSleInc    0.6166593  0.117670035082161
    ## 11        A2 =~  PopFatig    0.7737486  0.112219106210762
    ## 12        A2 =~   PopConc    0.8255673  0.107067934472661
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopDep ~~    PopDep    0.4135806  0.150162754061393
    ## 15    PopDep ~~    PopAnh    0.2303534  0.147032229678117
    ## 16    PopSui ~~    PopSui    0.6144242  0.258908717850755
    ## 17    PopAnh ~~    PopAnh    0.1799972  0.168008747298596
    ## 18 PopAppDec ~~ PopAppDec    0.9539851  0.249400112645498
    ## 19 PopAppDec ~~ PopAppInc   -0.1766863  0.142630825343451
    ## 20 PopAppInc ~~ PopAppInc    0.7807458   0.17357760463102
    ## 21 PopSleDec ~~ PopSleDec    0.5044582  0.334476985708858
    ## 22 PopSleDec ~~ PopSleInc   -0.4475303  0.213922104732558
    ## 23 PopSleInc ~~ PopSleInc    0.6197311  0.285097972386903
    ## 24  PopFatig ~~  PopFatig    0.4013152  0.320345293468049
    ## 25  PopGuilt ~~  PopGuilt    0.5453609  0.172076937939494
    ## 26   PopConc ~~   PopConc    0.3184394  0.277016654372463

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
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_psych_veg_bif.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##    1.13 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 14.94934 22 0.8644029 80.94934   1 0.06129222

``` r
pop_psych_veg_bif.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1          A =~    PopDep    0.6051552 0.0917451014321342
    ## 2          A =~    PopSui    0.4232915  0.144262921645641
    ## 3          A =~    PopAnh    0.8114164  0.082159094410916
    ## 4          A =~ PopAppDec    0.2181404  0.102564290342008
    ## 5          A =~ PopAppInc    0.4727039 0.0960798646943193
    ## 6          A =~ PopSleDec    0.6949500  0.119999703987343
    ## 7          A =~ PopSleInc    0.6388870  0.120650042479437
    ## 8          A =~  PopFatig    0.8224599  0.125396141747029
    ## 9          A =~  PopGuilt    0.6039098  0.105840249846791
    ## 10         A =~   PopConc    0.8389293  0.105747165709026
    ## 11         A ~~         A    1.0000000                   
    ## 12        A1 =~    PopDep    0.4767507  0.222613046661377
    ## 13        A1 =~    PopSui    0.7687473  0.368381433450844
    ## 14        A1 =~    PopAnh    0.2001627  0.183960854991122
    ## 15        A1 =~  PopGuilt    0.3384348  0.178284865148621
    ## 16        A1 ~~        A1    1.0000000                   
    ## 17        A2 =~ PopAppDec    0.1803789  0.224797307526843
    ## 18        A2 =~ PopAppInc    0.3059287  0.251207708177112
    ## 19        A2 =~ PopSleDec    0.5852114  0.465534048899667
    ## 20        A2 =~ PopSleInc   -0.4919360  0.350045328254747
    ## 21        A2 =~  PopFatig   -0.4026424  0.336093613490183
    ## 22        A2 =~   PopConc    0.2227521  0.195589764735015
    ## 23        A2 ~~        A2    1.0000000                   
    ## 24    PopDep ~~    PopDep    0.4064961  0.201842263297389
    ## 25    PopDep ~~    PopAnh    0.3373382    0.1231832354606
    ## 26    PopSui ~~    PopSui    0.2298516  0.567590000215953
    ## 27    PopAnh ~~    PopAnh    0.3015381  0.118758932620607
    ## 28 PopAppDec ~~ PopAppDec    0.9198776   0.24607797270643
    ## 29 PopAppDec ~~ PopAppInc   -0.2345410   0.15803323295662
    ## 30 PopAppInc ~~ PopAppInc    0.6829582  0.243979376974387
    ## 31 PopSleDec ~~ PopSleDec    0.1745724  0.579307447121603
    ## 32 PopSleDec ~~ PopSleInc   -0.1695422  0.434108876939409
    ## 33 PopSleInc ~~ PopSleInc    0.3498221  0.463495046440098
    ## 34  PopFatig ~~  PopFatig    0.1614394  0.475109222849474
    ## 35  PopGuilt ~~  PopGuilt    0.5207548  0.180035639392019
    ## 36   PopConc ~~   PopConc    0.2465792  0.296262976462272

### Affective-Neurovegetative (Elhai Model 2c)

``` r
pop_affect_veg.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
PopAppInc ~~ PopAppDec
PopSleInc ~~ PopSleDec
"
pop_affect_veg.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   1.209 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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

    ##                                              chisq df  AIC
    ## df The follow-up chi-square model did not converge 31 <NA>
    ##                                                                         CFI
    ## df Either the chi-square or null (i.e. independence) model did not converge
    ##                  SRMR
    ## df 0.0955728303687932

``` r
pop_affect_veg.fit$results
```

    ##          lhs op       rhs  Unstand_Est          Unstand_SE STD_Genotype
    ## 1         A1 =~    PopDep  0.219637297  0.0227211312227577    0.7741338
    ## 2         A1 =~    PopSui  0.123781881  0.0209528281945554    0.6885031
    ## 3         A1 =~  PopGuilt  0.177917525  0.0236584570286192    0.7430254
    ## 4         A1 ~~        A1  1.000000000                        1.0000000
    ## 5         A1 ~~        A2  0.818684317  0.0647665731437532    0.8186843
    ## 6         A2 =~    PopAnh  0.254704694  0.0209575343600385    0.8639359
    ## 7         A2 =~ PopAppDec  0.039192137  0.0186157609018021    0.2089550
    ## 8         A2 =~ PopAppInc  0.128509748  0.0250747050429693    0.4596915
    ## 9         A2 =~ PopSleDec  0.139632298  0.0219408618724336    0.6812670
    ## 10        A2 =~ PopSleInc  0.127907205  0.0239537116787179    0.6038650
    ## 11        A2 =~  PopFatig  0.179135764  0.0245324230188376    0.7601852
    ## 12        A2 =~   PopConc  0.182461053  0.0227945837234133    0.8039762
    ## 13        A2 ~~        A2  1.000000000                        1.0000000
    ## 14    PopDep ~~    PopDep  0.032256553  0.0112806352554952    0.4007169
    ## 15    PopDep ~~    PopAnh  0.031472736 0.00865613809502582    0.3762610
    ## 16    PopSui ~~    PopSui  0.017000283 0.00860593715052202    0.5259653
    ## 17    PopAnh ~~    PopAnh  0.022043697  0.0108421522408869    0.2536148
    ## 18 PopAppDec ~~ PopAppDec  0.033643754 0.00878210168237376    0.9563391
    ## 19 PopAppInc ~~ PopAppDec -0.009034245 0.00745100880265082   -0.1722972
    ## 20 PopAppInc ~~ PopAppInc  0.061636872  0.0134569190935185    0.7886832
    ## 21 PopSleDec ~~ PopSleDec  0.022511265  0.0140228528161719    0.5358757
    ## 22 PopSleInc ~~ PopSleDec -0.018443217 0.00921507839072616   -0.4248272
    ## 23 PopSleInc ~~ PopSleInc  0.028505018  0.0127151904655436    0.6353453
    ## 24  PopFatig ~~  PopFatig  0.023440203  0.0178172615195873    0.4221185
    ## 25  PopGuilt ~~  PopGuilt  0.025681685   0.010561418809187    0.4479133
    ## 26   PopConc ~~   PopConc  0.018213506  0.0138497560201302    0.3536241
    ##       STD_Genotype_SE    STD_All      p_value
    ## 1  0.0800829390913145  0.7741338 4.178183e-22
    ## 2   0.116544539796899  0.6885025 3.470318e-09
    ## 3  0.0988033090404851  0.7430254 5.467137e-14
    ## 4                      1.0000000           NA
    ## 5  0.0647665877959047  0.8186843 1.261878e-36
    ## 6  0.0710861354122772  0.8639359 5.504801e-34
    ## 7  0.0992508361926856  0.2089549 3.526345e-02
    ## 8  0.0896947053664451  0.4596916 2.974191e-07
    ## 9   0.107049697775534  0.6812668 1.965280e-10
    ## 10  0.113088263952322  0.6038655 9.306684e-08
    ## 11  0.104106382993355  0.7601852 2.835203e-13
    ## 12  0.100439525330741  0.8039755 1.198758e-15
    ## 13                     1.0000000           NA
    ## 14  0.140137219156871  0.4007168 4.243597e-03
    ## 15  0.103485332505363  1.1802747 2.770263e-04
    ## 16  0.266254294936349  0.5259643 4.822126e-02
    ## 17  0.124739730246954  0.2536148 4.203760e-02
    ## 18  0.249634837469471  0.9563379 1.276517e-04
    ## 19  0.142101723227113 -0.1983905 2.253263e-01
    ## 20  0.172189801294967  0.7886836 4.642855e-06
    ## 21  0.333810347541146  0.5358755 1.084218e-01
    ## 22  0.212263535278602 -0.7280738 4.534745e-02
    ## 23   0.28340806891535  0.6353464 2.497377e-02
    ## 24   0.32085926423244  0.4221185 1.883118e-01
    ## 25  0.184201131830184  0.4479133 1.503018e-02
    ## 26  0.268898259806467  0.3536235 1.884838e-01

``` r
pop_affect_veg_constr.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopDep ~~ PopAnh
c2 > 0.001
PopAnh ~~ c2*PopAnh
c7 > 0.001
PopGuilt ~~ c7*PopGuilt
"
pop_affect_veg_constr.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_affect_veg_constr.model, fix_resid=FALSE)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   7.536 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_constr.model, : A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_affect_veg_constr.fit$modelfit
```

    ##       chisq df   p_chisq      AIC CFI      SRMR
    ## df 30.35264 33 0.5995866 74.35264   1 0.1123385

``` r
pop_affect_veg_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1         A1 =~    PopDep    0.7816760 0.0817349861503451
    ## 2         A1 =~    PopSui    0.6837610  0.116024414698104
    ## 3         A1 =~  PopGuilt    0.7371764 0.0983352213533668
    ## 4         A1 ~~        A1    1.0000000                   
    ## 5         A1 ~~        A2    0.8289461 0.0648763361954196
    ## 6         A2 =~    PopAnh    0.8801193 0.0732211968083118
    ## 7         A2 =~ PopAppDec    0.1885544 0.0966207299598514
    ## 8         A2 =~ PopAppInc    0.4522553 0.0902475988336744
    ## 9         A2 =~ PopSleDec    0.6438605  0.107671199224865
    ## 10        A2 =~ PopSleInc    0.5659901  0.111444757773795
    ## 11        A2 =~  PopFatig    0.7562465  0.104205048640434
    ## 12        A2 =~   PopConc    0.7982988  0.100661113105576
    ## 13        A2 ~~        A2    1.0000000                   
    ## 14    PopDep ~~    PopDep    0.3889824  0.143675618942531
    ## 15    PopDep ~~    PopAnh    0.3535102  0.108409767390845
    ## 16    PopSui ~~    PopSui    0.5324699  0.265339093640685
    ## 17    PopAnh ~~    PopAnh    0.2253902  0.130825392449184
    ## 18 PopAppDec ~~ PopAppDec    0.9644479   0.25131757017448
    ## 19 PopAppInc ~~ PopAppInc    0.7954641    0.1732316899317
    ## 20 PopSleDec ~~ PopSleDec    0.5854457   0.33241074917095
    ## 21 PopSleInc ~~ PopSleInc    0.6796531  0.285924264376715
    ## 22  PopFatig ~~  PopFatig    0.4280916  0.321064659002741
    ## 23  PopGuilt ~~  PopGuilt    0.4565706  0.182906241960823
    ## 24   PopConc ~~   PopConc    0.3627199  0.269143537182512

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
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_affect_veg_bif.model): CFI estimates below 0 should not be trusted, and
    ## indicate that the other model fit estimates should be interpreted with caution.
    ## A negative CFI estimates typically appears due to negative residual variances.

    ## elapsed 
    ##   1.102 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##   6.127 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.00241743805150268 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  0.215273668975602 . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

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
    ## df 17.14824 24 0.8422709 79.14824   1 0.06815131

``` r
pop_affect_veg_bif_constr.fit$results[c(1,2,3,6,7)]
```

    ##          lhs op       rhs STD_Genotype    STD_Genotype_SE
    ## 1          A =~    PopDep  0.627332287 0.0863700007367264
    ## 2          A =~    PopSui  0.510751966  0.108481053294169
    ## 3          A =~    PopAnh  0.856559079 0.0813452432490731
    ## 4          A =~ PopAppDec  0.203986003 0.0994409801322619
    ## 5          A =~ PopAppInc  0.467005284 0.0914896337698255
    ## 6          A =~ PopSleDec  0.764568865  0.125881594340571
    ## 7          A =~ PopSleInc  0.534946357  0.127271456512818
    ## 8          A =~  PopFatig  0.738520614  0.112053975395942
    ## 9          A =~  PopGuilt  0.641209445 0.0944833041121248
    ## 10         A =~   PopConc  0.823003795  0.105970218528199
    ## 11         A ~~         A  1.000000000                   
    ## 12        A1 =~    PopDep  0.451178030  0.265063627659287
    ## 13        A1 =~    PopSui  0.663759684   0.36788525267666
    ## 14        A1 =~  PopGuilt  0.283690739  0.165843630971553
    ## 15        A1 ~~        A1  1.000000000                   
    ## 16        A2 =~    PopAnh  0.210043179  0.111698031732063
    ## 17        A2 =~ PopAppDec -0.122773299  0.191515956721645
    ## 18        A2 =~ PopAppInc -0.112424095  0.175249331694881
    ## 19        A2 =~ PopSleDec -0.692900191  0.348172500388403
    ## 20        A2 =~ PopSleInc  0.642253672  0.256241685091402
    ## 21        A2 =~  PopFatig  0.378446165  0.224545052708205
    ## 22        A2 =~   PopConc -0.083851048  0.190901268363022
    ## 23        A2 ~~        A2  1.000000000                   
    ## 24    PopDep ~~    PopDep  0.402892601  0.227471265316279
    ## 25    PopDep ~~    PopAnh  0.386451141  0.111718449484741
    ## 26    PopSui ~~    PopSui  0.298555201  0.514634431948263
    ## 27    PopAnh ~~    PopAnh  0.222188339  0.140480896247316
    ## 28 PopAppDec ~~ PopAppDec  0.943316728  0.250727900584799
    ## 29 PopAppInc ~~ PopAppInc  0.769266355   0.19003535685593
    ## 30 PopSleDec ~~ PopSleDec  0.001000239  0.597977430761541
    ## 31 PopSleInc ~~ PopSleInc  0.301339046  0.387421642262966
    ## 32  PopFatig ~~  PopFatig  0.311365930  0.340687522675193
    ## 33  PopGuilt ~~  PopGuilt  0.508370094    0.1788505330553
    ## 34   PopConc ~~   PopConc  0.315633448  0.276353067832643

### Model comparisons

``` r
data.frame(Model=c('1a', '1b', '2a', '2a(ii)', '2b(i)', '2b(ii)', '2c(i)', '2c(ii)', '3'),
       Name=c('Common',
              'Common (gating)',
              'Psych-Somatic',
              'Psych-Somatic (BiF)',
              'Psych-Neuroveg',
              'Psych-Neuroveg (BiF)',
              'Affect-Neuroveg',
              'Affect-Neuroveg (BiF)',
              'Cog-Mood-Neuroveg')) %>%
bind_cols(
bind_rows(
lapply(list(pop_commonfactor.fit,
            pop_commonfactor_gating.fit,
            pop_psych_soma_constr.fit,
            pop_psych_soma_bif.fit,
            pop_psych_veg.fit,
            pop_psych_veg_bif.fit,
            pop_affect_veg_constr.fit,
            pop_affect_veg_bif_constr.fit,
            pop_kendler_neale_constr.fit),
       function(fit) signif(bind_cols(fit$modelfit, fit$results %>%
                     filter(lhs == 'A1' & rhs == 'A2') %>%
                     summarize(A1A2cor=mean(STD_Genotype), A1A2se=mean(as.numeric(STD_Genotype_SE)))), 3))
))
```

    ##         Model                  Name chisq df p_chisq   AIC   CFI   SRMR A1A2cor
    ## df...1     1a                Common  61.3 33 0.00198 105.0 0.992 0.1150     NaN
    ## df...2     1b       Common (gating)  38.2 32 0.21000  84.2 0.998 0.1040     NaN
    ## df...3     2a         Psych-Somatic  38.2 31 0.17600  86.2 0.998 0.1040   1.000
    ## df...4 2a(ii)   Psych-Somatic (BiF)  17.4 22 0.74200  83.4 1.000 0.0591     NaN
    ## df...5  2b(i)        Psych-Neuroveg  39.4 31 0.14400  87.4 0.998 0.0976   0.871
    ## df...6 2b(ii)  Psych-Neuroveg (BiF)  14.9 22 0.86400  80.9 1.000 0.0613     NaN
    ## df...7  2c(i)       Affect-Neuroveg  30.4 33 0.60000  74.4 1.000 0.1120   0.829
    ## df...8 2c(ii) Affect-Neuroveg (BiF)  17.1 24 0.84200  79.1 1.000 0.0682     NaN
    ## df...9      3     Cog-Mood-Neuroveg  38.5 28 0.08960  92.5 0.997 0.1030   1.540
    ##        A1A2se
    ## df...1    NaN
    ## df...2    NaN
    ## df...3 0.1380
    ## df...4    NaN
    ## df...5 0.1020
    ## df...6    NaN
    ## df...7 0.0649
    ## df...8    NaN
    ## df...9 0.7360

# Regression models

Regress clinical symptoms on the population symptom factors

``` r
pop_psych_veg_bif_clin_reg_a.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
A ~ ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c1 > 0.001
PopDep ~~ c1*PopDep
"
pop_psych_veg_bif_clin_reg_a.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif_clin_reg_a.model)
```

    ## [1] "Running primary model"
    ## [1] "Calculating model chi-square"
    ## [1] "Calculating CFI"
    ## [1] "Calculating Standardized Results"
    ## [1] "Calculating SRMR"
    ## elapsed 
    ##  31.765 
    ## [1] "The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was  0.0371957759493139 As a result of the smoothing, the largest Z-statistic change for the genetic covariances was  Inf . We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS."

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg_bif_clin_reg_a.model): A difference greater than .025 was observed
    ## pre- and post-smoothing in the genetic covariance matrix. This reflects a large
    ## difference and results should be interpreted with caution!! This can often
    ## result from including low powered traits, and you might consider removing those
    ## traits from the model. If you are going to run a multivariate GWAS we strongly
    ## recommend setting the smooth_check argument to true to check smoothing for each
    ## SNP.

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model =
    ## pop_psych_veg_bif_clin_reg_a.model): A difference greater than .025 was observed
    ## pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This
    ## reflects a large difference and results should be interpreted with caution!!
    ## This can often result from including low powered traits, and you might consider
    ## removing those traits from the model. If you are going to run a multivariate
    ## GWAS we strongly recommend setting the smooth_check argument to true to check
    ## smoothing for each SNP.

``` r
pop_psych_veg_bif_clin_reg_a.fit$modelfit
```

    ##       chisq df      p_chisq      AIC       CFI      SRMR
    ## df 325.4651 92 8.014794e-28 413.4651 0.9739853 0.1697943

``` r
pop_psych_veg_bif_clin_reg_a.fit$results[c(1:3,6,7,9)] %>%
filter(op == '~')
```

    ##   lhs op         rhs STD_Genotype   STD_Genotype_SE   p_value
    ## 1   A  ~  ClinAppDec   -0.1572080 0.220790433662912 0.4765732
    ## 2   A  ~  ClinAppInc    0.4392811  0.40460264386533 0.2776314
    ## 3   A  ~  ClinSleDec    0.1257140 0.260378340917268 0.6292028
    ## 4   A  ~  ClinSleInc    0.3921574 0.576628130467399 0.4964694
    ## 5   A  ~ ClinPsycInc    0.4384564   0.4219321254146 0.2987260
    ## 6   A  ~     ClinSui    1.3947971  1.22376987172374 0.2543979

``` r
pop_psych_veg_bif_clin_reg_a1.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
A1 ~ ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c1 > 0.001
PopDep ~~ c1*PopDep
"
pop_psych_veg_bif_clin_reg_a1.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif_clin_reg_a1.model)
```

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = pop_psych_veg_bif_clin_reg_a1.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

``` r
pop_psych_veg_bif_clin_reg_a2.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A  =~ NA*PopDep + PopAnh + PopGuilt + PopSui + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A  ~~ 1*A
A ~~ 0*A1
A ~~ 0*A2
A1 ~~ 0*A2
A2 ~ ClinAppDec + ClinAppInc + ClinSleDec + ClinSleInc + ClinPsycInc + ClinSui
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
c1 > 0.001
PopDep ~~ c1*PopDep
"
pop_psych_veg_bif_clin_reg_a2.fit <- usermodel(symptoms_covstruct, estimation='DWLS', model=pop_psych_veg_bif_clin_reg_a2.model)
```

    ## Warning in usermodel(symptoms_covstruct, estimation = "DWLS", model = pop_psych_veg_bif_clin_reg_a2.model): Error: The primary model did not converge! Additional warnings or errors are likely being printed by lavaan. 
    ##             The model output is also printed below (without standard errors) in case this is helpful for troubleshooting. Please note
    ##             that these results should not be interpreted.

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

    ##  [1]  5.141297642  1.829910381  1.213248010  1.018928407  0.859131851
    ##  [6]  0.702233999  0.388599832  0.293480958  0.007868141 -0.028475750
    ## [11] -0.426223470

``` r
S_sim_pos <- S_sim[,sim_cov_keep,sim_cov_keep]

S_sim_pos_ev <- plyr::aaply(S_sim_pos, 1, function(x) eigen(cov2cor(x))$values)

summary(S_sim_pos_ev)
```

    ##        1               2               3               4         
    ##  Min.   :4.298   Min.   :1.564   Min.   :1.127   Min.   :0.8009  
    ##  1st Qu.:4.937   1st Qu.:1.920   1st Qu.:1.349   1st Qu.:1.0508  
    ##  Median :5.289   Median :2.080   Median :1.480   Median :1.1451  
    ##  Mean   :5.490   Mean   :2.104   Mean   :1.477   Mean   :1.1533  
    ##  3rd Qu.:5.843   3rd Qu.:2.253   3rd Qu.:1.596   3rd Qu.:1.2494  
    ##  Max.   :9.541   Max.   :3.431   Max.   :1.867   Max.   :1.5283  
    ##        5                6                7                8           
    ##  Min.   :0.6220   Min.   :0.3968   Min.   :0.1447   Min.   :-0.03045  
    ##  1st Qu.:0.8046   1st Qu.:0.5608   1st Qu.:0.3170   1st Qu.: 0.07870  
    ##  Median :0.8677   Median :0.6339   Median :0.4178   Median : 0.17195  
    ##  Mean   :0.8762   Mean   :0.6320   Mean   :0.4063   Mean   : 0.16820  
    ##  3rd Qu.:0.9343   3rd Qu.:0.7219   3rd Qu.:0.4771   3rd Qu.: 0.24652  
    ##  Max.   :1.1353   Max.   :1.0520   Max.   :0.6923   Max.   : 0.46514  
    ##        9                  10                  11         
    ##  Min.   :-0.27138   Min.   :-0.914782   Min.   :-5.2784  
    ##  1st Qu.:-0.06645   1st Qu.:-0.416622   1st Qu.:-1.2814  
    ##  Median :-0.01120   Median :-0.193797   Median :-0.8558  
    ##  Mean   :-0.02187   Mean   :-0.280802   Mean   :-1.0052  
    ##  3rd Qu.: 0.03430   3rd Qu.:-0.116528   3rd Qu.:-0.5088  
    ##  Max.   : 0.17940   Max.   :-0.003535   Max.   :-0.2225

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
    ##       0.327       0.413       0.005       0.253       0.901       0.157 
    ##      PopDep      PopAnh   PopAppDec   PopAppInc   PopSleDec   PopSleInc 
    ##       0.089       0.047       0.451       0.467       0.408       0.712 
    ##  PopPsycInc  PopPsycDec    PopFatig    PopGuilt     PopConc      PopSui 
    ##       0.005       0.005       0.494       0.506       0.497       0.628 
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
    ## PopAnh       0.915  -0.215           0.254 
    ## PopAppDec            0.396  -0.173   0.596 
    ## PopAppInc    0.326  -0.649                 
    ## PopSleDec    0.409           0.621   0.198 
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
