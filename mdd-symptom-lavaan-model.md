CFA of MDD symptoms phenotypes
================
Mark Adams, Bradley Jermy, Jackson Thorp, Andrew Grotzinger, Michel Nivard

Phenotype level CFA models of MDD symptoms.

``` {.r}
library(dplyr)
library(lavaan)
```

# UKB CIDI

Load the phenotype file

``` {.r}
ukb_cidi <- readRDS("pheno/ukb_mhq_cidi_symptoms_screen_absent_present_contrasts.rds")
```

Rename symptoms to match the GenomicSEM analysis: UKB is one of the *Pop*ulation samples.

``` {.r}
ab0pr1 <- function(symptom) case_when(symptom %in% c('Absent', 'Not present') ~ 0L,
                                      symptom == 'Present' ~ 1L,
                                      TRUE ~ NA_integer_)

ukb_symp <-
as_tibble(ukb_cidi) %>%
transmute(f.eid,
PopDep=ab0pr1(cidi_mood),
PopAnh=ab0pr1(cidi_anhed),
PopAppDec=ab0pr1(cidi_wloss),
PopAppInc=ab0pr1(cidi_wgain),
PopSleDec=ab0pr1(cidi_insomn),
PopSleInc=ab0pr1(cidi_hypersom),
PopFatig=ab0pr1(cidi_tired),
PopGuilt=ab0pr1(cidi_worth),
PopConc=ab0pr1(cidi_cnctr),
PopSui=ab0pr1(cidi_death)) %>%
filter_at(vars(starts_with('Pop')), any_vars(!is.na(.)))
```

# Models

## Common factor

``` {.r}
pop_commonfactor.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
"

pop_commonfactor.fit <- cfa(model=pop_commonfactor.model, data=ukb_symp,  missing='ML')
                               
summary(pop_commonfactor.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 43 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         30
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                        
    ##   Test statistic                              13744.396
    ##   Degrees of freedom                                 35
    ##   P-value (Chi-square)                            0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.911
    ##   Tucker-Lewis Index (TLI)                       0.886
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -526424.400
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1052908.801
    ##   Bayesian (BIC)                           1053207.758
    ##   Sample-size adjusted Bayesian (BIC)      1053112.417
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.050
    ##   90 Percent confidence interval - lower         0.049
    ##   90 Percent confidence interval - upper         0.051
    ##   P-value RMSEA <= 0.05                          0.575
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.134
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.401    0.001  345.056    0.000
    ##     PopAnh            0.389    0.001  350.594    0.000
    ##     PopAppDec         0.143    0.003   42.490    0.000
    ##     PopAppInc         0.150    0.003   51.878    0.000
    ##     PopSleDec         0.227    0.003   75.905    0.000
    ##     PopSleInc         0.155    0.003   61.392    0.000
    ##     PopFatig          0.325    0.003  125.552    0.000
    ##     PopGuilt          0.350    0.003  108.035    0.000
    ##     PopConc           0.371    0.003  131.309    0.000
    ##     PopSui            0.200    0.003   61.969    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.921    0.000
    ##    .PopAnh            0.396    0.001  320.699    0.000
    ##    .PopAppDec         0.320    0.003  106.203    0.000
    ##    .PopAppInc         0.132    0.003   50.805    0.000
    ##    .PopSleDec         0.584    0.003  213.266    0.000
    ##    .PopSleInc         0.053    0.002   23.101    0.000
    ##    .PopFatig          0.580    0.002  238.662    0.000
    ##    .PopGuilt          0.261    0.003   87.638    0.000
    ##    .PopConc           0.515    0.003  194.614    0.000
    ##    .PopSui            0.379    0.003  130.765    0.000
    ##     A1                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##    .PopDep            0.087    0.001  163.810    0.000
    ##    .PopAnh            0.088    0.000  191.440    0.000
    ##    .PopAppDec         0.236    0.001  191.297    0.000
    ##    .PopAppInc         0.173    0.001  189.768    0.000
    ##    .PopSleDec         0.166    0.001  183.400    0.000
    ##    .PopSleInc         0.129    0.001  188.497    0.000
    ##    .PopFatig          0.105    0.001  162.521    0.000
    ##    .PopGuilt          0.198    0.001  180.225    0.000
    ##    .PopConc           0.111    0.001  152.363    0.000
    ##    .PopSui            0.233    0.001  198.067    0.000

Add negative residual correlations to directional symptoms

``` {.r}
pop_commonfactor_dir.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"

pop_commonfactor_dir.fit <- cfa(model=pop_commonfactor_dir.model, data=ukb_symp,  missing='ML')
                               
summary(pop_commonfactor_dir.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 50 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         32
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              8053.091
    ##   Degrees of freedom                                33
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.948
    ##   Tucker-Lewis Index (TLI)                       0.929
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -523578.748
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1047221.496
    ##   Bayesian (BIC)                           1047540.384
    ##   Sample-size adjusted Bayesian (BIC)      1047438.687
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.039
    ##   90 Percent confidence interval - lower         0.039
    ##   90 Percent confidence interval - upper         0.040
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.140
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.402    0.001  347.367    0.000
    ##     PopAnh            0.387    0.001  350.597    0.000
    ##     PopAppDec         0.163    0.003   48.409    0.000
    ##     PopAppInc         0.165    0.003   57.135    0.000
    ##     PopSleDec         0.244    0.003   80.878    0.000
    ##     PopSleInc         0.169    0.003   66.192    0.000
    ##     PopFatig          0.327    0.003  127.441    0.000
    ##     PopGuilt          0.350    0.003  108.572    0.000
    ##     PopConc           0.373    0.003  133.879    0.000
    ##     PopSui            0.203    0.003   62.848    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##  .PopAppDec ~~                                        
    ##    .PopAppInc        -0.048    0.001  -62.924    0.000
    ##  .PopSleDec ~~                                        
    ##    .PopSleInc        -0.021    0.001  -38.579    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.921    0.000
    ##    .PopAnh            0.396    0.001  320.698    0.000
    ##    .PopAppDec         0.306    0.003  100.862    0.000
    ##    .PopAppInc         0.121    0.003   46.285    0.000
    ##    .PopSleDec         0.572    0.003  206.972    0.000
    ##    .PopSleInc         0.043    0.002   18.224    0.000
    ##    .PopFatig          0.578    0.002  239.220    0.000
    ##    .PopGuilt          0.260    0.003   87.669    0.000
    ##    .PopConc           0.513    0.003  195.673    0.000
    ##    .PopSui            0.377    0.003  129.793    0.000
    ##     A1                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##    .PopDep            0.086    0.001  165.070    0.000
    ##    .PopAnh            0.089    0.000  197.998    0.000
    ##    .PopAppDec         0.233    0.001  190.450    0.000
    ##    .PopAppInc         0.171    0.001  188.957    0.000
    ##    .PopSleDec         0.163    0.001  181.532    0.000
    ##    .PopSleInc         0.127    0.001  187.001    0.000
    ##    .PopFatig          0.105    0.001  164.075    0.000
    ##    .PopGuilt          0.199    0.001  181.396    0.000
    ##    .PopConc           0.111    0.001  154.606    0.000
    ##    .PopSui            0.233    0.001  198.167    0.000

## Two-factor models

### Psychological-Somatic (Elhai Model 2a)

``` {.r}
pop_psych_soma.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui 
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_psych_soma.fit <- cfa(model=pop_psych_soma.model, data=ukb_symp, missing="ML")

summary(pop_psych_soma.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 51 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         33
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              7558.880
    ##   Degrees of freedom                                32
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.951
    ##   Tucker-Lewis Index (TLI)                       0.932
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -523331.643
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1046729.285
    ##   Bayesian (BIC)                           1047058.138
    ##   Sample-size adjusted Bayesian (BIC)      1046953.263
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.039
    ##   90 Percent confidence interval - lower         0.038
    ##   90 Percent confidence interval - upper         0.039
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.126
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.401    0.001  344.505    0.000
    ##     PopAnh            0.389    0.001  349.253    0.000
    ##     PopGuilt          0.352    0.003  108.128    0.000
    ##     PopConc           0.368    0.003  129.497    0.000
    ##     PopSui            0.203    0.003   62.643    0.000
    ##   A2 =~                                               
    ##     PopAppDec         0.158    0.003   48.061    0.000
    ##     PopAppInc         0.166    0.003   59.085    0.000
    ##     PopSleDec         0.243    0.003   82.591    0.000
    ##     PopSleInc         0.167    0.002   67.025    0.000
    ##     PopFatig          0.324    0.003  128.135    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##  .PopAppDec ~~                                        
    ##    .PopAppInc        -0.049    0.001  -64.424    0.000
    ##  .PopSleDec ~~                                        
    ##    .PopSleInc        -0.024    0.001  -42.008    0.000
    ##   A1 ~~                                               
    ##     A2                0.948    0.003  371.044    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.921    0.000
    ##    .PopAnh            0.396    0.001  320.699    0.000
    ##    .PopGuilt          0.260    0.003   86.821    0.000
    ##    .PopConc           0.517    0.003  194.477    0.000
    ##    .PopSui            0.377    0.003  129.729    0.000
    ##    .PopAppDec         0.315    0.003  107.779    0.000
    ##    .PopAppInc         0.127    0.003   50.606    0.000
    ##    .PopSleDec         0.582    0.003  214.950    0.000
    ##    .PopSleInc         0.050    0.002   22.011    0.000
    ##    .PopFatig          0.593    0.002  239.605    0.000
    ##     A1                0.000                           
    ##     A2                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##    .PopDep            0.087    0.001  162.581    0.000
    ##    .PopAnh            0.088    0.000  187.174    0.000
    ##    .PopGuilt          0.198    0.001  179.426    0.000
    ##    .PopConc           0.112    0.001  152.319    0.000
    ##    .PopSui            0.232    0.001  197.737    0.000
    ##    .PopAppDec         0.232    0.001  188.679    0.000
    ##    .PopAppInc         0.169    0.001  185.688    0.000
    ##    .PopSleDec         0.159    0.001  173.121    0.000
    ##    .PopSleInc         0.126    0.001  182.186    0.000
    ##    .PopFatig          0.098    0.001  138.569    0.000

``` {.r}
pop_psych_soma_bif.model <- "
A =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopConc + PopSui 
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig
A1 ~~ 1*A1
A2 ~~ 1*A2
A ~~ 0*A1 + 0*A2
A1 ~~ 0*A2
"
pop_psych_soma_bif.fit <- cfa(model=pop_psych_soma_bif.model, data=ukb_symp, missing="ML")

summary(pop_psych_soma_bif.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 62 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         41
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              4095.589
    ##   Degrees of freedom                                24
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.974
    ##   Tucker-Lewis Index (TLI)                       0.951
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -521599.997
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1043281.994
    ##   Bayesian (BIC)                           1043690.569
    ##   Sample-size adjusted Bayesian (BIC)      1043560.269
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.033
    ##   90 Percent confidence interval - lower         0.032
    ##   90 Percent confidence interval - upper         0.034
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.043
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A =~                                                
    ##     PopDep            0.796    4.890    0.163    0.871
    ##     PopAnh            0.900    5.529    0.163    0.871
    ##     PopAppDec         0.348    2.140    0.163    0.871
    ##     PopAppInc         0.339    2.081    0.163    0.871
    ##     PopSleDec         0.519    3.191    0.163    0.871
    ##     PopSleInc         0.319    1.958    0.163    0.871
    ##     PopFatig          0.737    4.530    0.163    0.871
    ##     PopGuilt          0.851    5.232    0.163    0.871
    ##     PopConc           0.833    5.120    0.163    0.871
    ##     PopSui            0.505    3.101    0.163    0.871
    ##   A1 =~                                               
    ##     PopDep            0.313    0.005   66.156    0.000
    ##     PopAnh            0.208    0.004   54.945    0.000
    ##     PopGuilt          0.159    0.005   30.399    0.000
    ##     PopConc          -0.008    0.005   -1.626    0.104
    ##     PopSui            0.113    0.005   22.504    0.000
    ##   A2 =~                                               
    ##     PopAppDec         0.264    0.006   46.989    0.000
    ##     PopAppInc        -0.176    0.004  -47.710    0.000
    ##     PopSleDec         0.070    0.003   26.409    0.000
    ##     PopSleInc        -0.068    0.003  -26.652    0.000
    ##     PopFatig         -0.020    0.002   -9.410    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A ~~                                                
    ##     A1                0.000                           
    ##     A2                0.000                           
    ##   A1 ~~                                               
    ##     A2                0.000                           
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.920    0.000
    ##    .PopAnh            0.396    0.001  320.706    0.000
    ##    .PopAppDec         0.356    0.003  135.960    0.000
    ##    .PopAppInc         0.176    0.002   77.753    0.000
    ##    .PopSleDec         0.648    0.003  238.687    0.000
    ##    .PopSleInc         0.104    0.002   50.605    0.000
    ##    .PopFatig          0.673    0.003  225.488    0.000
    ##    .PopGuilt          0.266    0.004   62.586    0.000
    ##    .PopConc           0.626    0.004  176.845    0.000
    ##    .PopSui            0.366    0.004   85.845    0.000
    ##     A                 0.000                           
    ##     A1                0.000                           
    ##     A2                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##    .PopDep            0.066    0.001   47.699    0.000
    ##    .PopAnh            0.089    0.001  138.022    0.000
    ##    .PopAppDec         0.163    0.003   54.847    0.000
    ##    .PopAppInc         0.140    0.001  100.114    0.000
    ##    .PopSleDec         0.156    0.001  166.262    0.000
    ##    .PopSleInc         0.124    0.001  171.782    0.000
    ##    .PopFatig          0.095    0.001  134.113    0.000
    ##    .PopGuilt          0.197    0.001  174.914    0.000
    ##    .PopConc           0.098    0.001  108.970    0.000
    ##    .PopSui            0.231    0.001  193.245    0.000
    ##     A                 0.132    1.622    0.081    0.935

### Psychological-Neurovegetative (Elhai Model 2b)

``` {.r}
pop_psych_veg.model <- "
A1 =~ NA*PopDep + PopAnh + PopGuilt + PopSui
A2 =~ NA*PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_psych_veg.fit <- cfa(model=pop_psych_veg.model, data=ukb_symp, missing="ML")

summary(pop_psych_veg.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 50 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         33
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              4929.103
    ##   Degrees of freedom                                32
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.968
    ##   Tucker-Lewis Index (TLI)                       0.956
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -522016.754
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1044099.508
    ##   Bayesian (BIC)                           1044428.361
    ##   Sample-size adjusted Bayesian (BIC)      1044323.485
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.031
    ##   90 Percent confidence interval - lower         0.030
    ##   90 Percent confidence interval - upper         0.032
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.077
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.398    0.001  327.347    0.000
    ##     PopAnh            0.398    0.001  341.083    0.000
    ##     PopGuilt          0.348    0.003   99.611    0.000
    ##     PopSui            0.200    0.003   60.422    0.000
    ##   A2 =~                                               
    ##     PopAppDec         0.143    0.003   48.611    0.000
    ##     PopAppInc         0.143    0.003   55.822    0.000
    ##     PopSleDec         0.228    0.003   86.203    0.000
    ##     PopSleInc         0.144    0.002   63.046    0.000
    ##     PopFatig          0.301    0.002  129.491    0.000
    ##     PopConc           0.344    0.003  135.651    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##  .PopAppDec ~~                                        
    ##    .PopAppInc        -0.049    0.001  -63.651    0.000
    ##  .PopSleDec ~~                                        
    ##    .PopSleInc        -0.023    0.001  -41.789    0.000
    ##   A1 ~~                                               
    ##     A2                0.875    0.003  286.206    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.921    0.000
    ##    .PopAnh            0.396    0.001  320.705    0.000
    ##    .PopGuilt          0.265    0.003   84.093    0.000
    ##    .PopSui            0.381    0.003  129.327    0.000
    ##    .PopAppDec         0.334    0.003  126.710    0.000
    ##    .PopAppInc         0.151    0.002   65.632    0.000
    ##    .PopSleDec         0.605    0.002  243.391    0.000
    ##    .PopSleInc         0.074    0.002   35.570    0.000
    ##    .PopFatig          0.626    0.002  264.652    0.000
    ##    .PopConc           0.567    0.003  218.345    0.000
    ##     A1                0.000                           
    ##     A2                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##    .PopDep            0.089    0.001  148.981    0.000
    ##    .PopAnh            0.081    0.001  148.020    0.000
    ##    .PopGuilt          0.196    0.001  172.831    0.000
    ##    .PopSui            0.232    0.001  196.227    0.000
    ##    .PopAppDec         0.232    0.001  189.298    0.000
    ##    .PopAppInc         0.170    0.001  187.579    0.000
    ##    .PopSleDec         0.158    0.001  175.023    0.000
    ##    .PopSleInc         0.127    0.001  184.920    0.000
    ##    .PopFatig          0.097    0.001  145.635    0.000
    ##    .PopConc           0.100    0.001  131.301    0.000

``` {.r}
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
"

pop_psych_veg_bif.fit <- cfa(model=pop_psych_veg_bif.model, data=ukb_symp, missing="ML")

summary(pop_psych_veg_bif.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 54 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         40
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              3790.216
    ##   Degrees of freedom                                25
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.976
    ##   Tucker-Lewis Index (TLI)                       0.956
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -521447.310
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1042974.621
    ##   Bayesian (BIC)                           1043373.231
    ##   Sample-size adjusted Bayesian (BIC)      1043246.109
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.031
    ##   90 Percent confidence interval - lower         0.030
    ##   90 Percent confidence interval - upper         0.032
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.050
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.324    0.004   78.343    0.000
    ##     PopAnh            0.212    0.003   65.004    0.000
    ##     PopGuilt          0.159    0.005   32.380    0.000
    ##     PopSui            0.119    0.005   23.921    0.000
    ##   A2 =~                                               
    ##     PopAppDec         0.243    0.004   56.912    0.000
    ##     PopAppInc        -0.179    0.004  -50.653    0.000
    ##     PopSleDec         0.093    0.003   30.179    0.000
    ##     PopSleInc        -0.067    0.002  -28.025    0.000
    ##     PopFatig          0.006    0.003    2.376    0.017
    ##     PopConc           0.053    0.003   17.616    0.000
    ##   A =~                                                
    ##     PopDep            0.281    0.004   78.015    0.000
    ##     PopAnh            0.323    0.002  145.966    0.000
    ##     PopGuilt          0.309    0.004   86.580    0.000
    ##     PopSui            0.181    0.003   54.373    0.000
    ##     PopAppDec         0.094    0.003   28.926    0.000
    ##     PopAppInc         0.144    0.003   54.325    0.000
    ##     PopSleDec         0.178    0.003   70.203    0.000
    ##     PopSleInc         0.122    0.002   58.909    0.000
    ##     PopFatig          0.265    0.002  114.484    0.000
    ##     PopConc           0.299    0.003  114.861    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 ~~                                               
    ##     A                 0.000                           
    ##   A2 ~~                                               
    ##     A                 0.000                           
    ##   A1 ~~                                               
    ##     A2                0.000                           
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.920    0.000
    ##    .PopAnh            0.396    0.001  320.706    0.000
    ##    .PopGuilt          0.267    0.004   62.606    0.000
    ##    .PopSui            0.365    0.004   84.694    0.000
    ##    .PopAppDec         0.375    0.003  144.940    0.000
    ##    .PopAppInc         0.167    0.002   71.518    0.000
    ##    .PopSleDec         0.656    0.003  259.810    0.000
    ##    .PopSleInc         0.102    0.002   51.136    0.000
    ##    .PopFatig          0.678    0.003  243.102    0.000
    ##    .PopConc           0.628    0.003  201.251    0.000
    ##     A1                0.000                           
    ##     A2                0.000                           
    ##     A                 0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##     A                 1.000                           
    ##    .PopDep            0.064    0.001   44.929    0.000
    ##    .PopAnh            0.090    0.001  140.033    0.000
    ##    .PopGuilt          0.196    0.001  175.369    0.000
    ##    .PopSui            0.231    0.001  193.222    0.000
    ##    .PopAppDec         0.178    0.002   82.948    0.000
    ##    .PopAppInc         0.135    0.001   94.558    0.000
    ##    .PopSleDec         0.155    0.001  161.801    0.000
    ##    .PopSleInc         0.123    0.001  167.750    0.000
    ##    .PopFatig          0.095    0.001  142.935    0.000
    ##    .PopConc           0.097    0.001  122.971    0.000

### Affective-Neurovegetative (Elhai Model 2c)

``` {.r}
pop_affect_veg.model <- "
A1 =~ NA*PopDep + PopGuilt + PopSui
A2 =~ NA*PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"
pop_affect_veg.fit <- cfa(model=pop_affect_veg.model, data=ukb_symp, missing="ML")

summary(pop_affect_veg.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 52 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         33
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              7872.353
    ##   Degrees of freedom                                32
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.949
    ##   Tucker-Lewis Index (TLI)                       0.929
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -523488.379
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1047042.758
    ##   Bayesian (BIC)                           1047371.612
    ##   Sample-size adjusted Bayesian (BIC)      1047266.736
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.039
    ##   90 Percent confidence interval - lower         0.039
    ##   90 Percent confidence interval - upper         0.040
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.144
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.412    0.001  315.781    0.000
    ##     PopGuilt          0.392    0.004   87.059    0.000
    ##     PopSui            0.229    0.004   58.218    0.000
    ##   A2 =~                                               
    ##     PopAnh            0.388    0.001  350.192    0.000
    ##     PopAppInc         0.164    0.003   57.145    0.000
    ##     PopAppDec         0.162    0.003   48.328    0.000
    ##     PopSleInc         0.168    0.003   66.173    0.000
    ##     PopSleDec         0.242    0.003   80.622    0.000
    ##     PopFatig          0.327    0.003  127.527    0.000
    ##     PopConc           0.371    0.003  133.529    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##  .PopAppInc ~~                                        
    ##    .PopAppDec        -0.048    0.001  -62.917    0.000
    ##  .PopSleInc ~~                                        
    ##    .PopSleDec        -0.021    0.001  -38.572    0.000
    ##   A1 ~~                                               
    ##     A2                0.972    0.002  506.580    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.920    0.000
    ##    .PopGuilt          0.226    0.004   56.735    0.000
    ##    .PopSui            0.356    0.003  103.452    0.000
    ##    .PopAnh            0.396    0.001  320.697    0.000
    ##    .PopAppInc         0.122    0.003   46.892    0.000
    ##    .PopAppDec         0.307    0.003  101.765    0.000
    ##    .PopSleInc         0.044    0.002   18.726    0.000
    ##    .PopSleDec         0.573    0.003  208.227    0.000
    ##    .PopFatig          0.579    0.002  240.165    0.000
    ##    .PopConc           0.514    0.003  196.548    0.000
    ##     A1                0.000                           
    ##     A2                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##    .PopDep            0.078    0.001  107.717    0.000
    ##    .PopGuilt          0.192    0.001  159.801    0.000
    ##    .PopSui            0.230    0.001  193.241    0.000
    ##    .PopAnh            0.089    0.000  193.984    0.000
    ##    .PopAppInc         0.171    0.001  188.950    0.000
    ##    .PopAppDec         0.233    0.001  190.460    0.000
    ##    .PopSleInc         0.127    0.001  186.991    0.000
    ##    .PopSleDec         0.163    0.001  181.487    0.000
    ##    .PopFatig          0.104    0.001  163.470    0.000
    ##    .PopConc           0.111    0.001  154.048    0.000

Bifactor model

``` {.r}
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
"
pop_affect_veg_bif.fit <- cfa(model=pop_affect_veg_bif.model, data=ukb_symp, missing="ML")

summary(pop_affect_veg_bif.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 68 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         40
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              7497.803
    ##   Degrees of freedom                                25
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.952
    ##   Tucker-Lewis Index (TLI)                       0.913
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -523301.104
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1046682.208
    ##   Bayesian (BIC)                           1047080.818
    ##   Sample-size adjusted Bayesian (BIC)      1046953.697
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.044
    ##   90 Percent confidence interval - lower         0.043
    ##   90 Percent confidence interval - upper         0.044
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.138
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.051    0.006    8.328    0.000
    ##     PopGuilt          0.075    0.009    8.541    0.000
    ##     PopSui            0.172    0.019    8.933    0.000
    ##   A2 =~                                               
    ##     PopAnh           -0.000    0.002   -0.121    0.904
    ##     PopAppInc        -0.163    0.003  -48.216    0.000
    ##     PopAppDec         0.263    0.005   54.128    0.000
    ##     PopSleInc        -0.064    0.003  -24.761    0.000
    ##     PopSleDec         0.092    0.003   27.790    0.000
    ##     PopFatig          0.007    0.003    2.126    0.034
    ##     PopConc           0.046    0.003   13.300    0.000
    ##   A =~                                                
    ##     PopDep            0.400    0.001  341.380    0.000
    ##     PopGuilt          0.361    0.004   81.602    0.000
    ##     PopSui            0.225    0.005   49.314    0.000
    ##     PopAnh            0.389    0.001  349.618    0.000
    ##     PopAppInc         0.176    0.004   48.366    0.000
    ##     PopAppDec         0.134    0.005   26.804    0.000
    ##     PopSleInc         0.163    0.003   61.174    0.000
    ##     PopSleDec         0.223    0.003   69.082    0.000
    ##     PopFatig          0.326    0.003  125.695    0.000
    ##     PopConc           0.368    0.003  127.991    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 ~~                                               
    ##     A                 0.000                           
    ##   A2 ~~                                               
    ##     A                 0.000                           
    ##   A1 ~~                                               
    ##     A2                0.000                           
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.547    0.001  435.921    0.000
    ##    .PopGuilt          0.248    0.004   55.593    0.000
    ##    .PopSui            0.349    0.005   76.196    0.000
    ##    .PopAnh            0.396    0.001  320.698    0.000
    ##    .PopAppInc         0.114    0.003   37.638    0.000
    ##    .PopAppDec         0.327    0.004   83.075    0.000
    ##    .PopSleInc         0.047    0.002   19.834    0.000
    ##    .PopSleDec         0.587    0.003  205.369    0.000
    ##    .PopFatig          0.580    0.002  238.574    0.000
    ##    .PopConc           0.517    0.003  193.995    0.000
    ##     A1                0.000                           
    ##     A2                0.000                           
    ##     A                 0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##     A                 1.000                           
    ##    .PopDep            0.085    0.001  116.327    0.000
    ##    .PopGuilt          0.194    0.002  120.536    0.000
    ##    .PopSui            0.205    0.007   30.915    0.000
    ##    .PopAnh            0.088    0.000  188.446    0.000
    ##    .PopAppInc         0.143    0.001  115.823    0.000
    ##    .PopAppDec         0.167    0.003   64.861    0.000
    ##    .PopSleInc         0.124    0.001  170.180    0.000
    ##    .PopSleDec         0.158    0.001  162.836    0.000
    ##    .PopFatig          0.104    0.001  162.258    0.000
    ##    .PopConc           0.109    0.001  147.421    0.000

## 3 factor models

``` {.r}
pop_cog_mood_neuroveg.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
PopSleDec ~~ PopSleInc
PopAppDec ~~ PopAppInc
"
pop_cog_mood_neuroveg.fit <- cfa(model=pop_cog_mood_neuroveg.model, data=ukb_symp, missing="ML")
```

    ## Warning in lav_object_post_check(object): lavaan WARNING: covariance matrix of latent variables
    ##                 is not positive definite;
    ##                 use lavInspect(fit, "cov.lv") to investigate.

``` {.r}
summary(pop_cog_mood_neuroveg.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 57 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         36
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              4568.581
    ##   Degrees of freedom                                29
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.971
    ##   Tucker-Lewis Index (TLI)                       0.954
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -521836.493
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1043744.986
    ##   Bayesian (BIC)                           1044103.735
    ##   Sample-size adjusted Bayesian (BIC)      1043989.325
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.032
    ##   90 Percent confidence interval - lower         0.031
    ##   90 Percent confidence interval - upper         0.032
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.060
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopGuilt          0.167    0.005   30.739    0.000
    ##     PopConc           0.309    0.003  114.194    0.000
    ##     PopSui            0.155    0.003   56.723    0.000
    ##   A2 =~                                               
    ##     PopDep            0.386    0.001  264.707    0.000
    ##     PopAnh            0.414    0.001  281.333    0.000
    ##     PopGuilt          0.135    0.007   19.733    0.000
    ##   A3 =~                                               
    ##     PopSleDec         0.207    0.003   81.401    0.000
    ##     PopSleInc         0.134    0.002   63.118    0.000
    ##     PopFatig          0.274    0.002  117.753    0.000
    ##     PopAppDec         0.131    0.003   47.581    0.000
    ##     PopAppInc         0.135    0.002   56.998    0.000
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##  .PopSleDec ~~                                        
    ##    .PopSleInc        -0.023    0.001  -41.642    0.000
    ##  .PopAppDec ~~                                        
    ##    .PopAppInc        -0.049    0.001  -64.027    0.000
    ##   A1 ~~                                               
    ##     A2                0.840    0.007  119.498    0.000
    ##     A3                1.013    0.005  185.055    0.000
    ##   A2 ~~                                               
    ##     A3                0.802    0.006  141.638    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopGuilt          0.321    0.004   89.272    0.000
    ##    .PopConc           0.601    0.003  197.072    0.000
    ##    .PopSui            0.432    0.003  163.154    0.000
    ##    .PopDep            0.547    0.001  435.921    0.000
    ##    .PopAnh            0.396    0.001  320.707    0.000
    ##    .PopSleDec         0.632    0.003  249.041    0.000
    ##    .PopSleInc         0.089    0.002   43.952    0.000
    ##    .PopFatig          0.661    0.003  253.147    0.000
    ##    .PopAppDec         0.350    0.003  138.930    0.000
    ##    .PopAppInc         0.165    0.002   74.505    0.000
    ##     A1                0.000                           
    ##     A2                0.000                           
    ##     A3                0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##     A3                1.000                           
    ##    .PopGuilt          0.203    0.001  184.262    0.000
    ##    .PopConc           0.105    0.001  102.739    0.000
    ##    .PopSui            0.234    0.001  197.700    0.000
    ##    .PopDep            0.099    0.001  115.273    0.000
    ##    .PopAnh            0.068    0.001   72.906    0.000
    ##    .PopSleDec         0.158    0.001  173.937    0.000
    ##    .PopSleInc         0.127    0.001  184.086    0.000
    ##    .PopFatig          0.097    0.001  139.480    0.000
    ##    .PopAppDec         0.232    0.001  189.288    0.000
    ##    .PopAppInc         0.170    0.001  186.714    0.000

``` {.r}
pop_cog_mood_neuroveg_bif.model <- "
A1 =~ NA*PopGuilt + PopConc + PopSui
A2 =~ NA*PopDep + PopAnh + PopGuilt
A3 =~ NA*PopSleDec + PopSleInc + PopFatig + PopAppDec + PopAppInc
A =~ NA*PopDep + PopGuilt + PopSui + PopAnh + PopAppInc + PopAppDec + PopSleInc + PopSleDec + PopFatig + PopConc
A1 ~~ 1*A1
A2 ~~ 1*A2
A3 ~~ 1*A3
A ~~ 0*A1 + 0*A2 + 0*A3
A1 ~~ 0*A2 + 0*A3
A2 ~~ 0*A3
"
pop_cog_mood_neuroveg_bif.fit <- cfa(model=pop_cog_mood_neuroveg_bif.model, data=ukb_symp, missing="ML")

summary(pop_cog_mood_neuroveg_bif.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 113 iterations
    ## 
    ##   Estimator                                         ML
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         42
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                       
    ##   Test statistic                              4244.884
    ##   Degrees of freedom                                23
    ##   P-value (Chi-square)                           0.000
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            154800.316
    ##   Degrees of freedom                                45
    ##   P-value                                        0.000
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.973
    ##   Tucker-Lewis Index (TLI)                       0.947
    ## 
    ## Loglikelihood and Information Criteria:
    ## 
    ##   Loglikelihood user model (H0)            -521674.644
    ##   Loglikelihood unrestricted model (H1)    -519552.202
    ##                                                       
    ##   Akaike (AIC)                             1043433.289
    ##   Bayesian (BIC)                           1043851.829
    ##   Sample-size adjusted Bayesian (BIC)      1043718.351
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.034
    ##   90 Percent confidence interval - lower         0.033
    ##   90 Percent confidence interval - upper         0.035
    ##   P-value RMSEA <= 0.05                          1.000
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.054
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                             Standard
    ##   Information                                 Observed
    ##   Observed information based on                Hessian
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopGuilt          0.037    0.020    1.839    0.066
    ##     PopConc          -0.019    0.011   -1.753    0.080
    ##     PopSui            0.359    0.199    1.806    0.071
    ##   A2 =~                                               
    ##     PopDep            0.292    0.005   59.900    0.000
    ##     PopAnh            0.211    0.004   59.614    0.000
    ##     PopGuilt          0.134    0.005   26.356    0.000
    ##   A3 =~                                               
    ##     PopSleDec         0.069    0.003   25.966    0.000
    ##     PopSleInc        -0.068    0.003  -27.097    0.000
    ##     PopFatig         -0.021    0.002   -9.962    0.000
    ##     PopAppDec         0.262    0.006   47.362    0.000
    ##     PopAppInc        -0.178    0.004  -47.981    0.000
    ##   A =~                                                
    ##     PopDep            0.822   10.355    0.079    0.937
    ##     PopGuilt          0.820   10.326    0.079    0.937
    ##     PopSui            0.417    5.254    0.079    0.937
    ##     PopAnh            0.916   11.545    0.079    0.937
    ##     PopAppInc         0.343    4.326    0.079    0.937
    ##     PopAppDec         0.362    4.555    0.079    0.937
    ##     PopSleInc         0.323    4.073    0.079    0.937
    ##     PopSleDec         0.530    6.677    0.079    0.937
    ##     PopFatig          0.742    9.351    0.079    0.937
    ##     PopConc           0.859   10.824    0.079    0.937
    ## 
    ## Covariances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 ~~                                               
    ##     A                 0.000                           
    ##   A2 ~~                                               
    ##     A                 0.000                           
    ##   A3 ~~                                               
    ##     A                 0.000                           
    ##   A1 ~~                                               
    ##     A2                0.000                           
    ##     A3                0.000                           
    ##   A2 ~~                                               
    ##     A3                0.000                           
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopGuilt          0.286    0.004   63.809    0.000
    ##    .PopConc           0.614    0.003  198.129    0.000
    ##    .PopSui            0.442    0.003  170.425    0.000
    ##    .PopDep            0.547    0.001  435.921    0.000
    ##    .PopAnh            0.396    0.001  320.707    0.000
    ##    .PopSleDec         0.644    0.003  253.473    0.000
    ##    .PopSleInc         0.102    0.002   51.416    0.000
    ##    .PopFatig          0.669    0.003  245.117    0.000
    ##    .PopAppDec         0.353    0.003  137.285    0.000
    ##    .PopAppInc         0.174    0.002   78.495    0.000
    ##     A1                0.000                           
    ##     A2                0.000                           
    ##     A3                0.000                           
    ##     A                 0.000                           
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##     A2                1.000                           
    ##     A3                1.000                           
    ##    .PopGuilt          0.199    0.002   96.321    0.000
    ##    .PopConc           0.098    0.001  100.018    0.000
    ##    .PopSui            0.105    0.142    0.734    0.463
    ##    .PopDep            0.075    0.002   40.510    0.000
    ##    .PopAnh            0.085    0.001   88.047    0.000
    ##    .PopSleDec         0.156    0.001  167.074    0.000
    ##    .PopSleInc         0.124    0.001  171.941    0.000
    ##    .PopFatig          0.096    0.001  141.715    0.000
    ##    .PopAppDec         0.163    0.003   55.799    0.000
    ##    .PopAppInc         0.139    0.001   99.066    0.000
    ##     A                 0.130    3.282    0.040    0.968

# Model comparison

``` {.r}
fits <- list(
pop_commonfactor.fit,
pop_commonfactor_dir.fit,

pop_psych_soma_bif.fit,
pop_psych_soma.fit,

pop_psych_veg_bif.fit,
pop_psych_veg.fit,

pop_affect_veg_bif.fit,
pop_affect_veg.fit,

pop_cog_mood_neuroveg.fit,
pop_cog_mood_neuroveg_bif.fit)

model_fits <- 
data.frame(Model=c('1a', '1b', '2a', '2a(ii)', '2b(i)', '2b(ii)', '2c(i)', '2c(ii)', '3', '3(ii)'),
   Name=c('Common',
          'Common (direction)',
          'Psych-Somatic',
          'Psych-Somatic (BiF)',
          'Psych-Neuroveg',
          'Psych-Neuroveg (BiF)',
          'Affect-Neuroveg',
          'Affect-Neuroveg (BiF)',
          'Cog-Mood-Neuroveg',
          'Cog-Mood-Neuroveg (BiF)')) %>%
bind_cols(
bind_rows(lapply(fits, function(fit) fitMeasures(fit)[c('aic', 'bic', 'cfi', 'srmr')]))
)

model_fits %>%
mutate(dAIC=aic-min(aic)) %>%
select(-bic)
```

    ##     Model                    Name     aic       cfi       srmr      dAIC
    ## 1      1a                  Common 1052909 0.9114124 0.13443218 9934.1803
    ## 2      1b      Common (direction) 1047221 0.9481757 0.14027941 4246.8754
    ## 3      2a           Psych-Somatic 1043282 0.9736902 0.04269307  307.3729
    ## 4  2a(ii)     Psych-Somatic (BiF) 1046729 0.9513627 0.12594659 3754.6644
    ## 5   2b(i)          Psych-Neuroveg 1042975 0.9756699 0.05019574    0.0000
    ## 6  2b(ii)    Psych-Neuroveg (BiF) 1044100 0.9683558 0.07676342 1124.8868
    ## 7   2c(i)         Affect-Neuroveg 1046682 0.9517121 0.13760998 3707.5875
    ## 8  2c(ii)   Affect-Neuroveg (BiF) 1047043 0.9493371 0.14436362 4068.1377
    ## 9       3       Cog-Mood-Neuroveg 1043745 0.9706661 0.05950532  770.3651
    ## 10  3(ii) Cog-Mood-Neuroveg (BiF) 1043433 0.9727190 0.05434575  458.6680
