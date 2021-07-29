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
                               
standardizedSolution(pop_commonfactor.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~    PopDep   0.806 0.001 578.784      0    0.803    0.809
    ## 2         A1 =~    PopAnh   0.795 0.001 619.588      0    0.792    0.797
    ## 3         A1 =~ PopAppDec   0.282 0.006  45.376      0    0.270    0.294
    ## 4         A1 =~ PopAppInc   0.339 0.006  57.140      0    0.327    0.350
    ## 5         A1 =~ PopSleDec   0.487 0.005  93.871      0    0.477    0.498
    ## 6         A1 =~ PopSleInc   0.395 0.006  70.256      0    0.384    0.406
    ## 7         A1 =~  PopFatig   0.709 0.003 213.080      0    0.703    0.716
    ## 8         A1 =~  PopGuilt   0.618 0.004 156.640      0    0.610    0.626
    ## 9         A1 =~   PopConc   0.744 0.003 238.986      0    0.738    0.750
    ## 10        A1 =~    PopSui   0.383 0.005  70.255      0    0.372    0.394
    ## 11        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 12    PopDep ~~    PopDep   0.350 0.002 156.038      0    0.346    0.355
    ## 13    PopAnh ~~    PopAnh   0.368 0.002 180.671      0    0.364    0.372
    ## 14 PopAppDec ~~ PopAppDec   0.920 0.004 262.561      0    0.914    0.927
    ## 15 PopAppInc ~~ PopAppInc   0.885 0.004 220.268      0    0.877    0.893
    ## 16 PopSleDec ~~ PopSleDec   0.762 0.005 150.600      0    0.752    0.772
    ## 17 PopSleInc ~~ PopSleInc   0.844 0.004 189.530      0    0.835    0.852
    ## 18  PopFatig ~~  PopFatig   0.497 0.005 105.143      0    0.487    0.506
    ## 19  PopGuilt ~~  PopGuilt   0.618 0.005 126.629      0    0.608    0.627
    ## 20   PopConc ~~   PopConc   0.446 0.005  96.343      0    0.437    0.455
    ## 21    PopSui ~~    PopSui   0.853 0.004 204.431      0    0.845    0.862
    ## 22    PopDep ~1             1.100 0.003 344.014      0    1.094    1.106
    ## 23    PopAnh ~1             0.809 0.003 278.380      0    0.804    0.815
    ## 24 PopAppDec ~1             0.633 0.007  90.720      0    0.620    0.647
    ## 25 PopAppInc ~1             0.299 0.006  46.783      0    0.287    0.312
    ## 26 PopSleDec ~1             1.253 0.009 136.066      0    1.235    1.271
    ## 27 PopSleInc ~1             0.136 0.006  22.150      0    0.124    0.148
    ## 28  PopFatig ~1             1.264 0.009 134.983      0    1.246    1.283
    ## 29  PopGuilt ~1             0.460 0.006  71.171      0    0.448    0.473
    ## 30   PopConc ~1             1.033 0.009 118.548      0    1.016    1.050
    ## 31    PopSui ~1             0.726 0.007 104.103      0    0.712    0.739
    ## 32        A1 ~1             0.000 0.000      NA     NA    0.000    0.000

Add negative residual correlations to directional symptoms

``` {.r}
pop_commonfactor_dir.model <- "
A1 =~ NA*PopDep + PopAnh + PopAppDec + PopAppInc + PopSleDec + PopSleInc + PopFatig + PopGuilt + PopConc + PopSui
A1 ~~ 1*A1
PopAppDec ~~ PopAppInc
PopSleDec ~~ PopSleInc
"

pop_commonfactor_dir.fit <- cfa(model=pop_commonfactor_dir.model, data=ukb_symp,  missing='ML')
                               
standardizedSolution(pop_commonfactor_dir.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~    PopDep   0.807 0.001 588.386      0    0.805    0.810
    ## 2         A1 =~    PopAnh   0.791 0.001 624.105      0    0.789    0.794
    ## 3         A1 =~ PopAppDec   0.320 0.006  52.754      0    0.308    0.332
    ## 4         A1 =~ PopAppInc   0.371 0.006  64.283      0    0.360    0.382
    ## 5         A1 =~ PopSleDec   0.517 0.005 103.204      0    0.507    0.527
    ## 6         A1 =~ PopSleInc   0.428 0.006  77.687      0    0.417    0.439
    ## 7         A1 =~  PopFatig   0.711 0.003 218.104      0    0.705    0.718
    ## 8         A1 =~  PopGuilt   0.618 0.004 157.763      0    0.610    0.625
    ## 9         A1 =~   PopConc   0.746 0.003 246.103      0    0.740    0.751
    ## 10        A1 =~    PopSui   0.388 0.005  71.555      0    0.377    0.398
    ## 11        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 12 PopAppDec ~~ PopAppInc  -0.240 0.004 -67.644      0   -0.247   -0.233
    ## 13 PopSleDec ~~ PopSleInc  -0.148 0.004 -39.150      0   -0.156   -0.141
    ## 14    PopDep ~~    PopDep   0.348 0.002 157.011      0    0.344    0.352
    ## 15    PopAnh ~~    PopAnh   0.374 0.002 186.639      0    0.370    0.378
    ## 16 PopAppDec ~~ PopAppDec   0.898 0.004 231.403      0    0.890    0.905
    ## 17 PopAppInc ~~ PopAppInc   0.862 0.004 201.319      0    0.854    0.871
    ## 18 PopSleDec ~~ PopSleDec   0.732 0.005 141.221      0    0.722    0.743
    ## 19 PopSleInc ~~ PopSleInc   0.817 0.005 173.282      0    0.808    0.826
    ## 20  PopFatig ~~  PopFatig   0.494 0.005 106.545      0    0.485    0.503
    ## 21  PopGuilt ~~  PopGuilt   0.618 0.005 127.800      0    0.609    0.628
    ## 22   PopConc ~~   PopConc   0.444 0.005  98.341      0    0.435    0.453
    ## 23    PopSui ~~    PopSui   0.850 0.004 202.091      0    0.841    0.858
    ## 24    PopDep ~1             1.100 0.003 344.014      0    1.094    1.106
    ## 25    PopAnh ~1             0.809 0.003 278.379      0    0.804    0.815
    ## 26 PopAppDec ~1             0.600 0.007  85.611      0    0.586    0.614
    ## 27 PopAppInc ~1             0.272 0.006  42.680      0    0.259    0.284
    ## 28 PopSleDec ~1             1.213 0.009 131.446      0    1.195    1.231
    ## 29 PopSleInc ~1             0.108 0.006  17.583      0    0.096    0.120
    ## 30  PopFatig ~1             1.257 0.009 135.463      0    1.239    1.275
    ## 31  PopGuilt ~1             0.459 0.006  71.214      0    0.446    0.472
    ## 32   PopConc ~1             1.026 0.009 119.409      0    1.009    1.043
    ## 33    PopSui ~1             0.720 0.007 103.214      0    0.707    0.734
    ## 34        A1 ~1             0.000 0.000      NA     NA    0.000    0.000

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

standardizedSolution(pop_psych_soma.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~    PopDep   0.806 0.001 575.920      0    0.804    0.809
    ## 2         A1 =~    PopAnh   0.796 0.001 611.091      0    0.793    0.798
    ## 3         A1 =~  PopGuilt   0.621 0.004 157.299      0    0.613    0.629
    ## 4         A1 =~   PopConc   0.740 0.003 233.398      0    0.734    0.747
    ## 5         A1 =~    PopSui   0.388 0.005  71.253      0    0.377    0.399
    ## 6         A2 =~ PopAppDec   0.312 0.006  51.880      0    0.300    0.324
    ## 7         A2 =~ PopAppInc   0.374 0.006  66.098      0    0.362    0.385
    ## 8         A2 =~ PopSleDec   0.520 0.005 104.047      0    0.510    0.530
    ## 9         A2 =~ PopSleInc   0.427 0.005  77.808      0    0.416    0.438
    ## 10        A2 =~  PopFatig   0.719 0.003 213.801      0    0.712    0.725
    ## 11        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 12        A2 ~~        A2   1.000 0.000      NA     NA    1.000    1.000
    ## 13 PopAppDec ~~ PopAppInc  -0.249 0.004 -69.065      0   -0.256   -0.241
    ## 14 PopSleDec ~~ PopSleInc  -0.168 0.004 -42.263      0   -0.176   -0.160
    ## 15    PopDep ~~    PopDep   0.350 0.002 154.953      0    0.345    0.354
    ## 16    PopAnh ~~    PopAnh   0.367 0.002 176.957      0    0.363    0.371
    ## 17  PopGuilt ~~  PopGuilt   0.614 0.005 125.331      0    0.605    0.624
    ## 18   PopConc ~~   PopConc   0.452 0.005  96.136      0    0.442    0.461
    ## 19    PopSui ~~    PopSui   0.850 0.004 201.139      0    0.841    0.858
    ## 20 PopAppDec ~~ PopAppDec   0.903 0.004 240.966      0    0.895    0.910
    ## 21 PopAppInc ~~ PopAppInc   0.860 0.004 203.800      0    0.852    0.869
    ## 22 PopSleDec ~~ PopSleDec   0.730 0.005 140.497      0    0.720    0.740
    ## 23 PopSleInc ~~ PopSleInc   0.818 0.005 174.352      0    0.808    0.827
    ## 24  PopFatig ~~  PopFatig   0.483 0.005 100.028      0    0.474    0.493
    ## 25        A1 ~~        A2   0.948 0.003 371.044      0    0.943    0.953
    ## 26    PopDep ~1             1.100 0.003 344.014      0    1.094    1.106
    ## 27    PopAnh ~1             0.809 0.003 278.380      0    0.804    0.815
    ## 28  PopGuilt ~1             0.457 0.006  70.545      0    0.445    0.470
    ## 29   PopConc ~1             1.040 0.009 118.306      0    1.022    1.057
    ## 30    PopSui ~1             0.721 0.007 103.226      0    0.707    0.735
    ## 31 PopAppDec ~1             0.622 0.007  91.723      0    0.609    0.635
    ## 32 PopAppInc ~1             0.287 0.006  46.613      0    0.275    0.299
    ## 33 PopSleDec ~1             1.246 0.009 136.849      0    1.228    1.264
    ## 34 PopSleInc ~1             0.128 0.006  21.143      0    0.116    0.139
    ## 35  PopFatig ~1             1.313 0.010 134.952      0    1.294    1.332
    ## 36        A1 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 37        A2 ~1             0.000 0.000      NA     NA    0.000    0.000

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

standardizedSolution(pop_psych_soma_bif.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1          A =~    PopDep   0.581 0.026  21.949  0.000    0.529    0.633
    ## 2          A =~    PopAnh   0.668 0.027  24.833  0.000    0.616    0.721
    ## 3          A =~ PopAppDec   0.254 0.016  15.636  0.000    0.222    0.286
    ## 4          A =~ PopAppInc   0.285 0.019  14.875  0.000    0.248    0.323
    ## 5          A =~ PopSleDec   0.425 0.024  18.061  0.000    0.379    0.472
    ## 6          A =~ PopSleInc   0.307 0.022  13.726  0.000    0.263    0.350
    ## 7          A =~  PopFatig   0.655 0.032  20.708  0.000    0.593    0.716
    ## 8          A =~  PopGuilt   0.548 0.023  24.211  0.000    0.504    0.593
    ## 9          A =~   PopConc   0.694 0.030  23.389  0.000    0.636    0.752
    ## 10         A =~    PopSui   0.348 0.019  18.660  0.000    0.311    0.385
    ## 11        A1 =~    PopDep   0.630 0.025  25.221  0.000    0.581    0.679
    ## 12        A1 =~    PopAnh   0.424 0.026  16.261  0.000    0.373    0.476
    ## 13        A1 =~  PopGuilt   0.282 0.022  12.763  0.000    0.238    0.325
    ## 14        A1 =~   PopConc  -0.018 0.029  -0.613  0.540   -0.076    0.040
    ## 15        A1 =~    PopSui   0.215 0.024   9.044  0.000    0.168    0.261
    ## 16        A2 =~ PopAppDec   0.529 0.023  22.922  0.000    0.484    0.574
    ## 17        A2 =~ PopAppInc  -0.409 0.029 -14.333  0.000   -0.464   -0.353
    ## 18        A2 =~ PopSleDec   0.158 0.029   5.517  0.000    0.102    0.214
    ## 19        A2 =~ PopSleInc  -0.180 0.034  -5.319  0.000   -0.246   -0.114
    ## 20        A2 =~  PopFatig  -0.050 0.031  -1.584  0.113   -0.111    0.012
    ## 21        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 22        A2 ~~        A2   1.000 0.000      NA     NA    1.000    1.000
    ## 23         A ~~        A1   0.000 0.000      NA     NA    0.000    0.000
    ## 24         A ~~        A2   0.000 0.000      NA     NA    0.000    0.000
    ## 25        A1 ~~        A2   0.000 0.000      NA     NA    0.000    0.000
    ## 26    PopDep ~~    PopDep   0.266 0.040   6.668  0.000    0.188    0.345
    ## 27    PopAnh ~~    PopAnh   0.373 0.038   9.700  0.000    0.298    0.448
    ## 28 PopAppDec ~~ PopAppDec   0.656 0.026  25.551  0.000    0.606    0.706
    ## 29 PopAppInc ~~ PopAppInc   0.752 0.026  28.824  0.000    0.701    0.803
    ## 30 PopSleDec ~~ PopSleDec   0.794 0.022  36.349  0.000    0.751    0.837
    ## 31 PopSleInc ~~ PopSleInc   0.874 0.019  46.760  0.000    0.837    0.910
    ## 32  PopFatig ~~  PopFatig   0.569 0.041  13.722  0.000    0.488    0.650
    ## 33  PopGuilt ~~  PopGuilt   0.620 0.026  24.237  0.000    0.570    0.670
    ## 34   PopConc ~~   PopConc   0.518 0.041  12.574  0.000    0.437    0.599
    ## 35    PopSui ~~    PopSui   0.833 0.016  51.971  0.000    0.801    0.864
    ## 36         A ~~         A   1.000 0.000      NA     NA    1.000    1.000
    ## 37    PopDep ~1             1.100 0.046  23.750  0.000    1.009    1.191
    ## 38    PopAnh ~1             0.809 0.040  20.430  0.000    0.732    0.887
    ## 39 PopAppDec ~1             0.715 0.033  21.546  0.000    0.650    0.780
    ## 40 PopAppInc ~1             0.407 0.033  12.251  0.000    0.342    0.472
    ## 41 PopSleDec ~1             1.462 0.058  25.250  0.000    1.348    1.575
    ## 42 PopSleInc ~1             0.275 0.036   7.613  0.000    0.204    0.346
    ## 43  PopFatig ~1             1.644 0.079  20.787  0.000    1.489    1.799
    ## 44  PopGuilt ~1             0.471 0.026  18.258  0.000    0.420    0.521
    ## 45   PopConc ~1             1.436 0.067  21.563  0.000    1.306    1.567
    ## 46    PopSui ~1             0.696 0.030  23.424  0.000    0.638    0.754
    ## 47         A ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 48        A1 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 49        A2 ~1             0.000 0.000      NA     NA    0.000    0.000

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

standardizedSolution(pop_psych_veg.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~    PopDep   0.800 0.002 509.319      0    0.797    0.803
    ## 2         A1 =~    PopAnh   0.813 0.001 554.241      0    0.810    0.816
    ## 3         A1 =~  PopGuilt   0.618 0.004 142.266      0    0.609    0.626
    ## 4         A1 =~    PopSui   0.383 0.006  68.226      0    0.372    0.394
    ## 5         A2 =~ PopAppDec   0.285 0.006  51.617      0    0.274    0.296
    ## 6         A2 =~ PopAppInc   0.326 0.005  60.471      0    0.316    0.337
    ## 7         A2 =~ PopSleDec   0.498 0.005 105.359      0    0.488    0.507
    ## 8         A2 =~ PopSleInc   0.375 0.005  70.158      0    0.365    0.386
    ## 9         A2 =~  PopFatig   0.696 0.003 205.015      0    0.689    0.702
    ## 10        A2 =~   PopConc   0.736 0.003 229.785      0    0.729    0.742
    ## 11        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 12        A2 ~~        A2   1.000 0.000      NA     NA    1.000    1.000
    ## 13 PopAppDec ~~ PopAppInc  -0.244 0.004 -68.329      0   -0.251   -0.237
    ## 14 PopSleDec ~~ PopSleInc  -0.165 0.004 -42.301      0   -0.172   -0.157
    ## 15    PopDep ~~    PopDep   0.360 0.003 143.439      0    0.355    0.365
    ## 16    PopAnh ~~    PopAnh   0.338 0.002 141.770      0    0.334    0.343
    ## 17  PopGuilt ~~  PopGuilt   0.618 0.005 115.314      0    0.608    0.629
    ## 18    PopSui ~~    PopSui   0.853 0.004 198.441      0    0.845    0.862
    ## 19 PopAppDec ~~ PopAppDec   0.919 0.003 291.827      0    0.913    0.925
    ## 20 PopAppInc ~~ PopAppInc   0.893 0.004 253.535      0    0.887    0.900
    ## 21 PopSleDec ~~ PopSleDec   0.752 0.005 159.994      0    0.743    0.762
    ## 22 PopSleInc ~~ PopSleInc   0.859 0.004 214.300      0    0.851    0.867
    ## 23  PopFatig ~~  PopFatig   0.516 0.005 109.220      0    0.507    0.525
    ## 24   PopConc ~~   PopConc   0.459 0.005  97.371      0    0.449    0.468
    ## 25        A1 ~~        A2   0.875 0.003 286.206      0    0.869    0.881
    ## 26    PopDep ~1             1.100 0.003 344.017      0    1.094    1.106
    ## 27    PopAnh ~1             0.809 0.003 278.386      0    0.804    0.815
    ## 28  PopGuilt ~1             0.470 0.007  68.035      0    0.457    0.484
    ## 29    PopSui ~1             0.730 0.007 103.211      0    0.716    0.744
    ## 30 PopAppDec ~1             0.664 0.006 108.303      0    0.652    0.676
    ## 31 PopAppInc ~1             0.346 0.006  60.382      0    0.335    0.358
    ## 32 PopSleDec ~1             1.323 0.008 155.678      0    1.306    1.340
    ## 33 PopSleInc ~1             0.194 0.006  33.822      0    0.182    0.205
    ## 34  PopFatig ~1             1.447 0.010 148.174      0    1.428    1.466
    ## 35   PopConc ~1             1.215 0.009 130.197      0    1.196    1.233
    ## 36        A1 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 37        A2 ~1             0.000 0.000      NA     NA    0.000    0.000

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

standardizedSolution(pop_psych_veg_bif.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~    PopDep   0.651 0.008  79.023  0.000    0.635    0.668
    ## 2         A1 =~    PopAnh   0.434 0.007  65.392  0.000    0.421    0.447
    ## 3         A1 =~  PopGuilt   0.282 0.008  35.645  0.000    0.267    0.298
    ## 4         A1 =~    PopSui   0.225 0.009  25.495  0.000    0.208    0.242
    ## 5         A2 =~ PopAppDec   0.489 0.008  58.454  0.000    0.473    0.506
    ## 6         A2 =~ PopAppInc  -0.413 0.008 -51.666  0.000   -0.428   -0.397
    ## 7         A2 =~ PopSleDec   0.210 0.007  30.238  0.000    0.197    0.224
    ## 8         A2 =~ PopSleInc  -0.179 0.006 -28.219  0.000   -0.191   -0.166
    ## 9         A2 =~  PopFatig   0.015 0.006   2.375  0.018    0.003    0.028
    ## 10        A2 =~   PopConc   0.121 0.007  17.518  0.000    0.108    0.135
    ## 11         A =~    PopDep   0.564 0.007  79.864  0.000    0.550    0.578
    ## 12         A =~    PopAnh   0.660 0.004 157.798  0.000    0.652    0.669
    ## 13         A =~  PopGuilt   0.549 0.005 120.887  0.000    0.540    0.558
    ## 14         A =~    PopSui   0.344 0.005  63.048  0.000    0.333    0.355
    ## 15         A =~ PopAppDec   0.189 0.006  29.443  0.000    0.176    0.202
    ## 16         A =~ PopAppInc   0.332 0.006  57.763  0.000    0.321    0.343
    ## 17         A =~ PopSleDec   0.403 0.005  78.378  0.000    0.393    0.413
    ## 18         A =~ PopSleInc   0.324 0.005  63.057  0.000    0.314    0.334
    ## 19         A =~  PopFatig   0.651 0.004 166.508  0.000    0.643    0.658
    ## 20         A =~   PopConc   0.686 0.004 174.263  0.000    0.678    0.694
    ## 21        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 22        A2 ~~        A2   1.000 0.000      NA     NA    1.000    1.000
    ## 23         A ~~         A   1.000 0.000      NA     NA    1.000    1.000
    ## 24        A1 ~~         A   0.000 0.000      NA     NA    0.000    0.000
    ## 25        A2 ~~         A   0.000 0.000      NA     NA    0.000    0.000
    ## 26        A1 ~~        A2   0.000 0.000      NA     NA    0.000    0.000
    ## 27    PopDep ~~    PopDep   0.257 0.006  44.652  0.000    0.246    0.268
    ## 28    PopAnh ~~    PopAnh   0.375 0.003 135.874  0.000    0.370    0.381
    ## 29  PopGuilt ~~  PopGuilt   0.619 0.006  95.450  0.000    0.606    0.632
    ## 30    PopSui ~~    PopSui   0.831 0.006 131.634  0.000    0.819    0.843
    ## 31 PopAppDec ~~ PopAppDec   0.725 0.008  85.439  0.000    0.708    0.741
    ## 32 PopAppInc ~~ PopAppInc   0.719 0.007  95.973  0.000    0.705    0.734
    ## 33 PopSleDec ~~ PopSleDec   0.793 0.004 179.473  0.000    0.785    0.802
    ## 34 PopSleInc ~~ PopSleInc   0.863 0.004 214.312  0.000    0.855    0.871
    ## 35  PopFatig ~~  PopFatig   0.576 0.005 113.838  0.000    0.566    0.586
    ## 36   PopConc ~~   PopConc   0.515 0.005  97.225  0.000    0.504    0.525
    ## 37    PopDep ~1             1.100 0.003 344.018  0.000    1.094    1.106
    ## 38    PopAnh ~1             0.809 0.003 278.388  0.000    0.804    0.815
    ## 39  PopGuilt ~1             0.474 0.009  50.201  0.000    0.456    0.493
    ## 40    PopSui ~1             0.692 0.010  67.214  0.000    0.672    0.712
    ## 41 PopAppDec ~1             0.756 0.006 127.772  0.000    0.744    0.768
    ## 42 PopAppInc ~1             0.385 0.006  66.053  0.000    0.374    0.396
    ## 43 PopSleDec ~1             1.486 0.008 174.932  0.000    1.470    1.503
    ## 44 PopSleInc ~1             0.270 0.006  48.480  0.000    0.259    0.281
    ## 45  PopFatig ~1             1.666 0.012 139.993  0.000    1.643    1.689
    ## 46   PopConc ~1             1.444 0.012 121.654  0.000    1.421    1.467
    ## 47        A1 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 48        A2 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 49         A ~1             0.000 0.000      NA     NA    0.000    0.000

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

standardizedSolution(pop_affect_veg.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~    PopDep   0.827 0.002 454.756      0    0.824    0.831
    ## 2         A1 =~  PopGuilt   0.666 0.005 136.026      0    0.656    0.676
    ## 3         A1 =~    PopSui   0.430 0.006  68.467      0    0.418    0.443
    ## 4         A2 =~    PopAnh   0.793 0.001 619.642      0    0.790    0.795
    ## 5         A2 =~ PopAppInc   0.369 0.006  64.197      0    0.358    0.381
    ## 6         A2 =~ PopAppDec   0.318 0.006  52.591      0    0.306    0.330
    ## 7         A2 =~ PopSleInc   0.426 0.005  77.509      0    0.415    0.437
    ## 8         A2 =~ PopSleDec   0.515 0.005 102.522      0    0.505    0.525
    ## 9         A2 =~  PopFatig   0.711 0.003 217.859      0    0.705    0.718
    ## 10        A2 =~   PopConc   0.745 0.003 244.470      0    0.739    0.751
    ## 11        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 12        A2 ~~        A2   1.000 0.000      NA     NA    1.000    1.000
    ## 13 PopAppInc ~~ PopAppDec  -0.240 0.004 -67.635      0   -0.247   -0.233
    ## 14 PopSleInc ~~ PopSleDec  -0.148 0.004 -39.139      0   -0.156   -0.141
    ## 15    PopDep ~~    PopDep   0.316 0.003 104.906      0    0.310    0.322
    ## 16  PopGuilt ~~  PopGuilt   0.556 0.007  85.339      0    0.544    0.569
    ## 17    PopSui ~~    PopSui   0.815 0.005 150.609      0    0.804    0.825
    ## 18    PopAnh ~~    PopAnh   0.371 0.002 183.039      0    0.367    0.375
    ## 19 PopAppInc ~~ PopAppInc   0.864 0.004 203.257      0    0.855    0.872
    ## 20 PopAppDec ~~ PopAppDec   0.899 0.004 234.101      0    0.891    0.907
    ## 21 PopSleInc ~~ PopSleInc   0.819 0.005 174.784      0    0.809    0.828
    ## 22 PopSleDec ~~ PopSleDec   0.735 0.005 142.125      0    0.725    0.745
    ## 23  PopFatig ~~  PopFatig   0.494 0.005 106.392      0    0.485    0.503
    ## 24   PopConc ~~   PopConc   0.445 0.005  98.121      0    0.436    0.454
    ## 25        A1 ~~        A2   0.972 0.002 506.580      0    0.968    0.976
    ## 26    PopDep ~1             1.100 0.003 344.013      0    1.094    1.106
    ## 27  PopGuilt ~1             0.384 0.008  46.444      0    0.368    0.400
    ## 28    PopSui ~1             0.669 0.008  81.733      0    0.653    0.685
    ## 29    PopAnh ~1             0.809 0.003 278.378      0    0.804    0.815
    ## 30 PopAppInc ~1             0.274 0.006  43.237      0    0.261    0.286
    ## 31 PopAppDec ~1             0.603 0.007  86.417      0    0.589    0.616
    ## 32 PopSleInc ~1             0.110 0.006  18.057      0    0.098    0.122
    ## 33 PopSleDec ~1             1.218 0.009 132.313      0    1.200    1.236
    ## 34  PopFatig ~1             1.261 0.009 135.915      0    1.243    1.279
    ## 35   PopConc ~1             1.032 0.009 119.792      0    1.015    1.049
    ## 36        A1 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 37        A2 ~1             0.000 0.000      NA     NA    0.000    0.000

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

standardizedSolution(pop_affect_veg_bif.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~    PopDep   0.103 0.012   8.328  0.000    0.079    0.127
    ## 2         A1 =~  PopGuilt   0.131 0.015   8.764  0.000    0.102    0.161
    ## 3         A1 =~    PopSui   0.323 0.036   9.008  0.000    0.252    0.393
    ## 4         A2 =~    PopAnh  -0.001 0.004  -0.121  0.904   -0.009    0.008
    ## 5         A2 =~ PopAppInc  -0.363 0.008 -48.194  0.000   -0.378   -0.349
    ## 6         A2 =~ PopAppDec   0.522 0.010  54.903  0.000    0.503    0.540
    ## 7         A2 =~ PopSleInc  -0.162 0.007 -24.816  0.000   -0.175   -0.150
    ## 8         A2 =~ PopSleDec   0.199 0.007  27.513  0.000    0.185    0.213
    ## 9         A2 =~  PopFatig   0.014 0.007   2.125  0.034    0.001    0.028
    ## 10        A2 =~   PopConc   0.093 0.007  13.217  0.000    0.079    0.107
    ## 11         A =~    PopDep   0.803 0.001 564.425  0.000    0.801    0.806
    ## 12         A =~  PopGuilt   0.628 0.005 128.129  0.000    0.618    0.638
    ## 13         A =~    PopSui   0.421 0.007  59.569  0.000    0.407    0.434
    ## 14         A =~    PopAnh   0.795 0.001 613.477  0.000    0.793    0.798
    ## 15         A =~ PopAppInc   0.393 0.007  54.569  0.000    0.379    0.407
    ## 16         A =~ PopAppDec   0.265 0.009  28.176  0.000    0.247    0.283
    ## 17         A =~ PopSleInc   0.416 0.006  70.755  0.000    0.404    0.427
    ## 18         A =~ PopSleDec   0.481 0.006  84.140  0.000    0.469    0.492
    ## 19         A =~  PopFatig   0.711 0.003 213.476  0.000    0.704    0.717
    ## 20         A =~   PopConc   0.741 0.003 229.042  0.000    0.734    0.747
    ## 21        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 22        A2 ~~        A2   1.000 0.000      NA     NA    1.000    1.000
    ## 23         A ~~         A   1.000 0.000      NA     NA    1.000    1.000
    ## 24        A1 ~~         A   0.000 0.000      NA     NA    0.000    0.000
    ## 25        A2 ~~         A   0.000 0.000      NA     NA    0.000    0.000
    ## 26        A1 ~~        A2   0.000 0.000      NA     NA    0.000    0.000
    ## 27    PopDep ~~    PopDep   0.344 0.003 113.332  0.000    0.338    0.350
    ## 28  PopGuilt ~~  PopGuilt   0.588 0.009  66.873  0.000    0.571    0.606
    ## 29    PopSui ~~    PopSui   0.719 0.025  28.919  0.000    0.670    0.768
    ## 30    PopAnh ~~    PopAnh   0.367 0.002 178.073  0.000    0.363    0.371
    ## 31 PopAppInc ~~ PopAppInc   0.713 0.007 105.246  0.000    0.700    0.727
    ## 32 PopAppDec ~~ PopAppDec   0.658 0.010  63.009  0.000    0.637    0.678
    ## 33 PopSleInc ~~ PopSleInc   0.801 0.005 155.485  0.000    0.791    0.811
    ## 34 PopSleDec ~~ PopSleDec   0.729 0.005 137.661  0.000    0.719    0.740
    ## 35  PopFatig ~~  PopFatig   0.495 0.005 104.969  0.000    0.486    0.504
    ## 36   PopConc ~~   PopConc   0.443 0.005  95.540  0.000    0.434    0.452
    ## 37    PopDep ~1             1.100 0.003 344.014  0.000    1.094    1.106
    ## 38  PopGuilt ~1             0.432 0.010  44.788  0.000    0.413    0.451
    ## 39    PopSui ~1             0.652 0.011  60.045  0.000    0.631    0.673
    ## 40    PopAnh ~1             0.809 0.003 278.379  0.000    0.804    0.815
    ## 41 PopAppInc ~1             0.254 0.007  34.709  0.000    0.240    0.268
    ## 42 PopAppDec ~1             0.649 0.009  71.580  0.000    0.631    0.666
    ## 43 PopSleInc ~1             0.120 0.006  19.075  0.000    0.108    0.133
    ## 44 PopSleDec ~1             1.262 0.010 131.331  0.000    1.243    1.281
    ## 45  PopFatig ~1             1.265 0.009 134.947  0.000    1.246    1.283
    ## 46   PopConc ~1             1.041 0.009 117.897  0.000    1.024    1.059
    ## 47        A1 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 48        A2 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 49         A ~1             0.000 0.000      NA     NA    0.000    0.000

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
standardizedSolution(pop_cog_mood_neuroveg.fit)
```

    ##          lhs op       rhs est.std    se       z pvalue ci.lower ci.upper
    ## 1         A1 =~  PopGuilt   0.312 0.010  30.116      0    0.291    0.332
    ## 2         A1 =~   PopConc   0.690 0.004 165.382      0    0.681    0.698
    ## 3         A1 =~    PopSui   0.305 0.005  60.764      0    0.295    0.315
    ## 4         A2 =~    PopDep   0.775 0.002 341.975      0    0.771    0.780
    ## 5         A2 =~    PopAnh   0.846 0.002 361.818      0    0.841    0.851
    ## 6         A2 =~  PopGuilt   0.252 0.012  20.478      0    0.228    0.276
    ## 7         A3 =~ PopSleDec   0.461 0.005  95.532      0    0.452    0.470
    ## 8         A3 =~ PopSleInc   0.353 0.005  68.941      0    0.343    0.363
    ## 9         A3 =~  PopFatig   0.660 0.004 172.764      0    0.652    0.667
    ## 10        A3 =~ PopAppDec   0.262 0.005  49.868      0    0.251    0.272
    ## 11        A3 =~ PopAppInc   0.311 0.005  60.983      0    0.301    0.321
    ## 12        A1 ~~        A1   1.000 0.000      NA     NA    1.000    1.000
    ## 13        A2 ~~        A2   1.000 0.000      NA     NA    1.000    1.000
    ## 14        A3 ~~        A3   1.000 0.000      NA     NA    1.000    1.000
    ## 15 PopSleDec ~~ PopSleInc  -0.166 0.004 -41.998      0   -0.173   -0.158
    ## 16 PopAppDec ~~ PopAppInc  -0.246 0.004 -68.704      0   -0.253   -0.239
    ## 17  PopGuilt ~~  PopGuilt   0.708 0.005 131.638      0    0.697    0.718
    ## 18   PopConc ~~   PopConc   0.525 0.006  91.223      0    0.513    0.536
    ## 19    PopSui ~~    PopSui   0.907 0.003 296.572      0    0.901    0.913
    ## 20    PopDep ~~    PopDep   0.399 0.004 113.355      0    0.392    0.406
    ## 21    PopAnh ~~    PopAnh   0.284 0.004  71.864      0    0.277    0.292
    ## 22 PopSleDec ~~ PopSleDec   0.787 0.004 176.973      0    0.779    0.796
    ## 23 PopSleInc ~~ PopSleInc   0.875 0.004 242.034      0    0.868    0.882
    ## 24  PopFatig ~~  PopFatig   0.565 0.005 112.142      0    0.555    0.575
    ## 25 PopAppDec ~~ PopAppDec   0.932 0.003 339.398      0    0.926    0.937
    ## 26 PopAppInc ~~ PopAppInc   0.903 0.003 284.753      0    0.897    0.909
    ## 27        A1 ~~        A2   0.840 0.007 119.498      0    0.826    0.854
    ## 28        A1 ~~        A3   1.013 0.005 185.055      0    1.002    1.024
    ## 29        A2 ~~        A3   0.802 0.006 141.638      0    0.791    0.813
    ## 30  PopGuilt ~1             0.600 0.008  71.148      0    0.584    0.617
    ## 31   PopConc ~1             1.339 0.011 117.641      0    1.317    1.362
    ## 32    PopSui ~1             0.851 0.006 133.199      0    0.839    0.864
    ## 33    PopDep ~1             1.100 0.003 344.019      0    1.094    1.106
    ## 34    PopAnh ~1             0.809 0.003 278.388      0    0.804    0.815
    ## 35 PopSleDec ~1             1.410 0.009 162.960      0    1.393    1.427
    ## 36 PopSleInc ~1             0.235 0.006  41.677      0    0.224    0.246
    ## 37  PopFatig ~1             1.592 0.011 143.611      0    1.570    1.614
    ## 38 PopAppDec ~1             0.702 0.006 119.637      0    0.690    0.713
    ## 39 PopAppInc ~1             0.380 0.006  68.664      0    0.369    0.391
    ## 40        A1 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 41        A2 ~1             0.000 0.000      NA     NA    0.000    0.000
    ## 42        A3 ~1             0.000 0.000      NA     NA    0.000    0.000

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

standardizedSolution(pop_cog_mood_neuroveg_bif.fit)
```

    ##          lhs op       rhs est.std    se      z pvalue ci.lower ci.upper
    ## 1         A1 =~  PopGuilt   0.067 0.059  1.121  0.262   -0.050    0.183
    ## 2         A1 =~   PopConc  -0.042 0.064 -0.660  0.509   -0.168    0.083
    ## 3         A1 =~    PopSui   0.709 0.393  1.805  0.071   -0.061    1.479
    ## 4         A2 =~    PopDep   0.586 0.052 11.372  0.000    0.485    0.687
    ## 5         A2 =~    PopAnh   0.431 0.054  8.001  0.000    0.326    0.537
    ## 6         A2 =~  PopGuilt   0.242 0.047  5.186  0.000    0.150    0.333
    ## 7         A3 =~ PopSleDec   0.155 0.059  2.636  0.008    0.040    0.270
    ## 8         A3 =~ PopSleInc  -0.181 0.069 -2.613  0.009   -0.317   -0.045
    ## 9         A3 =~  PopFatig  -0.052 0.064 -0.818  0.413   -0.178    0.073
    ## 10        A3 =~ PopAppDec   0.526 0.047 11.103  0.000    0.433    0.619
    ## 11        A3 =~ PopAppInc  -0.412 0.058 -7.044  0.000   -0.526   -0.297
    ## 12         A =~    PopDep   0.596 0.055 10.893  0.000    0.489    0.703
    ## 13         A =~  PopGuilt   0.535 0.047 11.389  0.000    0.443    0.628
    ## 14         A =~    PopSui   0.298 0.037  7.968  0.000    0.224    0.371
    ## 15         A =~    PopAnh   0.677 0.056 12.131  0.000    0.567    0.786
    ## 16         A =~ PopAppInc   0.287 0.040  7.236  0.000    0.209    0.365
    ## 17         A =~ PopAppDec   0.262 0.034  7.675  0.000    0.195    0.328
    ## 18         A =~ PopSleInc   0.309 0.046  6.687  0.000    0.218    0.399
    ## 19         A =~ PopSleDec   0.431 0.049  8.804  0.000    0.335    0.526
    ## 20         A =~  PopFatig   0.652 0.065 10.072  0.000    0.525    0.779
    ## 21         A =~   PopConc   0.703 0.060 11.618  0.000    0.584    0.821
    ## 22        A1 ~~        A1   1.000 0.000     NA     NA    1.000    1.000
    ## 23        A2 ~~        A2   1.000 0.000     NA     NA    1.000    1.000
    ## 24        A3 ~~        A3   1.000 0.000     NA     NA    1.000    1.000
    ## 25        A1 ~~         A   0.000 0.000     NA     NA    0.000    0.000
    ## 26        A2 ~~         A   0.000 0.000     NA     NA    0.000    0.000
    ## 27        A3 ~~         A   0.000 0.000     NA     NA    0.000    0.000
    ## 28        A1 ~~        A2   0.000 0.000     NA     NA    0.000    0.000
    ## 29        A1 ~~        A3   0.000 0.000     NA     NA    0.000    0.000
    ## 30        A2 ~~        A3   0.000 0.000     NA     NA    0.000    0.000
    ## 31  PopGuilt ~~  PopGuilt   0.650 0.052 12.432  0.000    0.548    0.753
    ## 32   PopConc ~~   PopConc   0.504 0.085  5.925  0.000    0.338    0.671
    ## 33    PopSui ~~    PopSui   0.409 0.558  0.733  0.464   -0.685    1.502
    ## 34    PopDep ~~    PopDep   0.301 0.080  3.774  0.000    0.145    0.457
    ## 35    PopAnh ~~    PopAnh   0.356 0.080  4.435  0.000    0.199    0.514
    ## 36 PopSleDec ~~ PopSleDec   0.791 0.046 17.367  0.000    0.701    0.880
    ## 37 PopSleInc ~~ PopSleInc   0.872 0.039 22.487  0.000    0.796    0.948
    ## 38  PopFatig ~~  PopFatig   0.572 0.085  6.746  0.000    0.406    0.738
    ## 39 PopAppDec ~~ PopAppDec   0.655 0.053 12.440  0.000    0.552    0.758
    ## 40 PopAppInc ~~ PopAppInc   0.748 0.054 13.864  0.000    0.642    0.854
    ## 41         A ~~         A   1.000 0.000     NA     NA    1.000    1.000
    ## 42  PopGuilt ~1             0.518 0.055  9.393  0.000    0.410    0.627
    ## 43   PopConc ~1             1.392 0.133 10.483  0.000    1.132    1.652
    ## 44    PopSui ~1             0.874 0.076 11.439  0.000    0.724    1.024
    ## 45    PopDep ~1             1.100 0.095 11.553  0.000    0.913    1.287
    ## 46    PopAnh ~1             0.809 0.082  9.865  0.000    0.648    0.970
    ## 47 PopSleDec ~1             1.450 0.118 12.296  0.000    1.219    1.681
    ## 48 PopSleInc ~1             0.270 0.074  3.644  0.000    0.125    0.415
    ## 49  PopFatig ~1             1.630 0.161 10.154  0.000    1.315    1.944
    ## 50 PopAppDec ~1             0.707 0.068 10.432  0.000    0.574    0.840
    ## 51 PopAppInc ~1             0.402 0.068  5.910  0.000    0.269    0.535
    ## 52        A1 ~1             0.000 0.000     NA     NA    0.000    0.000
    ## 53        A2 ~1             0.000 0.000     NA     NA    0.000    0.000
    ## 54        A3 ~1             0.000 0.000     NA     NA    0.000    0.000
    ## 55         A ~1             0.000 0.000     NA     NA    0.000    0.000

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

## Categorical models

Treat symptoms as ordered categorical variables

``` {.r}
pop_commonfactor_ord.fit <- cfa(model=pop_commonfactor.model, data=ukb_symp, ordered=TRUE, missing='pairwise')
```

    ## Warning in lav_object_post_check(object): lavaan WARNING: some estimated ov
    ## variances are negative

``` {.r}
summary(pop_commonfactor_ord.fit, fit.measures=TRUE)
```

    ## lavaan 0.6-7 ended normally after 15 iterations
    ## 
    ##   Estimator                                       DWLS
    ##   Optimization method                           NLMINB
    ##   Number of free parameters                         20
    ##                                                       
    ##   Number of observations                        157196
    ##   Number of missing patterns                       136
    ##                                                       
    ## Model Test User Model:
    ##                                                Standard      Robust
    ##   Test Statistic                              61711.408   59905.839
    ##   Degrees of freedom                                 35          35
    ##   P-value (Chi-square)                            0.000       0.000
    ##   Scaling correction factor                                   1.030
    ##   Shift parameter                                             2.302
    ##        simple second-order correction                              
    ## 
    ## Model Test Baseline Model:
    ## 
    ##   Test statistic                            587591.796  516633.465
    ##   Degrees of freedom                                45          45
    ##   P-value                                        0.000       0.000
    ##   Scaling correction factor                                  1.137
    ## 
    ## User Model versus Baseline Model:
    ## 
    ##   Comparative Fit Index (CFI)                    0.895       0.884
    ##   Tucker-Lewis Index (TLI)                       0.865       0.851
    ##                                                                   
    ##   Robust Comparative Fit Index (CFI)                            NA
    ##   Robust Tucker-Lewis Index (TLI)                               NA
    ## 
    ## Root Mean Square Error of Approximation:
    ## 
    ##   RMSEA                                          0.106       0.104
    ##   90 Percent confidence interval - lower         0.105       0.104
    ##   90 Percent confidence interval - upper         0.107       0.105
    ##   P-value RMSEA <= 0.05                          0.000       0.000
    ##                                                                   
    ##   Robust RMSEA                                                  NA
    ##   90 Percent confidence interval - lower                        NA
    ##   90 Percent confidence interval - upper                        NA
    ## 
    ## Standardized Root Mean Square Residual:
    ## 
    ##   SRMR                                           0.200       0.200
    ## 
    ## Parameter Estimates:
    ## 
    ##   Standard errors                           Robust.sem
    ##   Information                                 Expected
    ##   Information saturated (h1) model        Unstructured
    ## 
    ## Latent Variables:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##   A1 =~                                               
    ##     PopDep            0.486    0.004  118.050    0.000
    ##     PopAnh            1.833    0.016  115.618    0.000
    ##     PopAppDec         0.088    0.002   39.492    0.000
    ##     PopAppInc         0.107    0.002   44.089    0.000
    ##     PopSleDec         0.137    0.002   56.303    0.000
    ##     PopSleInc         0.153    0.003   59.775    0.000
    ##     PopFatig          0.234    0.002   96.639    0.000
    ##     PopGuilt          0.200    0.002   90.744    0.000
    ##     PopConc           0.242    0.002   99.752    0.000
    ##     PopSui            0.117    0.002   54.203    0.000
    ## 
    ## Intercepts:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##    .PopDep            0.000                           
    ##    .PopAnh            0.000                           
    ##    .PopAppDec         0.000                           
    ##    .PopAppInc         0.000                           
    ##    .PopSleDec         0.000                           
    ##    .PopSleInc         0.000                           
    ##    .PopFatig          0.000                           
    ##    .PopGuilt          0.000                           
    ##    .PopConc           0.000                           
    ##    .PopSui            0.000                           
    ##     A1                0.000                           
    ## 
    ## Thresholds:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     PopDep|t1        -0.120    0.003  -37.752    0.000
    ##     PopAnh|t1         0.265    0.003   82.799    0.000
    ##     PopAppDec|t1      0.195    0.005   42.504    0.000
    ##     PopAppInc|t1      0.708    0.005  141.690    0.000
    ##     PopSleDec|t1     -0.678    0.005 -137.119    0.000
    ##     PopSleInc|t1      0.967    0.005  178.945    0.000
    ##     PopFatig|t1      -0.912    0.005 -176.127    0.000
    ##     PopGuilt|t1      -0.024    0.004   -5.556    0.000
    ##     PopConc|t1       -0.795    0.005 -158.068    0.000
    ##     PopSui|t1        -0.052    0.004  -12.105    0.000
    ## 
    ## Variances:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     A1                1.000                           
    ##    .PopDep            0.764                           
    ##    .PopAnh           -2.359                           
    ##    .PopAppDec         0.992                           
    ##    .PopAppInc         0.988                           
    ##    .PopSleDec         0.981                           
    ##    .PopSleInc         0.977                           
    ##    .PopFatig          0.945                           
    ##    .PopGuilt          0.960                           
    ##    .PopConc           0.941                           
    ##    .PopSui            0.986                           
    ## 
    ## Scales y*:
    ##                    Estimate  Std.Err  z-value  P(>|z|)
    ##     PopDep            1.000                           
    ##     PopAnh            1.000                           
    ##     PopAppDec         1.000                           
    ##     PopAppInc         1.000                           
    ##     PopSleDec         1.000                           
    ##     PopSleInc         1.000                           
    ##     PopFatig          1.000                           
    ##     PopGuilt          1.000                           
    ##     PopConc           1.000                           
    ##     PopSui            1.000
