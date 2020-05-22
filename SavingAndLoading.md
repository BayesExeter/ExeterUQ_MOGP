Saving and loading
------------------

`ExeterUQ_mogp` emulators are not standard objects. They are lists with
2 elements: the first of which is a python mogp\_emulator object, and
the second is a list of statistical elements of the fit such as prior
choices, mean function choices, elements useful for diagnostics and
other things we would like for transparent inference.

The usual `save()` and `load()` functions will seem to work, but will
only save a pointer for the python object (which once removed/reloaded
won’t work). Here is a MWE

``` r
mogp_dir <- "~/Dropbox/BayesExeter/mogp_emulator"
```

``` r
setwd('..')
source('BuildEmulator/BuildEmulator.R')
```

``` r
load("ConvectionModelExample.Rdata")
TestEm <- BuildNewEmulators(tData, HowManyEmulators = 2, meanFun="fitted")
```

    ## [1] "Max reduction is 0.141510566268957 using A_EPSILON"
    ## [1] "Max reduction is 0.0719279431890675 using A_U"
    ## [1] "Max reduction is 0.0593934956606317 using A_T"
    ## [1] "Max reduction is 0.0841015596079704 using A_EPSILON"
    ## [1] "Max reduction is 0.0215815680504046 using A_EPSILON"
    ## [1] "Max reduction is 0.0221974794878812 using A_EPSILON"
    ## [1] "Noise fitted, stopping algorithm"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_zav.400.600.theta_5_6 ~ A_EPSILON + 
    ##     I(A_EPSILON^2) + I(A_EPSILON^3) + I(A_EPSILON^4) + A_U + 
    ##     A_T + I(A_U * A_EPSILON) + I(A_T * A_EPSILON) + I(A_T * A_U), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15059 -0.06624  0.01615  0.05931  0.15216 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)        306.35982    0.03461 8852.437  < 2e-16 ***
    ## A_EPSILON           -0.31275    0.08180   -3.824  0.00106 ** 
    ## I(A_EPSILON^2)      -0.19500    0.21788   -0.895  0.38145    
    ## I(A_EPSILON^3)      -0.43915    0.12490   -3.516  0.00217 ** 
    ## I(A_EPSILON^4)       0.83401    0.24582    3.393  0.00289 ** 
    ## A_U                  0.30423    0.03234    9.409 8.73e-09 ***
    ## A_T                  0.33990    0.03172   10.715 9.78e-10 ***
    ## I(A_U * A_EPSILON)  -0.08042    0.05404   -1.488  0.15233    
    ## I(A_T * A_EPSILON)  -0.08069    0.06159   -1.310  0.20500    
    ## I(A_T * A_U)        -0.03824    0.06560   -0.583  0.56648    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09867 on 20 degrees of freedom
    ## Multiple R-squared:  0.9731, Adjusted R-squared:  0.961 
    ## F-statistic: 80.31 on 9 and 20 DF,  p-value: 1.004e-13
    ## 
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_zav.400.600.theta_5_6 ~ A_EPSILON + 
    ##     I(A_EPSILON^2) + I(A_EPSILON^3) + I(A_EPSILON^4) + A_U + 
    ##     A_T + I(A_U * A_EPSILON) + I(A_T * A_EPSILON) + I(A_T * A_U), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15059 -0.06624  0.01615  0.05931  0.15216 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)        306.35982    0.03461 8852.437  < 2e-16 ***
    ## A_EPSILON           -0.31275    0.08180   -3.824  0.00106 ** 
    ## I(A_EPSILON^2)      -0.19500    0.21788   -0.895  0.38145    
    ## I(A_EPSILON^3)      -0.43915    0.12490   -3.516  0.00217 ** 
    ## I(A_EPSILON^4)       0.83401    0.24582    3.393  0.00289 ** 
    ## A_U                  0.30423    0.03234    9.409 8.73e-09 ***
    ## A_T                  0.33990    0.03172   10.715 9.78e-10 ***
    ## I(A_U * A_EPSILON)  -0.08042    0.05404   -1.488  0.15233    
    ## I(A_T * A_EPSILON)  -0.08069    0.06159   -1.310  0.20500    
    ## I(A_T * A_U)        -0.03824    0.06560   -0.583  0.56648    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09867 on 20 degrees of freedom
    ## Multiple R-squared:  0.9731, Adjusted R-squared:  0.961 
    ## F-statistic: 80.31 on 9 and 20 DF,  p-value: 1.004e-13
    ## 
    ## [1] "1 I(A_T * A_U)"
    ## [1] "2 I(A_T * A_EPSILON)"
    ## [1] "3 I(A_U * A_EPSILON)"
    ## [1] "4 I(A_EPSILON^4)"
    ## [1] "5 I(A_EPSILON^3)"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_zav.400.600.theta_5_6 ~ A_EPSILON + 
    ##     I(A_EPSILON^2) + A_U + A_T, data = tData)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.31141 -0.07230  0.00132  0.06857  0.49418 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)    306.28286    0.04193 7304.864  < 2e-16 ***
    ## A_EPSILON       -0.57626    0.04796  -12.015 7.00e-12 ***
    ## I(A_EPSILON^2)   0.53877    0.09315    5.784 4.98e-06 ***
    ## A_U              0.32558    0.04794    6.792 4.06e-07 ***
    ## A_T              0.34208    0.04790    7.142 1.74e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.152 on 25 degrees of freedom
    ## Multiple R-squared:  0.9202, Adjusted R-squared:  0.9074 
    ## F-statistic: 72.06 on 4 and 25 DF,  p-value: 2.36e-13
    ## 
    ## [1] "Max reduction is 0.0516838339369975 using A_EPSILON"
    ## [1] "Max reduction is 0.0181512292708917 using A_T"
    ## [1] "Max reduction is 0.00908627733663028 using A_U"
    ## [1] "Noise fitted, stopping algorithm"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_Ay.theta_5_6 ~ A_EPSILON + A_T + 
    ##     A_U + I(A_T * A_EPSILON) + I(A_U * A_EPSILON) + I(A_U * A_T), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.122972 -0.017118  0.003162  0.026879  0.108563 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -0.466705   0.009108 -51.243  < 2e-16 ***
    ## A_EPSILON           0.173565   0.015924  10.900 1.47e-10 ***
    ## A_T                -0.078502   0.015793  -4.971 5.01e-05 ***
    ## A_U                -0.050460   0.015915  -3.171  0.00427 ** 
    ## I(A_T * A_EPSILON)  0.066547   0.030606   2.174  0.04022 *  
    ## I(A_U * A_EPSILON) -0.025424   0.025766  -0.987  0.33405    
    ## I(A_U * A_T)       -0.035683   0.031953  -1.117  0.27564    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0496 on 23 degrees of freedom
    ## Multiple R-squared:  0.8819, Adjusted R-squared:  0.8511 
    ## F-statistic: 28.62 on 6 and 23 DF,  p-value: 1.439e-09
    ## 
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_Ay.theta_5_6 ~ A_EPSILON + A_T + 
    ##     A_U + I(A_T * A_EPSILON) + I(A_U * A_EPSILON) + I(A_U * A_T), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.122972 -0.017118  0.003162  0.026879  0.108563 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -0.466705   0.009108 -51.243  < 2e-16 ***
    ## A_EPSILON           0.173565   0.015924  10.900 1.47e-10 ***
    ## A_T                -0.078502   0.015793  -4.971 5.01e-05 ***
    ## A_U                -0.050460   0.015915  -3.171  0.00427 ** 
    ## I(A_T * A_EPSILON)  0.066547   0.030606   2.174  0.04022 *  
    ## I(A_U * A_EPSILON) -0.025424   0.025766  -0.987  0.33405    
    ## I(A_U * A_T)       -0.035683   0.031953  -1.117  0.27564    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0496 on 23 degrees of freedom
    ## Multiple R-squared:  0.8819, Adjusted R-squared:  0.8511 
    ## F-statistic: 28.62 on 6 and 23 DF,  p-value: 1.439e-09
    ## 
    ## [1] "1 I(A_U * A_EPSILON)"
    ## [1] "2 I(A_U * A_T)"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_Ay.theta_5_6 ~ A_EPSILON + A_T + 
    ##     A_U + I(A_T * A_EPSILON), data = tData)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.125503 -0.019057  0.002386  0.026174  0.097011 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -0.466250   0.009099 -51.240  < 2e-16 ***
    ## A_EPSILON           0.175425   0.015790  11.110 3.68e-11 ***
    ## A_T                -0.076189   0.015710  -4.850 5.50e-05 ***
    ## A_U                -0.052910   0.015894  -3.329   0.0027 ** 
    ## I(A_T * A_EPSILON)  0.064706   0.030065   2.152   0.0412 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04981 on 25 degrees of freedom
    ## Multiple R-squared:  0.8705, Adjusted R-squared:  0.8498 
    ## F-statistic: 42.01 on 4 and 25 DF,  p-value: 9.526e-11

``` r
print(TestEm[[1]])
```

    ## Multi-Output Gaussian Process with:
    ## 2 emulators
    ## 30 training examples
    ## 3 input variables

``` r
save(TestEm, file="TestEm.RData")
rm(TestEm)
load("TestEm.RData")
print(TestEm[[1]])
```

    ## <pointer: 0x0>

The mogp part has gone on reload (it was never saved). To overcome this
we use the python functions `py_save_object` and `py_load_object` to
save out the mogp part. This means that our saving and loading of an
mogp emulator actually saves an RData file and a python object
separately. This is all handled automatically within our package using
the functions `saveExUQmogp` and `loadExUQmogp` First we build an
emulator again

``` r
TestEm <- BuildNewEmulators(tData, HowManyEmulators = 2, meanFun="fitted")
```

    ## [1] "Max reduction is 0.141510566268957 using A_EPSILON"
    ## [1] "Max reduction is 0.0719279431890675 using A_U"
    ## [1] "Max reduction is 0.0593934956606317 using A_T"
    ## [1] "Max reduction is 0.0841015596079704 using A_EPSILON"
    ## [1] "Max reduction is 0.0215815680504046 using A_EPSILON"
    ## [1] "Max reduction is 0.0221974794878812 using A_EPSILON"
    ## [1] "Noise fitted, stopping algorithm"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_zav.400.600.theta_5_6 ~ A_EPSILON + 
    ##     I(A_EPSILON^2) + I(A_EPSILON^3) + I(A_EPSILON^4) + A_U + 
    ##     A_T + I(A_U * A_EPSILON) + I(A_T * A_EPSILON) + I(A_T * A_U), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15059 -0.06624  0.01615  0.05931  0.15216 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)        306.35982    0.03461 8852.437  < 2e-16 ***
    ## A_EPSILON           -0.31275    0.08180   -3.824  0.00106 ** 
    ## I(A_EPSILON^2)      -0.19500    0.21788   -0.895  0.38145    
    ## I(A_EPSILON^3)      -0.43915    0.12490   -3.516  0.00217 ** 
    ## I(A_EPSILON^4)       0.83401    0.24582    3.393  0.00289 ** 
    ## A_U                  0.30423    0.03234    9.409 8.73e-09 ***
    ## A_T                  0.33990    0.03172   10.715 9.78e-10 ***
    ## I(A_U * A_EPSILON)  -0.08042    0.05404   -1.488  0.15233    
    ## I(A_T * A_EPSILON)  -0.08069    0.06159   -1.310  0.20500    
    ## I(A_T * A_U)        -0.03824    0.06560   -0.583  0.56648    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09867 on 20 degrees of freedom
    ## Multiple R-squared:  0.9731, Adjusted R-squared:  0.961 
    ## F-statistic: 80.31 on 9 and 20 DF,  p-value: 1.004e-13
    ## 
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_zav.400.600.theta_5_6 ~ A_EPSILON + 
    ##     I(A_EPSILON^2) + I(A_EPSILON^3) + I(A_EPSILON^4) + A_U + 
    ##     A_T + I(A_U * A_EPSILON) + I(A_T * A_EPSILON) + I(A_T * A_U), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.15059 -0.06624  0.01615  0.05931  0.15216 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)        306.35982    0.03461 8852.437  < 2e-16 ***
    ## A_EPSILON           -0.31275    0.08180   -3.824  0.00106 ** 
    ## I(A_EPSILON^2)      -0.19500    0.21788   -0.895  0.38145    
    ## I(A_EPSILON^3)      -0.43915    0.12490   -3.516  0.00217 ** 
    ## I(A_EPSILON^4)       0.83401    0.24582    3.393  0.00289 ** 
    ## A_U                  0.30423    0.03234    9.409 8.73e-09 ***
    ## A_T                  0.33990    0.03172   10.715 9.78e-10 ***
    ## I(A_U * A_EPSILON)  -0.08042    0.05404   -1.488  0.15233    
    ## I(A_T * A_EPSILON)  -0.08069    0.06159   -1.310  0.20500    
    ## I(A_T * A_U)        -0.03824    0.06560   -0.583  0.56648    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.09867 on 20 degrees of freedom
    ## Multiple R-squared:  0.9731, Adjusted R-squared:  0.961 
    ## F-statistic: 80.31 on 9 and 20 DF,  p-value: 1.004e-13
    ## 
    ## [1] "1 I(A_T * A_U)"
    ## [1] "2 I(A_T * A_EPSILON)"
    ## [1] "3 I(A_U * A_EPSILON)"
    ## [1] "4 I(A_EPSILON^4)"
    ## [1] "5 I(A_EPSILON^3)"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_zav.400.600.theta_5_6 ~ A_EPSILON + 
    ##     I(A_EPSILON^2) + A_U + A_T, data = tData)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.31141 -0.07230  0.00132  0.06857  0.49418 
    ## 
    ## Coefficients:
    ##                 Estimate Std. Error  t value Pr(>|t|)    
    ## (Intercept)    306.28286    0.04193 7304.864  < 2e-16 ***
    ## A_EPSILON       -0.57626    0.04796  -12.015 7.00e-12 ***
    ## I(A_EPSILON^2)   0.53877    0.09315    5.784 4.98e-06 ***
    ## A_U              0.32558    0.04794    6.792 4.06e-07 ***
    ## A_T              0.34208    0.04790    7.142 1.74e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.152 on 25 degrees of freedom
    ## Multiple R-squared:  0.9202, Adjusted R-squared:  0.9074 
    ## F-statistic: 72.06 on 4 and 25 DF,  p-value: 2.36e-13
    ## 
    ## [1] "Max reduction is 0.0516838339369975 using A_EPSILON"
    ## [1] "Max reduction is 0.0181512292708917 using A_T"
    ## [1] "Max reduction is 0.00908627733663028 using A_U"
    ## [1] "Noise fitted, stopping algorithm"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_Ay.theta_5_6 ~ A_EPSILON + A_T + 
    ##     A_U + I(A_T * A_EPSILON) + I(A_U * A_EPSILON) + I(A_U * A_T), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.122972 -0.017118  0.003162  0.026879  0.108563 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -0.466705   0.009108 -51.243  < 2e-16 ***
    ## A_EPSILON           0.173565   0.015924  10.900 1.47e-10 ***
    ## A_T                -0.078502   0.015793  -4.971 5.01e-05 ***
    ## A_U                -0.050460   0.015915  -3.171  0.00427 ** 
    ## I(A_T * A_EPSILON)  0.066547   0.030606   2.174  0.04022 *  
    ## I(A_U * A_EPSILON) -0.025424   0.025766  -0.987  0.33405    
    ## I(A_U * A_T)       -0.035683   0.031953  -1.117  0.27564    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0496 on 23 degrees of freedom
    ## Multiple R-squared:  0.8819, Adjusted R-squared:  0.8511 
    ## F-statistic: 28.62 on 6 and 23 DF,  p-value: 1.439e-09
    ## 
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_Ay.theta_5_6 ~ A_EPSILON + A_T + 
    ##     A_U + I(A_T * A_EPSILON) + I(A_U * A_EPSILON) + I(A_U * A_T), 
    ##     data = tData)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.122972 -0.017118  0.003162  0.026879  0.108563 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -0.466705   0.009108 -51.243  < 2e-16 ***
    ## A_EPSILON           0.173565   0.015924  10.900 1.47e-10 ***
    ## A_T                -0.078502   0.015793  -4.971 5.01e-05 ***
    ## A_U                -0.050460   0.015915  -3.171  0.00427 ** 
    ## I(A_T * A_EPSILON)  0.066547   0.030606   2.174  0.04022 *  
    ## I(A_U * A_EPSILON) -0.025424   0.025766  -0.987  0.33405    
    ## I(A_U * A_T)       -0.035683   0.031953  -1.117  0.27564    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0496 on 23 degrees of freedom
    ## Multiple R-squared:  0.8819, Adjusted R-squared:  0.8511 
    ## F-statistic: 28.62 on 6 and 23 DF,  p-value: 1.439e-09
    ## 
    ## [1] "1 I(A_U * A_EPSILON)"
    ## [1] "2 I(A_U * A_T)"
    ## 
    ## Call:
    ## lm(formula = WAVE1_AYOTTE_24SC_Ay.theta_5_6 ~ A_EPSILON + A_T + 
    ##     A_U + I(A_T * A_EPSILON), data = tData)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.125503 -0.019057  0.002386  0.026174  0.097011 
    ## 
    ## Coefficients:
    ##                     Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)        -0.466250   0.009099 -51.240  < 2e-16 ***
    ## A_EPSILON           0.175425   0.015790  11.110 3.68e-11 ***
    ## A_T                -0.076189   0.015710  -4.850 5.50e-05 ***
    ## A_U                -0.052910   0.015894  -3.329   0.0027 ** 
    ## I(A_T * A_EPSILON)  0.064706   0.030065   2.152   0.0412 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04981 on 25 degrees of freedom
    ## Multiple R-squared:  0.8705, Adjusted R-squared:  0.8498 
    ## F-statistic: 42.01 on 4 and 25 DF,  p-value: 9.526e-11

Now to save we call

``` r
save_ExUQmogp(TestEm, filename = "SavedEmulator")
```

It is important that the filename here has no extensions. Within your
directory, this function has now saved 2 files: `SavedEmulator.RData` is
the Rlist and `SavedEmulator_mogp` is the Python object.

To load an emulator, both the R and Python objects have to exist in the
directory (so you would be loading an emulator that was saved using
`save_ExUQmogp`). Removing `TestEm` to check the call we have

``` r
rm(TestEm)
TestEm <- load_ExUQmogp("SavedEmulator")
```

Checking successful loading:

``` r
newDesign <- 2*randomLHS(100,3)-1
preds <- TestEm$mogp$predict(newDesign, deriv=FALSE)
preds$mean
```

    ##             [,1]        [,2]       [,3]        [,4]        [,5]
    ## [1,] 306.0219210 305.6973448 306.144318 307.0975281 306.9691922
    ## [2,]  -0.4211346  -0.2650799  -0.367697  -0.5929042  -0.5568523
    ##             [,6]        [,7]        [,8]       [,9]       [,10]
    ## [1,] 306.4532875 306.0202856 306.0838060 306.622345 306.6490272
    ## [2,]  -0.4266806  -0.3056318  -0.3671993  -0.570047  -0.4978625
    ##            [,11]      [,12]       [,13]       [,14]       [,15]
    ## [1,] 305.6230907 306.705455 305.6879727 307.1867569 306.0159737
    ## [2,]  -0.3376158  -0.586262  -0.3585313  -0.7057138  -0.3890823
    ##            [,16]       [,17]      [,18]      [,19]       [,20]       [,21]
    ## [1,] 306.2738743 307.2678981 306.865798 306.371496 306.2564627 305.5319005
    ## [2,]  -0.3160374  -0.6698086  -0.550335  -0.426089  -0.3328492  -0.2768971
    ##           [,22]       [,23]      [,24]       [,25]       [,26]       [,27]
    ## [1,] 306.181522 305.9641676 307.500365 306.0692302 306.2758465 306.3419297
    ## [2,]  -0.428255  -0.3156126  -0.740213  -0.2843562  -0.4300289  -0.3501341
    ##            [,28]       [,29]       [,30]      [,31]       [,32]
    ## [1,] 307.5079830 306.3665314 306.7579515 306.063023 306.3223994
    ## [2,]  -0.6627713  -0.4396135  -0.5207699  -0.363739  -0.4159208
    ##            [,33]       [,34]       [,35]       [,36]       [,37]
    ## [1,] 307.2215527 306.4134038 306.4237424 306.1211470 306.3717343
    ## [2,]  -0.5843013  -0.3825297  -0.5000733  -0.3124263  -0.4045056
    ##            [,38]       [,39]       [,40]       [,41]       [,42]
    ## [1,] 306.4958289 306.3426264 305.8842076 306.9088159 306.2090776
    ## [2,]  -0.4911815  -0.4382054  -0.3229232  -0.6209635  -0.4545825
    ##           [,43]       [,44]       [,45]       [,46]       [,47]
    ## [1,] 306.772238 305.5587501 306.2015189 306.2637945 306.6474336
    ## [2,]  -0.611154  -0.3006546  -0.4098229  -0.3960148  -0.5600961
    ##            [,48]       [,49]       [,50]       [,51]       [,52]
    ## [1,] 306.3807224 305.7335969 306.9181139 306.9030389 306.5360954
    ## [2,]  -0.4576071  -0.3559677  -0.6199427  -0.5368094  -0.5198106
    ##            [,53]       [,54]       [,55]       [,56]       [,57]
    ## [1,] 306.1075841 306.1644704 306.3678451 306.2562300 306.4275328
    ## [2,]  -0.3967853  -0.4348353  -0.3354886  -0.4467822  -0.4722101
    ##            [,58]       [,59]       [,60]       [,61]       [,62]
    ## [1,] 305.7236535 306.6222794 306.2840862 306.1085230 306.4205846
    ## [2,]  -0.3650616  -0.5390255  -0.3439853  -0.3949218  -0.4543767
    ##            [,63]       [,64]       [,65]       [,66]       [,67]
    ## [1,] 306.1708736 306.3013217 306.6269890 306.1652252 307.1543438
    ## [2,]  -0.3480482  -0.4508516  -0.5426265  -0.3457983  -0.6266626
    ##            [,68]       [,69]       [,70]       [,71]       [,72]
    ## [1,] 306.2890326 307.1977145 306.5703311 306.7116652 306.7436590
    ## [2,]  -0.2781456  -0.6620267  -0.4098238  -0.5187069  -0.5675807
    ##            [,73]       [,74]       [,75]      [,76]       [,77]
    ## [1,] 307.1136562 306.6425476 307.5166046 306.810018 306.3892351
    ## [2,]  -0.5894745  -0.4921519  -0.6906846  -0.608987  -0.4602103
    ##            [,78]       [,79]       [,80]       [,81]       [,82]
    ## [1,] 306.1848345 306.3061682 307.3602148 306.8057425 306.0801601
    ## [2,]  -0.4103142  -0.4113795  -0.7373195  -0.6016241  -0.3136388
    ##            [,83]       [,84]       [,85]       [,86]       [,87]
    ## [1,] 307.5806088 306.4744420 306.1908663 306.9205447 306.5459928
    ## [2,]  -0.7210212  -0.5091084  -0.4284909  -0.6349606  -0.4854509
    ##            [,88]       [,89]       [,90]       [,91]       [,92]
    ## [1,] 305.4927634 306.4127937 306.1494338 306.6147497 306.1759965
    ## [2,]  -0.2819085  -0.3731011  -0.3669052  -0.5148394  -0.4108669
    ##           [,93]       [,94]       [,95]       [,96]       [,97]
    ## [1,] 306.460895 306.1331615 307.2662352 307.1437253 306.2739008
    ## [2,]  -0.469893  -0.4234563  -0.6167827  -0.7004072  -0.4690272
    ##            [,98]       [,99]      [,100]
    ## [1,] 306.3113668 306.5578412 306.5973617
    ## [2,]  -0.4799531  -0.4899127  -0.5066341
