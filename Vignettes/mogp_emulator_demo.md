`mogp_emulator`
---------------

In this Markdown, we demonstrate the functionalities of
[`mogp_emulator`](https://github.com/alan-turing-institute/mogp_emulator),
a Python package for fitting Gaussian Process Emulators to computer
simulator outputs. This demo is fully based on the `gp_demo.R` available
from [the Alan Turing Institute GitHub
repository](https://github.com/alan-turing-institute/mogp_emulator/tree/devel/mogp_emulator/demos)
and written by Eric Daub. For more detailed information about all the
functions of `mogp_emulator`, please visit [Multi-Output GP Emulator
0.2.0
documentation](https://mogp-emulator.readthedocs.io/en/latest/MultiOutputGP.html).

### Preliminaries

We start by specifying the directory where `mogp_emulator` is installed
so that the python is correctly imported. Please note that your
directory will be different from mine.

``` r
mogp_dir <- "~/Dropbox/BayesExeter/mogp_emulator"
```

``` r
setwd('..')
source('BuildEmulator/BuildEmulator.R')
```

`GaussianProcess` class
-----------------------

We start by considering how to construct a GP emulator for a single
output.

### Data

We proceed to create some data of size 10 with two inputs `x_1` and
`x_2` and output `y`. We store the inputs and output values in a data
frame `x`.

``` r
n_train <- 10
x_scale <- 2.
x1 <- runif(n_train)*x_scale
x2 <- runif(n_train)*x_scale
y1 <- exp(-x1**2 - x2**2)
x <- data.frame(x1, x2, y1)
head(x)
```

    ##          x1         x2         y1
    ## 1 1.9365627 0.42818763 0.01957269
    ## 2 0.4790881 0.06261196 0.79180060
    ## 3 0.4509908 1.58926842 0.06527366
    ## 4 1.2495470 0.09197784 0.20808103
    ## 5 1.6159036 0.08964692 0.07286252
    ## 6 1.8199133 0.71165914 0.02195977

To construct an emulator object, `mogp_emulator` requires data (inputs
and output(s) values) as a matrix, but often you may want to construct a
regression to obtain the form of the regression functions, i.e.
*h*(*x*<sub>*i*</sub>), *i* = 1, …, *q*, using a data frame in `R`. To
do this, we can split this data frame into `inputs`, `targets`, and a
`dictionary` mapping column names to integer indices (note in `Python`
starting from zero) using the function `extract_targets`.

``` r
target_list <- extract_targets(x, target_cols = c("y1"))

inputs <- target_list[[1]]
targets <- target_list[[2]]
inputdict <- target_list[[3]]
```

``` r
print(inputs)
```

    ##             [,1]       [,2]
    ##  [1,] 1.93656274 0.42818763
    ##  [2,] 0.47908812 0.06261196
    ##  [3,] 0.45099076 1.58926842
    ##  [4,] 1.24954704 0.09197784
    ##  [5,] 1.61590356 0.08964692
    ##  [6,] 1.81991327 0.71165914
    ##  [7,] 1.55828524 0.44048660
    ##  [8,] 0.09305291 0.91000597
    ##  [9,] 1.20664915 1.74182010
    ## [10,] 0.07556575 0.58592760

``` r
print(targets)
```

    ##  [1] 0.01957269 0.79180060 0.06527366 0.20808103 0.07286252 0.02195977
    ##  [7] 0.07263696 0.43310733 0.01122132 0.70537808

``` r
print(inputdict)
```

    ## {'x1': 0, 'x2': 1}

### Model specification

We can specify the mean function formula as a string. This formula
specification is adopted in `R`. Note that we include both inputs and
their first-order interaction inside the mean function.

``` r
mean_func <- "y1 ~ x1 + x2 + I(x1*x2)"
```

Priors are specified by giving a list of prior objects (or NULL if you
wish to use weak prior information). Each distribution has some
parameters to set `NormalPrior` is `(mean, std)`, `Gamma` is
`(shape, scale)`, and `InvGammaPrior` is `(shape, scale)`.

If you don’t know how many parameters you need to specify, it depends on
the mean function and the number of input dimensions. Mean functions
have a fixed number of parameters (though in some cases this can depend
on the dimension of the inputs as well), and then covariance functions
have one correlation length per input dimension plus a covariance scale
and a nugget parameter.

If in doubt, you can create the GP instance with no priors, use
`gp$n_params` to get the number, and then set the priors manually using
`gp$priors <- priors`

In this case, we have 4 mean function parameters (normal distribution on
a linear scale), 2 correlations lengths (normal distribution on a log
scale, so lognormal), a sigma^2 covariance parameter (inverse gamma) and
a nugget (Gamma). If you choose an adaptive or fixed nugget, the nugget
prior is ignored.

``` r
priors <- list(mogp_priors$NormalPrior(0., 1.),
               mogp_priors$NormalPrior(0., 1.),
               mogp_priors$NormalPrior(0., 1.),
               mogp_priors$NormalPrior(0., 1.),
               mogp_priors$NormalPrior(0., 1.),
               mogp_priors$NormalPrior(0., 1.),
               mogp_priors$InvGammaPrior(2., 1.), 
               mogp_priors$GammaPrior(1., 0.2))
```

### Create a GP instance

We proceed to create the GP instance. We start by fitting a single
output emulator and therefore we use `GaussianProcess` class.

Note that we provide a fixed value for the nugget. However, we have an
option to set `nugget="adaptive"`, to make the noise only as large as
necessary to successfully invert the covariance matrix, while
`nugget="fit"` to treat nugget as GP hyperparameter that we are
interested to estimate.

``` r
gp <- mogp_emulator$GaussianProcess(inputs, targets,
                                    mean=mean_func,
                                    priors=priors,
                                    nugget=0.0001,
                                    inputdict=inputdict)
```

`gp` is fit using `fit_GP_MAP` function. It accepts an object of
`GaussianProcess` class or `MultiOutputGP` class and returns the same
type of object with the hyperparameters fit via MAP (maximum a
posteriori) estimation.

``` r
gp <- mogp_emulator$fit_GP_MAP(gp)
```

``` r
names(gp)
```

    ##  [1] "current_logpost" "D"               "fit"            
    ##  [4] "get_K_matrix"    "inputs"          "invQt"          
    ##  [7] "kernel"          "L"               "logpost_deriv"  
    ## [10] "logpost_hessian" "logposterior"    "mean"           
    ## [13] "n"               "n_params"        "nugget"         
    ## [16] "nugget_type"     "predict"         "priors"         
    ## [19] "targets"         "theta"

We derive the information about the current value of log-posterior
(`current_logpost`) computed at the current values of GP model
hyperparameters (`theta`).

``` r
print(gp$current_logpost)
```

    ## [1] -3.36574

Note that the first four values of a vector `params` correspond to the
regression coefficients computed on a linear scale (no transformation).
Meanwhile, the remaining four values of a vector `params` correspond to
two correlation length parameters values, sigma^ and a nugget computed
on a log scale. We have to negate and then exponentiate these values to
obtain these values on a linear scale.

``` r
print(gp$theta)
```

    ## [1]  0.9735752 -0.4809720 -0.5553553  0.2968787 -0.3406861 -0.5183254
    ## [7] -1.6095716 -2.1550027

``` r
params <- gp$theta
# Regression coefficients (parameters on a linear scale)
print(params[1:4])
```

    ## [1]  0.9735752 -0.4809720 -0.5553553  0.2968787

``` r
# Correlation length, sigma^2 and nugget (parameters on a log scale)
print(exp(-params[5:length(params)]))
```

    ## [1] 1.405912 1.679213 5.000669 8.627913

Another important remark is that the nugget value provided in the
`params` cannot be trusted if you provided a fixed value in your GP
instance specification. If you want to check your nugget value, do the
following

``` r
print(gp$nugget)
```

    ## [1] 1e-04

### Predictions

Now create some test data to make predictions and compare with unknown
values

``` r
n_test <- 10000

x1_test <- runif(n_test)*x_scale
x2_test <- runif(n_test)*x_scale

x_test <- cbind(x1_test, x2_test)
y_actual <- exp(-x1_test**2 - x2_test**2)

y_predict <- gp$predict(x_test)
```

`y_predict` is an object holding the mean, variance and derivatives (if
computed). Users can assess the values via `y_predict$mean`,
`y_predict$unc`, and `y_predict$deriv`. We compute a root mean squared
value (RMSE). In our example, this value is small and close to zero,
which indicates the accuracy in prediction of obtained emulator.

``` r
print(sqrt(sum((y_actual - y_predict$mean)**2)/n_test))
```

    ## [1] 0.0308244

`MultiOutputGP` class
---------------------

We proceed to construct two emulators for two outputs using
`MultiOutputGP` class.

### Data

We start by generating data and making sure that it is in a right
format.

``` r
y2 <- sin(x1)+exp(x2)*2
x <- data.frame(x1, x2, y1, y2)
target_list <- extract_targets(x, target_cols = c("y1", "y2"))
inputs <- target_list[[1]]
targets <- target_list[[2]]
inputdict <- target_list[[3]]
```

As before, we split a data frame into `inputs`, which is a matrix with
its columns corresponding to inputs’ values.

``` r
head(inputs)
```

    ##           [,1]       [,2]
    ## [1,] 1.9365627 0.42818763
    ## [2,] 0.4790881 0.06261196
    ## [3,] 0.4509908 1.58926842
    ## [4,] 1.2495470 0.09197784
    ## [5,] 1.6159036 0.08964692
    ## [6,] 1.8199133 0.71165914

`targets` is matrix where each row corresponds to one of the outputs.

``` r
print(targets)
```

    ##            [,1]      [,2]        [,3]     [,4]       [,5]       [,6]
    ## [1,] 0.01957269 0.7918006  0.06527366 0.208081 0.07286252 0.02195977
    ## [2,] 4.00279785 2.5901974 10.23618291 3.141523 3.18655889 5.04386797
    ##            [,7]      [,8]        [,9]     [,10]
    ## [1,] 0.07263696 0.4331073  0.01122132 0.7053781
    ## [2,] 4.10684765 5.0615934 12.34987308 3.6688074

Finally, `dictionary` contains the mapping of column names to integer
indices. Note that indices start from 0, since `mogp_emulator` is a
Python package.

``` r
print(inputdict)
```

    ## {'x1': 0, 'x2': 1}

### Model specification

`MultiOutputGP` has an option to specify mean function, priors, kernel
type and nugget treatment for each emulator. For mean function
specification, we creat a list with two string formulaes.

``` r
mean_func <- list("y ~ x1 + x2 + I(x1*x2)", 
                  "y ~ x1 + x2 + I(x1^3)")
```

For prior specification, we can provide a list of lists of priors or
none.

``` r
priors_mult <- list(priors, priors)
```

We can specify different form of kernel for each emulator (squared
exponential or Matern52), i.e.

``` r
kernel_mult <- list(mogp_kernels$SquaredExponential(), 
                    mogp_kernels$Matern52())
```

Finally, we can treat nugget in a different way for each emulator.

``` r
nugget_mult <- list(0.0001, "adaptive")
```

### Create a GP instance

``` r
gp_multi <- mogp_emulator$MultiOutputGP(inputs, targets, 
                                        inputdict = inputdict, 
                                        mean = mean_func, 
                                        priors = priors_mult, 
                                        nugget=nugget_mult)
```

`gp_multi` is fit using `fit_GP_MAP` function. It accepts an object of
class `MultiOutputGP` class and returns the same type of object with
hyperparameters fit via MAP (maximum a posteriori) estimation

``` r
gp_multi <- mogp_emulator$fit_GP_MAP(gp_multi)
names(gp_multi)
```

    ## [1] "D"           "emulators"   "n"           "n_emulators" "predict"

All the information about each individual emulator is stored in
`gp_multi$emulators`, which is a list. For example, here is all the
information about an emulator for the first output, `y_1`.

``` r
names(gp_multi$emulators[[1]])
```

    ##  [1] "current_logpost" "D"               "fit"            
    ##  [4] "get_K_matrix"    "inputs"          "invQt"          
    ##  [7] "kernel"          "L"               "logpost_deriv"  
    ## [10] "logpost_hessian" "logposterior"    "mean"           
    ## [13] "n"               "n_params"        "nugget"         
    ## [16] "nugget_type"     "predict"         "priors"         
    ## [19] "targets"         "theta"

### Predictions

We can generate predictions using these two emulators. We create some
test data to make predictions and compate with unknown values.

``` r
n_test <- 10000

x1_test <- runif(n_test)*x_scale
x2_test <- runif(n_test)*x_scale

x_test <- cbind(x1_test, x2_test)
y1_actual <- exp(-x1_test**2 - x2_test**2)
y2_actual <- sin(x1_test) + exp(x2_test)*2

y_predict <- gp_multi$predict(x_test)
```

We proceed to compute RMSE for the first emulator and second emulator

``` r
print(sqrt(sum((y1_actual - y_predict$mean[1, ])**2)/n_test))
```

    ## [1] 0.03166501

``` r
print(sqrt(sum((y2_actual - y_predict$mean[2, ])**2)/n_test))
```

    ## [1] 0.5033606
