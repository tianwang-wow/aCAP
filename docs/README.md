aCAP
================

Fit adjustable combination of amplitude and phase (aCAP) model for
functional data.

## Installation

``` r
devtools::install_github("tianwang-wow/aCAP")
library(aCAP)
```

## Main functions

**`aCAPf1f2`**: Find an optimal warping function from curve 1 to curve
2.

- Argument
  - `func1`: A matrix for curve 1 with time (equidistant) in the first
    column and function values in the second column.
  - `func2`: A matrix for curve 2 with time (equidistant) in the first
    column and function values in the second column.
  - `alpha`: A vector of alpha values. Default = 0.5.
  - `plots`: Whether to plot the warping results. Default = TRUE.
  - `verbose`: Whether to output the detailed procedure and distances.
    Default = FALSE.
  - `n_digit`: Number of decimal points in the output. Default = 5.
  - `max_drv`: An integer specifying the maximum derivative allowed for
    the warping function. Default = 10. A smaller value decreases
    computation time.
- Value
  - A list of warping results for each alpha.
  - `warping_func`: A matrix for the warping function with time in the
    first column and function values in the second column.
  - `warped_func`: A matrix for the warped function with time in the
    first column and function values in the second column.
  - `d_A`: Amplitude shape distance.
  - `d_P`: Phase distance.
  - `d_C`: Optimal combined distance.
  - `alpha`: The alpha value(s) used.

**`aCAP`**: Fit aCAP to a set of functions.

- Argument
  - `fdata`: A matrix for functions with time (equidistant) in the first
    column and function values in the remaining columns.
  - `alpha`: A vector of alpha values. Default = 0.5.
  - `n_ite`: Maximum number of iterations.
  - `plots`: Whether to plot the warping results. Default = TRUE.
  - `verbose`: Whether to output the detailed procedure. Default =
    FALSE.
  - `n_digit`: Number of decimal points in the output. Default = 5.
  - `stopping`: Criterion for the stopping rule. Default = 0.02.
  - `col`: Color of plotted functions. Default = ‘gray60’.
  - `max_drv`: An integer specifying the maximum derivative allowed for
    the warping function. Default = 10. A smaller value decreases
    computation time.
- Value
  - A list of aCAP results for each alpha.
  - `warping_func`: A matrix for the warping functions with time in the
    first column and function values in the remaining columns.
  - `warped_func`: A matrix for the warped functions with time in the
    first column and function values in the remaining columns.
  - `V_A`: Warping-related amplitude variation.
  - `V_P`: Phase variation.
  - `alpha`: The alpha value(s) used.
  - `total_ite`: Number of iterations.

## Example of warping from one curve to another

<!-- Below is a step-by-step example of finding optimal warping function from one curve to another. -->

Simulate 2 curves:

``` r
t = seq(0, 1, length.out = 100+1) # time
fdata1 = cbind(t, 4*exp(-(t-0.3)^2 / 0.01) + 6*exp(-(t-0.7)^2 / 0.03)) # curve 1
fdata2 = cbind(t, 7*exp(-(t-0.4)^2 / 0.03) + 3*exp(-(t-0.8)^2 / 0.01)) # curve 2
```

Find an optimal warping from curve 1 to curve 2 under alpha = 0.5:

``` r
fit = aCAPf1f2(func1 = fdata1, func2 = fdata2, alpha = 0.5, verbose = F)
```

<img src="./docs/README_files/figure-gfm/find optimal warping-1.png" style="display: block; margin: auto;" />

In the right panel, curve 1 is the solid curve, curve 2 is the dashed
curve, and the warped curve from curve 1 to curve 2 is the dotted one.

## Example of fitting aCAP to a set of curves

<!-- Below is a step-by-step example of fitting aCAP for a set of curves. -->

Simulate a set of 20 curves:

``` r
set.seed(12345)
n = 20 # number of curves
t = seq(0, 1, length.out = 100+1) # time
fdata = matrix(t, nrow = length(t), ncol = n+1)
for (j in 1:n) {
  fdata[, j+1] =
    rnorm(1, 1.5, 0.2) * exp(-(t-runif(1, 0.2, 0.4))^2 / 0.01) +
    rnorm(1, 2, 0.2) * exp(-(t-runif(1, 0.6, 0.8))^2 / 0.01)
}
```

Fit aCAP to the simulated 20 curves using alpha = 0.5:

``` r
fit = aCAP(fdata = fdata, alpha = 0.5, verbose = F, plot = F)
```

Plot separation results:

``` r
par(mfrow = c(1,3))

plot(t, fdata[,2], main = 'observed functions', ylab = 'x(t)', col = 'gray60', type = 'l', ylim = c(0, 2.7))
for (j in 2:n) { lines(t, fdata[,j+1], type = 'l', col = 'gray60') }

plot(t, fit$"0.5"$warping_func[,2], main = 'warping functions', ylab = 'h(t)', col = 'gray60', type = 'l')
for (j in 2:n) { lines(t, fit$"0.5"$warping_func[,j+1], type = 'l', col = 'gray60') }

plot(t, fit$"0.5"$warped_func[,2], main = 'warped functions', ylab = 'x(h(t))', col = 'gray60', type = 'l', ylim = c(0, 2.7))
for (j in 2:n) { lines(t, fit$"0.5"$warped_func[,j+1], type = 'l', col = 'gray60') }
```

<img src="./docs/README_files/figure-gfm/plot results-1.png" style="display: block; margin: auto;" />
