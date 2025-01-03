% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func.R
\name{aCAP}
\alias{aCAP}
\title{Fitting aCAP model to a set of curves}
\usage{
aCAP(
  fdata,
  alpha = 0.5,
  n_ite = 100,
  plots = TRUE,
  verbose = FALSE,
  n_digit = 5,
  stopping = 0.02,
  col = "gray60",
  max_drv = 10
)
}
\arguments{
\item{fdata}{A matrix for functions with time (equidistant) in the first column and function values in the remaining columns.}

\item{alpha}{A vector of alpha values. Default = 0.5.}

\item{n_ite}{Maximum number of iterations.}

\item{plots}{Whether to plot the warping results. Default = TRUE.}

\item{verbose}{Whether to output the detailed procedure. Default = FALSE.}

\item{n_digit}{Number of decimal points in the output. Default = 5.}

\item{stopping}{Criterion for the stopping rule. Default = 0.02.}

\item{col}{Color of plotted functions. Default = 'gray60'.}

\item{max_drv}{An integer specifying the maximum derivative allowed for the warping function. Default = 10. A smaller value decreases computation time.}
}
\value{
A list of aCAP results for each alpha.
\itemize{
  \item \code{warping_func}: A matrix for the warping functions with time in the first column and function values in the remaining columns.
  \item \code{warped_func}: A matrix for the warped functions with time in the first column and function values in the remaining columns.
  \item \code{V_A}: Warping-related amplitude variation.
  \item \code{V_P}: Phase variation.
  \item \code{alpha}: The alpha value(s) used in model fitting.
  \item \code{total_ite}: Number of iterations.
}
}
\description{
Fitting aCAP model to a set of curves
}
\examples{
set.seed(12345)
n = 20 # number of curves
t = seq(0, 1, length.out = 100+1) # time
n_t = length(t)
fdata = matrix(t, nrow = n_t, ncol = n+1)
for (j in 1:n) {
  fdata[, j+1] =
    rnorm(1, 1.5, 0.2) * exp(-(t-runif(1, 0.2, 0.4))^2 / 0.01) +
   rnorm(1, 2, 0.2) * exp(-(t-runif(1, 0.6, 0.8))^2 / 0.01)
}
fit = aCAP(fdata = fdata, alpha = seq(0, 1, 0.25), verbose = FALSE, plot = FALSE, max_drv = 5)
}
