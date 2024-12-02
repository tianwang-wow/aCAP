\name{aCAP-package}
\alias{aCAP-package}
\docType{package}
\title{
  aCAP: Adjustable Combination of Amplitude and Phase Model for Functional Data
}
\description{
  Fit an adjustable combination of amplitude and phase model for functional data with given alpha value. The parameter alpha is a modelling choice for balancing the amplitude and phase components.
}
\author{
Tian Wang, tianwang@wustl.edu.

Jimin Ding, jmding@wustl.edu.

Maintainer: Tian Wang <tianwang@wustl.edu>
}
\references{
  https://doi.org/10.1093/jrsssb/qkae112
}
\keyword{ functional data, warping }
\examples{
  \dontrun{
  ### Do pairwise warping
  t = seq(0, 1, length.out = 100+1) # time
  fdata1 = cbind(t, 4*exp(-(t-0.3)^2 / 0.01) + 6*exp(-(t-0.7)^2 / 0.03)) # curve 1
  fdata2 = cbind(t, 7*exp(-(t-0.4)^2 / 0.03) + 3*exp(-(t-0.8)^2 / 0.01)) # curve 2
  fit = aCAPf1f2(func1 = fdata1, func2 = fdata2, alpha = 0.5)
  
  ### Fit aCAP for a set of curves
  set.seed(12345)
  n = 10 # number of curves
  t = seq(0, 1, length.out = 100+1) # time
  n_t = length(t)
  fdata = matrix(t, nrow = n_t, ncol = n+1)
  for (j in 1:n) {
    fdata[, j+1] =
      rnorm(1, 1.5, 0.2) * exp(-(t-runif(1, 0.2, 0.4))^2 / 0.01) +
      rnorm(1, 2, 0.2) * exp(-(t-runif(1, 0.6, 0.8))^2 / 0.01)
  }
  fit = aCAP(fdata = fdata, alpha = seq(0, 1, 0.25), plot = FALSE)
  
  }
}
