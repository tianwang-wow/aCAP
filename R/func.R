#' Pairwise warping using the optimal combined distance
#'
#' @param func1 A matrix for curve 1 with time (equidistant) in the first column and function values in the second column.
#' @param func2 A matrix for curve 2 with time (equidistant) in the first column and function values in the second column.
#' @param alpha A vector of alpha values. Default = 0.5.
#' @param plots Whether to plot the warping results. Default = TRUE.
#' @param verbose Whether to output the procedure. Default = FALSE.
#' @param n_digit Number of decimal points in the output. Default = 5.
#' @param max_drv An integer of the maximum derivative allowed for warping function. Default = 10. A smaller value decreases computation time.
#'
#' @return A list of warping results for each alpha.
#' \itemize{
#'   \item \code{warping_func}: A matrix for warping function with time in the first column and function values in the second column.
#'   \item \code{warped_func}: A matrix for warped function with time in the first column and function values in the second column.
#'   \item \code{d_A}: Amplitude shape distance.
#'   \item \code{d_P}: Phase distance.
#'   \item \code{d_C}: Optimal combined distance.
#'   \item \code{alpha}: The alpha values used in warping.
#' }
#' @export
#'
#' @examples
#' t = seq(0, 1, length.out = 100+1) # time
#' fdata1 = cbind(t, 4*exp(-(t-0.3)^2 / 0.01) + 6*exp(-(t-0.7)^2 / 0.03)) # curve 1
#' fdata2 = cbind(t, 7*exp(-(t-0.4)^2 / 0.03) + 3*exp(-(t-0.8)^2 / 0.01)) # curve 2
#' fit = aCAPf1f2(func1 = fdata1, func2 = fdata2, alpha = 0.5, verbose = FALSE)
aCAPf1f2 = function(func1,
                    func2,
                    alpha = 0.5,
                    plots = TRUE,
                    verbose = FALSE,
                    n_digit = 5,
                    max_drv = 10)
{
  if (max_drv > 10) { warning('The max of derivative may be too large due to computation intensity. Consider setting a smaller value.') }
  t = func1[,1]
  dfunc1 = drv(func1)
  dfunc2 = drv(func2)
  Jfunc1_ext_ls = est_Jfunc(dfdata = dfunc1,
                            subs = 1:max_drv)
  Jfunc2_ext_ls = est_Jfunc(dfdata = dfunc2,
                            subs = 1:max_drv)
  
  ### Do iteration
  output = vector('list', length(alpha))
  names(output) = alpha
  for (v in 1:length(alpha)) {
    current_alpha = alpha[v]
    d_A = d_P = d_C = NULL
    warping_func = aCAPf1f2_cpp(Jfunc1_ext_ls = Jfunc1_ext_ls,
                                Jfunc2_ext_ls = Jfunc2_ext_ls,
                                alpha = current_alpha,
                                max_drv = max_drv)
    warped_func = comp(func1, warping_func)
    
    ### Find d_A and d_P
    if (verbose) {
      d_A = Lp(cbind(t, Jfunc2_ext_ls$`1`[,2] - Jfunc(fdata = warped_func)[,2]), p = 2) / sqrt(2)
      d_P = Lp(cbind(t, sqrt(drv(warping_func)[,2])-1), p = 2) / sqrt(2)
      d_C = sqrt(current_alpha * d_A^2 + (1-current_alpha) * d_P^2)
      cat(paste('\n*** alpha =',current_alpha),'\n');
      cat(paste('d_A:', dec(d_A, n_digit),'  ')); cat(paste('d_P:', dec(d_P, n_digit)),'  ')
    }
    ### Draw plots
    if (plots) {
      graphics::par(mfrow = c(1,2))
      ### Plot warping function
      plot(t, warping_func[,2], type = 'l', ylab = 'h(t)', lwd = 1.5,
           main =  paste0('warping function (alpha=', current_alpha, ')'))
      graphics::lines(t, t, lty = 4, lwd = 0.5)
      ### Plot warped function
      plot(t, func1[,2], type = 'l', ylim = range(func1[,2], func2[,2]), ylab = 'x(h(t)) & y(t)', lwd = 1.5,
           main = paste0('warping results (alpha=', current_alpha, ')'))
      graphics::lines(t, func2[,2], lty = 2, lwd = 1.5)
      graphics::lines(t, warped_func[,2], lty = 3, lwd = 1.5)
    }
    rownames(warping_func) = rownames(warped_func) = rownames(func1)
    colnames(warping_func) = colnames(warped_func) = colnames(func1)
    
    output[[v]] = list(warping_func = warping_func,
                       warped_func = warped_func,
                       d_A = d_A,
                       d_P = d_P,
                       d_C = d_C,
                       alpha = current_alpha)
  }
  return(output)
}

#' Fitting aCAP model on a set of curves
#'
#' @param fdata A matrix for functions with time (equidistant) in the first column and function values in the remaining columns.
#' @param alpha A vector of alpha values. Default = 0.5.
#' @param n_ite Maximum number of iterations.
#' @param plots Whether to plot the warping results. Default = TRUE.
#' @param verbose Whether to output the procedure. Default = FALSE.
#' @param n_digit Number of decimal points in the output. Default = 5.
#' @param stopping Stopping rule criterion. Default = 0.02.
#' @param col Color of plotted functions. Default = 'gray60'.
#' @param max_drv An integer of the maximum derivative allowed for the warping function. Default = 10.
#'
#' @return A list of aCAP results for each alpha.
#' \itemize{
#'   \item \code{warping_func}: A matrix for warping functions with time in the first column and function values in the remaining columns.
#'   \item \code{warped_func}: A matrix for warped functions with time in the first column and function values in the remaining columns.
#'   \item \code{V_A}: Warping-related amplitude variation.
#'   \item \code{V_P}: Phase variation.
#'   \item \code{alpha}: The alpha values used in model fitting.
#'   \item \code{total_ite}: Number of iterations.
#' }
#' @export
#'
#' @examples
#' set.seed(12345)
#' n = 20 # number of curves
#' t = seq(0, 1, length.out = 100+1) # time
#' n_t = length(t)
#' fdata = matrix(t, nrow = n_t, ncol = n+1)
#' for (j in 1:n) {
#'   fdata[, j+1] =
#'     rnorm(1, 1.5, 0.2) * exp(-(t-runif(1, 0.2, 0.4))^2 / 0.01) +
#'    rnorm(1, 2, 0.2) * exp(-(t-runif(1, 0.6, 0.8))^2 / 0.01)
#' }
#' fit = aCAP(fdata = fdata, alpha = seq(0, 1, 0.25), verbose = FALSE, plot = FALSE, max_drv = 5)
aCAP = function(fdata,
                   alpha = 0.5,
                   n_ite = 100,
                   plots = TRUE,
                   verbose = FALSE,
                   n_digit = 5,
                   stopping = 0.02,
                   col = 'gray60',
                   max_drv = 10)
{
  if (max_drv > 10) { warning('The max of derivative may be too large due to computation intensity. Consider setting a smaller value.') }
  if (verbose) cat("In prepraration... \n")
  n_fdata = ncol(fdata)-1
  t = fdata[,1]
  dfdata = drv(fdata)
  Jfdata_ext_ls = est_Jfunc(dfdata = dfdata,
                            subs = 1:max_drv)
  
  ### Begin fitting aCAP
  output = vector('list', length(alpha))
  names(output) = alpha
  for (v in 1:length(alpha)) {
    current_alpha = alpha[v]
    if (verbose) { cat(paste('***** alpha =', current_alpha), '***** \n') }
    
    ### Step 1: Initialization
    k = 0 # iteration index
    stop_ite = FALSE
    V_A = V_P = NULL
    warping_func = matrix(rep(t, n_fdata+1), ncol = n_fdata+1)
    while (T) {
      if (verbose) { cat(paste('-- Iteration index:', k, '\n')) }
      
      ### Step 2: Find target mean function
      if (k == 0) {
        warped_func = fdata
      } else {
        warped_func = comp(fdata, warping_func)
      }
      ### Plot warping and warped functions
      if (plots) {
        graphics::par(mfrow = c(1,2))
        plotn(warping_func, xlab = 't', ylab = 'h(t)', col = col, n = n_fdata,
              main = paste0('warping functions (alpha=', current_alpha, ', l=', k, ')'))
        plotn(warped_func, xlab = 't', ylab = 'x(h(t))', col = col, n = n_fdata,
              main = paste0('warped functions (alpha=', current_alpha, ', l=', k, ')'))
      }
      ### Find Frechet mean of warped functions in d_V
      Jwarped_func_ext_ls = est_Jfunc(dfdata = drv(warped_func),
                                      subs = 1:max_drv)
      if (k > 0) { old_Jfmean_ext = Jfmean_ext }
      Jfmean_ext_ls = lapply(Jwarped_func_ext_ls, function(x) fmean_dv(Jfdata = x))
      Jfmean_ext = Jfmean_ext_ls[[1]]
      ### Check stopping rule
      diff_dist = ifelse(k == 0,
                         Inf,
                         Lp(cbind(Jfmean_ext[,1], Jfmean_ext[,2] - old_Jfmean_ext[,2]), 2))
      if (verbose) { cat(paste('Stopping rule:', dec(diff_dist, n_digit)), '\n') }
      if (diff_dist < stopping | k >= n_ite | current_alpha == 0) { stop_ite = T }
      
      ### Step 3: Find optimal warpings
      for (j in 1:n_fdata) {
        current_func_ls = lapply(Jfdata_ext_ls, function(x) x[, c(1, j+1)])
        warping_func[, j+1] = aCAPf1f2_cpp(Jfunc1_ext_ls = current_func_ls,
                                           Jfunc2_ext_ls = Jfmean_ext_ls,
                                           alpha = current_alpha,
                                           max_drv = max_drv)[,2]
      }
      ### Find V_A and V_P
      if (verbose & !stop_ite) {
        warped_func = comp(fdata, warping_func)
        V_A = VA(fdata = warped_func)
        V_P = VP(hdata = warping_func)
        cat(paste('V_A:', dec(V_A, n_digit), '   '))
        cat(paste('V_P:', dec(V_P, n_digit), '   '))
        cat(paste('alpha*V_A + (1-alpha)*V_P =',
                  dec(current_alpha * V_A + (1-current_alpha) * V_P, n_digit)), '\n')
      }
      if (stop_ite) {
        warped_func = comp(fdata, warping_func)
        V_A = VA(fdata = warped_func)
        V_P = VP(hdata = warping_func)
        if (verbose) {
          cat('Final variations: \n')
          cat(paste('V_A:', dec(V_A, n_digit), '   '))
          cat(paste('V_P:', dec(V_P, n_digit), '   '))
          cat(paste('alpha*V_A + (1-alpha)*V_P =',
                    dec(current_alpha * V_A + (1-current_alpha) * V_P, n_digit)), '\n')
          cat('Iterations stopped. \n \n')
        }
        break
      }
      
      ### Step 4: Center optimal warpings
      warping_func = centered_h(warping_func)
      
      k = k+1
    }
    rownames(warping_func) = rownames(warped_func) = rownames(fdata)
    colnames(warping_func) = colnames(warped_func) = colnames(fdata)
    
    output[[v]] = list(warping_func = warping_func,
                       warped_func = warped_func,
                       V_A = V_A,
                       V_P = V_P,
                       alpha = current_alpha,
                       total_ite = k)
  }
  return(output)
}

dec = function(x, k = 3) { return(trimws(format(round(x, k), nsmall = k))) }

plotn = function(fdata,
                 n = 20,
                 xlab = 't',
                 ylab = 'x(t)',
                 ylim = NULL,
                 main = 'functions',
                 col = 'gray60',
                 plotmean = TRUE,
                 lwd_mean = 3.5,
                 ...)
{
  n_func = ncol(fdata)-1
  if (is.null(ylim)) { ylim = range(fdata[,-1]) }
  if (length(col)==1) { col = rep(col, n_func) }
  
  plot(fdata[,1], fdata[,2],
       type = 'l', xlab = xlab, ylim = ylim, ylab = ylab, main = main, col = col[1], ...)
  if (min(n, n_func) >= 2) {
    for (j in 2:min(n_func, n)) { graphics::lines(fdata[,1], fdata[,j+1], type = 'l', col = col[j], ...) }
    if (plotmean) {
      fmean = fmean_L2(fdata)
      graphics::lines(fmean[,1], fmean[,2], lwd = lwd_mean, col = 'black')
    }
  }
}

drv = function(fdata)
{
  nr = nrow(fdata)
  dfdata = rbind(diff(fdata[1:2,]),
                 diff(fdata,2),
                 diff(fdata[(nr-1):nr,]))
  return(cbind(fdata[,1], dfdata[,-1]/dfdata[,1]))
}

inv = function(hdata)
{
  time = hdata[,1]
  inv_hdata = hdata
  for (j in 2:ncol(hdata)) { inv_hdata[,j] = stats::approxfun(hdata[,j], time)(time) }
  inv_hdata[nrow(inv_hdata),] = 1 # make sure h(1) = 1
  return(inv_hdata)
}

comp = function(fdata, hdata)
{
  time = fdata[,1]
  new_fdata = fdata
  n_col = ncol(fdata)
  if (ncol(hdata) == 2) {
    for (j in 2:ncol(new_fdata)) { new_fdata[,j] = stats::approxfun(time, fdata[,j])(hdata[,2]) }
  } else {
    for (j in 2:ncol(new_fdata)) { new_fdata[,j] = stats::approxfun(time, fdata[,j])(hdata[,j]) }
  }
  if (sum(is.na(new_fdata[1, ])) > 0) {
    new_fdata[1, 2:n_col] = fdata[1, 2:n_col]
  }
  nt = nrow(new_fdata)
  if (sum(is.na(new_fdata[nt, ])) > 0) {
    new_fdata[nt, 2:n_col] = fdata[nt, 2:n_col]
  }
  return(new_fdata)
}

int = function(fdata)
{
  nr = nrow(fdata)
  gap = fdata[2,1]
  return(as.vector(rep(gap, nr-1) %*% ((fdata[1:(nr-1),-1] + fdata[2:nr,-1])/2)))
}

int_eacht = function(fdata)
{
  nr = nrow(fdata)
  integral = fdata
  integral[1,] = 0
  diff_t = diff(fdata[,1])
  cha = as.matrix((fdata[1:(nr-1),-1] + fdata[2:nr,-1])/2)
  
  gap = fdata[2,1]
  temp_int = as.matrix(cumsum(data.frame(cha)))
  integral[-1,-1] = temp_int * gap
  
  return(integral)
}

LVD = function(fdata)
{
  time = fdata[,1]
  len_t = length(time)
  abs_dfdata = cbind(time, abs(drv(fdata)[,-1]))
  lvd_fdata = fdata
  lvd_fdata = int_eacht(abs_dfdata)
  lvd_fdata[,-1] = t(t(lvd_fdata[,-1]) / lvd_fdata[len_t,-1])
  return(lvd_fdata)
}

LVQ = function(fdata) { return(inv(LVD(fdata))) }

Lp = function(fdata, p = 2) { return(int(cbind(fdata[,1], abs(fdata[,-1])^p))^(1/p)) }

fmean_L2 = function(fdata) { return(cbind(fdata[,1], rowMeans(fdata[,-1]))) }

Jfunc = function(fdata,
                 dfdata,
                 L1norm_drv)
{
  if (missing(dfdata)) { dfdata = drv(fdata) }
  abs_dfdata = cbind(dfdata[,1], abs(dfdata[,-1]))
  if (missing(L1norm_drv)) { L1norm_drv = int(abs_dfdata) }
  Jfdata = dfdata
  Jfdata[,-1] = sign(dfdata[,-1]) * sqrt(t(t(abs_dfdata[,-1]) / L1norm_drv))
  return(Jfdata)
}

fmean_dv = function(fdata,
                    Jfdata,
                    L1norm_drv)
{
  if (missing(Jfdata)) {
    if (missing(L1norm_drv)) {
      L1norm_drv = int(cbind(fdata[,1], abs(drv(fdata)[,-1])))
    }
    Jfdata = Jfunc(fdata = fdata, L1norm_drv = L1norm_drv)
  }
  top = cbind(Jfdata[,1], rowSums(Jfdata[,-1]))
  bottom = Lp(top, 2)
  return(cbind(Jfdata[,1], top[,2] / bottom))
}

hmean = function(hdata)
{
  temp = cbind(hdata[,1], rowSums(sqrt(drv(hdata)[,-1]))^2)
  top = int_eacht(temp)
  bottom = top[length(top)]
  return(cbind(hdata[,1], top[,2] / bottom))
}

centered_h = function(hdata) { return(comp(hdata, inv(hmean(hdata)))) }

est = function(fdata, new_t)
{
  t = fdata[,1]
  est_new_t = matrix(nrow = length(new_t), ncol = ncol(fdata))
  est_new_t[,1] = new_t
  for (j in 2:ncol(fdata)) { est_new_t[,j] = stats::approxfun(t, fdata[,j])(new_t) }
  return(est_new_t)
}

VA = function(fdata, Jfdata, Jfmean)
{
  if (missing(Jfdata)) { Jfdata = Jfunc(fdata) }
  if (missing(Jfmean)) { Jfmean = fmean_dv(Jfdata = Jfdata) }
  return(mean(int(cbind(Jfdata[,1], (Jfdata[,-1]-Jfmean[,2])^2))/2))
}

VP = function(hdata) { return(mean(int(cbind(hdata[,1], (sqrt(drv(hdata)[,-1])-1)^2))/2)) }

L1_drv = function(fdata) { return(colSums(abs(diff(fdata[, -1, drop = FALSE])))) }

est_Jfunc = function(dfdata,
                     L1norm_drv = NULL,
                     subs)
{
  t = dfdata[,1]
  gap = t[2]
  n_subs = length(subs)
  Jfunc_ext_ls = vector('list', length = n_subs)
  for (b in 1:n_subs) {
    t_ext = seq(0,1, gap * 1/subs[b])
    dfdata_ext = est(dfdata, t_ext)
    if (is.null(L1norm_drv)) {
      Jfunc_ext_ls[[b]] = Jfunc(dfdata = dfdata_ext)
    } else {
      Jfunc_ext_ls[[b]] = Jfunc(dfdata = dfdata_ext, L1norm_drv = L1norm_drv)
    }
  }
  names(Jfunc_ext_ls) = subs
  return(Jfunc_ext_ls)
}

