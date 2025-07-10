#' @title Dynamic Modeling of Sparse Longitudinal Data and Functional Snippets
#' @description This function performs dynamic modeling of sparse longitudinal 
#' data and functional snippets using stochastic differential equations. 
#' The method employs local linear regression to estimate both the conditional 
#' mean and the conditional variance.
#' @param Ly a list of \eqn{n} vectors containing the observed values
#' for each individual.
#' @param Lt a list of \eqn{n} vectors containing the observation time points
#' for each individual corresponding to \code{Ly}. Each vector should be
#' sorted in ascending order.
#' @param z0 initial condition (starting value and starting time) for
#' the underlying stochastic process.
#' @param tp a vector of length \eqn{K} of discretized time points at which
#' the stochastic process is evaluated
#' @param optns a list of options control parameters specified by
#' \code{list(name = value)}. See `Details'.
#' @details Available control options are
#' \describe{
#' \item{M}{a scalar holding the number of Monte Carlo simulations to run. Default is 100.}
#' \item{regular}{whether to assume regular observation time points
#' (time spacing is the same for all individuals). Default is TRUE.}
#' \item{kernel}{smoothing kernel choice, common for mean and variance.
#' Available options are 'gauss', 'rect', 'epan', 'gausvar' and 'quar'.
#' Default is 'gauss'.}
#' \item{bw1}{bandwidth for conditional mean estimation, if not entered
#' it would be chosen from cross validation.}
#' \item{bw2}{bandwidth for conditional variance estimation, if not entered
#' it would be chosen from cross validation.}
#' \item{bm}{Brownian motion path used.}
#' @return A \code{dm} object --- a list containing the following fields:
#' \item{path}{a \eqn{M} by \eqn{K} matrix holding the \eqn{M} estimated sample
#' paths of the stochastic process at the \eqn{K} discretized time points.}
#' \item{Ly}{the original \code{Ly} used.}
#' \item{Lt}{the original \code{Lt} used.}
#' \item{z0}{the initial condition used.}
#' \item{t}{the discretized time points used.}
#' \item{optns}{the control options used.}
#' @examples
#' @references
#' \itemize{
#' \item \cite{Zhou, Y. and Müller, H.G., 2023. Dynamic Modeling of Sparse
#' Longitudinal Data and Functional Snippets With Stochastic Differential
#' Equations. arXiv preprint arXiv:2306.10221.}
#' }
#' @export

ldm = function(Ly = NULL,
               Lt = NULL,
               z0 = NULL,
               tp = NULL,
               optns = list()) {
  start_time = Sys.time()
  
  if (is.null(Ly) | is.null(Lt) | is.null(z0) | is.null(tp)) stop("require the input of Ly, Lt, z0, and tp")
  if (!is.list(Ly) | !is.list(Lt)) stop("Ly and Lt must be lists")
  if (length(Ly) != length(Lt)) stop("Ly and Lt must have the same length")
  if (any(sapply(Ly, length) - sapply(Lt, length))) stop("each individual's observation/time points must match")
  
  if (is.null(optns$M)) optns$M = 100
  if (is.null(optns$cores)) optns$cores = max(1, parallel::detectCores() - 1)
  if (is.null(optns$regular)) optns$regular = TRUE
  if (is.null(optns$kernel)) optns$kernel = "gauss"
  
  cl = parallel::makeCluster(optns$cores)
  on.exit(parallel::stopCluster(cl))
  
  if (is.vector(z0)) {
    if (z0[2] != tp[1]) stop("starting time must match first grid point")
    z0 = matrix(rep(z0, each = optns$M), ncol = 2)
  } else {
    if (any(unique(z0[, 2]) != tp[1])) stop("starting time must match first grid point")
  }
  
  K = length(tp)
  if (is.null(optns$bm)) {
    optns$bm = matrix(rnorm((K - 1) * optns$M), nrow = optns$M)
  } else if (is.vector(optns$bm)) {
    optns$bm = t(optns$bm)
  }
  
  n = length(Ly)
  y = NULL
  z = NULL
  for (i in 1:n) {
    Ni = length(Ly[[i]])
    y = c(y, Ly[[i]][-1])
    if (optns$regular) {
      z = rbind(z, cbind(Ly[[i]][-Ni], Lt[[i]][-Ni]))
    } else {
      z = rbind(z, cbind(Ly[[i]][-Ni], Lt[[i]][-Ni], Lt[[i]][-1]))
    }
  }
  p = ncol(z)
  colnames(z) = c("y1", "t", "s")[1:p]
  p0 = min(p, 2)
  zt = z[, 1:p0]
  N = length(y)
  
  kern = kerFctn(optns$kernel)
  KF = function(x, h) {
    prod(sapply(seq_along(h), function(i) kern(x[i] / h[i])))
  }
  
  if (is.null(optns$bw1)) {
    hs = matrix(0, p0, 20)
    for (l in 1:p0) {
      hs[l, ] = exp(seq(
        log(N^(-1 / (1 + p0)) * diff(range(zt[, l])) / 10),
        log(5 * N^(-1 / (1 + p0)) * diff(range(zt[, l]))),
        length.out = 20
      ))
    }
    
    parallel::clusterExport(cl, varlist = c("zt", "y", "KF", "kern", "p0", "N", "hs"), envir = environment())
    
    cv = unlist(parallel::parLapply(cl, 0:(20^p0 - 1), function(k) {
      h = numeric(p0)
      for (l in 1:p0) {
        kl = floor((k %% (20^l)) / (20^(l - 1))) + 1
        h[l] = hs[l, kl]
      }
      err = 0
      for (j in 1:N) {
        a = zt[j, ]
        if (p0 > 1) {
          mu1 = rowMeans(apply(zt[-j, , drop = FALSE], 1, function(zi) KF(zi - a, h) * (zi - a)))
          mu2 = matrix(rowMeans(apply(zt[-j, , drop = FALSE], 1, function(zi) KF(zi - a, h) * tcrossprod(zi - a))), ncol = p0)
        } else {
          mu1 = mean(sapply(zt[-j], function(zi) KF(zi - a, h) * (zi - a)))
          mu2 = mean(sapply(zt[-j], function(zi) KF(zi - a, h) * (zi - a)^2))
        }
        if (any(is.na(mu2)) || any(!is.finite(mu2))) return(Inf)
        if (inherits(try(solve(mu2), silent = TRUE), "try-error")) return(Inf)
        wc = t(mu1) %*% solve(mu2)
        w = apply(zt[-j, , drop = FALSE], 1, function(zi) KF(zi - a, h) * (1 - wc %*% (zi - a)))
        yj = weighted.mean(y[-j], w)
        err = err + (yj - y[j])^2 / N
      }
      err
    }))
    
    bwi = which.min(cv) - 1
    optns$bw1 = numeric(p0)
    for (l in 1:p0) {
      kl = floor((bwi %% (20^l)) / (20^(l - 1))) + 1
      optns$bw1[l] = hs[l, kl]
    }
  }
  
  if (is.null(optns$bw2)) optns$bw2 = optns$bw1
  if (!optns$regular) {
    optns$bw1 = optns$bw1[c(1, 2, 2)]
    optns$bw2 = optns$bw2[c(1, 2, 2)]
  }
  
  parallel::clusterExport(cl, varlist = c("z", "y", "p", "KF", "optns"), envir = environment())
  
  cm_list = parallel::parLapply(cl, 1:N, function(i) {
    a = z[i, ]
    if (p > 1) {
      mu1 = rowMeans(apply(z, 1, function(zi) KF(zi - a, optns$bw1) * (zi - a)))
      mu2 = matrix(rowMeans(apply(z, 1, function(zi) KF(zi - a, optns$bw1) * tcrossprod(zi - a))), ncol = p)
    } else {
      mu1 = mean(sapply(z, function(zi) KF(zi - a, optns$bw1) * (zi - a)))
      mu2 = mean(sapply(z, function(zi) KF(zi - a, optns$bw1) * (zi - a)^2))
    }
    if (inherits(try(solve(mu2), silent = TRUE), "try-error")) return(y[i])
    wc = t(mu1) %*% solve(mu2)
    w = apply(z, 1, function(zi) KF(zi - a, optns$bw1) * (1 - wc %*% (zi - a)))
    weighted.mean(y, w)
  })
  
  cm = unlist(cm_list)
  rcov = (y - cm)^2
  
  parallel::clusterExport(cl, varlist = c("rcov", "z", "y", "KF", "optns", "tp", "K", "p", "z0"), envir = environment())
  
  simulate_one_path = function(q) {
    path_q = rep(NA_real_, K)
    path_q[1] = z0[q, 1]
    for (i in 2:K) {
      a = c(path_q[i - 1], tp[i - 1])
      if (!optns$regular) a = c(a, tp[i])
      if (p > 1) {
        mu1 = rowMeans(apply(z, 1, function(zi) KF(zi - a, optns$bw1) * (zi - a)))
        mu2 = matrix(rowMeans(apply(z, 1, function(zi) KF(zi - a, optns$bw1) * tcrossprod(zi - a))), ncol = p)
      } else {
        mu1 = mean(sapply(z, function(zi) KF(zi - a, optns$bw1) * (zi - a)))
        mu2 = mean(sapply(z, function(zi) KF(zi - a, optns$bw1) * (zi - a)^2))
      }
      if (inherits(try(solve(mu2), silent = TRUE), "try-error")) break
      wc = t(mu1) %*% solve(mu2)
      w1 = apply(z, 1, function(zi) KF(zi - a, optns$bw1) * (1 - wc %*% (zi - a)))
      cmi = weighted.mean(y, w1)
      w2 = apply(z, 1, function(zi) KF(zi - a, optns$bw2))
      cci = max(weighted.mean(rcov, w2), 0)
      path_q[i] = optns$bm[q, i - 1] * sqrt(cci) + cmi
    }
    path_q
  }
  
  path_list = parallel::parLapply(cl, 1:optns$M, simulate_one_path)
  path = do.call(rbind, path_list)
  
  end_time = Sys.time()
  message("Total runtime: ", round(difftime(end_time, start_time, units = "secs"), 2), " seconds")
  
  res = list(
    path = path,
    Ly = Ly,
    Lt = Lt,
    z0 = z0,
    tp = tp,
    optns = optns
  )
  class(res) = "dm"
  res
}