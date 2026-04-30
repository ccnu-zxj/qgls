#' @noRd
qgls_single_base_Ok <- function(input_x, p_naive, n, q_func, d_func, c_F) {
  k <- length(input_x)
  r_i <- p_naive * (n + 1)
  p <- (r_i - c_F) / (n - 2 * c_F + 1)
  
  q_val <- q_func(p)
  f_val <- d_func(q_val)
  qf_val <- q_val * f_val
  
  p_pad <- c(0, p, 1)
  f_pad <- c(0, f_val, 0)
  qf_pad <- c(0, qf_val, 0)
  
  dp <- diff(p_pad)
  # 【代数防御层】：处理极值节点的浮点下溢
  diff_f_p <- suppressWarnings(diff(f_pad) / dp)
  diff_f_p[is.nan(diff_f_p) | is.infinite(diff_f_p)] <- 0
  diff_qf_p <- suppressWarnings(diff(qf_pad) / dp)
  diff_qf_p[is.nan(diff_qf_p) | is.infinite(diff_qf_p)] <- 0
  
  u_w <- f_val * (diff_f_p[1:k] - diff_f_p[2:(k+1)])
  u_z <- f_val * (diff_qf_p[1:k] - diff_qf_p[2:(k+1)])
  
  I_11 <- sum(u_w)
  I_22 <- sum(u_z * q_val)
  I_12 <- sum(u_z)
  det_I <- I_11 * I_22 - I_12^2
  
  w <- (I_22 * u_w - I_12 * u_z) / det_I
  v <- (I_11 * u_z - I_12 * u_w) / det_I
  
  loc_est <- sum(w * input_x)
  scl_est <- sum(v * input_x)
  
  residual <- input_x - (loc_est + scl_est * q_val)
  g_pad <- c(0, f_val * residual, 0)
  
  diff_g_p <- suppressWarnings(diff(g_pad)^2 / dp)
  diff_g_p[is.nan(diff_g_p) | is.infinite(diff_g_p)] <- 0
  J_n_raw <- sum(diff_g_p)
  J_n <- (n / scl_est^2) * J_n_raw
  
  list(location = loc_est, scale = scl_est, J_n = J_n)
}

#' Adaptive Quantile Generalized Least Squares
#' @export
qgls <- function(x, p, n) {
  if (length(x) != length(p)) stop('Dimension mismatch.')
  if (n <= 0) stop('Sample size must be positive.')
  if (any(p <= 0 | p >= 1)) stop('Probabilities must be in (0,1).')
  if (any(diff(p) <= 0)) stop('Probabilities must be strictly increasing.')
  if (any(diff(x) < 0)) stop('Quantiles must be monotonically non-decreasing.')
  if (length(x) < 3) stop('Degrees of freedom insufficient (k >= 3 required).')
  
  euler_gamma <- -digamma(1)
  bases <- list(
    Normal = list(is_log=FALSE, q=qnorm, d=dnorm, c_F=0.375, to_mean=function(l,s) l, to_sd=function(l,s) s),
    Logistic = list(is_log=FALSE, q=qlogis, d=dlogis, c_F=0.314, to_mean=function(l,s) l, to_sd=function(l,s) s*pi/sqrt(3)),
    Laplace = list(is_log=FALSE, q=function(p) ifelse(p<0.5, log(2*p), -log(2*(1-p))), d=function(x) 0.5*exp(-abs(x)), c_F=0.0, to_mean=function(l,s) l, to_sd=function(l,s) s*sqrt(2)),
    Gumbel = list(is_log=FALSE, q=function(p) -log(-log(p)), d=function(x) exp(-(x+exp(-x))), c_F=0.44, to_mean=function(l,s) l+s*euler_gamma, to_sd=function(l,s) s*pi/sqrt(6)),
    LogNormal = list(is_log=TRUE, q=qnorm, d=dnorm, c_F=0.375, to_mean=function(l,s) exp(l+s^2/2), to_sd=function(l,s) sqrt((exp(s^2)-1)*exp(2*l+s^2)))
  )
  
  results <- list()
  for (name in names(bases)) {
    b <- bases[[name]]
    if (b$is_log && any(x <= 0)) next
    input_x <- if (b$is_log) log(x) else x
    fit <- tryCatch(qgls_single_base_Ok(input_x, p, n, b$q, b$d, b$c_F), error = function(e) NULL)
    if (!is.null(fit)) {
      fit$True_Mean <- b$to_mean(fit$location, fit$scale)
      fit$True_SD <- b$to_sd(fit$location, fit$scale)
      results[[name]] <- fit
    }
  }
  
  if (length(results) == 0) stop('System failed: Algebraic singularity.')
  J_values <- sapply(results, function(r) r$J_n)
  opt_name <- names(which.min(J_values))
  opt_fit <- results[[opt_name]]
  
  out <- list(Optimal_Distribution = opt_name, Estimated_Mean = opt_fit$True_Mean, Estimated_SD = opt_fit$True_SD, Goodness_of_Fit = opt_fit$J_n, All_Distances = J_values)
  class(out) <- 'qgls'
  return(out)
}

#' @export
print.qgls <- function(x, ...) {
  cat('\nAdaptive Quantile Generalized Least Squares (Q-GLS)\n')
  cat('===================================================\n')
  cat(sprintf('Optimal Base Measure : %s\n', x$Optimal_Distribution))
  cat(sprintf('Estimated Mean       : %f\n', x$Estimated_Mean))
  cat(sprintf('Estimated SD         : %f\n', x$Estimated_SD))
  cat(sprintf('Generalized Distance : %f\n\n', x$Goodness_of_Fit))
  invisible(x)
}
