#' @title Calculate survival probability
#' @description A unction to calculate survival probability for a single subject beyond a 
#'              single time point for a single posterior draw of the parameters
#' @param t scalar beyond which survival probabilities to be calculated
#' @param x design matrix of the fixed effects
#' @param d design matrix of the random effects
#' @param c design matrix of the survival sub-model
#' @param log_lambda matrix of posterior draws of log_lambda
#' @param log_nu matrix of posterior draws of log_nu
#' @param omega matrix of posterior draws of omega
#' @param eta matrix of posterior draws of eta
#' @param alpha matrix of posterior draws of alpha
#' @param B list of matrices for posterior draws of B
#' @param Q scalar for number of quadratures
#' @param wt vector of weights
#' @param pt vector of abscissa

 surv_prob_calc <- function(t, x, d, c, timeVar, log_lambda, log_nu, omega, 
                            eta, alpha, B, wt, pt, Q, bh){
   
   t_quad <- 0.5 * t * (1 + pt)
   
   c_quad <- apply(c, 2, function(i) rep(i, each = Q))
   
   #x_base <- x[1, , drop = FALSE]
   x_quad <- x[rep(1, Q), ]
   x_quad[, timeVar] <- t_quad
   
   #d_base <- d[1, , drop = FALSE]
   d_quad <- d[rep(1, Q), ]
   d_quad[, timeVar] <- t_quad
   
   h <- exp(log_lambda + log_nu + (exp(log_nu) - 1) * log(t_quad) + 
            c_quad %*% omega + 
            eta * (x_quad %*% alpha + d_quad %*% B))
   
   out <- exp(- 0.5 * t * sum(wt * h))
   return(out)
 }
 