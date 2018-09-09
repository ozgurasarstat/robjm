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

 surv_prob_calc <- function(t, ft_batch_data_quad_i, c_quad_i, 
                            timeVar, log_lambda, log_nu, omega, 
                            eta, alpha, B, wt, pt, Q, bh, deriv,
                            fixed_long, random_long){
   
     t_quad <- 0.5 * t * (1 + pt)
     
     ft_batch_data_quad_i[, timeVar] <- t_quad
     
     x_quad <- model.matrix(fixed_long, ft_batch_data_quad_i)
     d_quad <- model.matrix(random_long, ft_batch_data_quad_i)
     
     if(!is.null(deriv)){
       x_quad_deriv <- model.matrix(deriv$deriv_fixed_formula, ft_batch_data_quad_i)
       d_quad_deriv <- model.matrix(deriv$deriv_random_formula, ft_batch_data_quad_i)
       
       alpha_deriv <- alpha[deriv$deriv_alpha_ind]
       B_deriv     <- B[deriv$deriv_B_ind]
     }
     
     if(is.null(deriv)){
       h <- exp(log_lambda + log_nu + (exp(log_nu) - 1) * log(t_quad) + 
                  c_quad_i %*% omega + 
                  eta * (x_quad %*% alpha + d_quad %*% B))       
     }else{
       h <- exp(log_lambda + log_nu + (exp(log_nu) - 1) * log(t_quad) + 
                  c_quad_i %*% omega + 
                  eta[1] * (x_quad %*% alpha + d_quad %*% B) + 
                  eta[2] * (x_quad_deriv %*% alpha_deriv + d_quad_deriv %*% B_deriv))
     }

     
     out <- exp(- 0.5 * t * sum(wt * h))
     return(out)
 }
 