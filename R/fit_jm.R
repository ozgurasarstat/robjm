#' @title Bayesian inference for joint models
#' 
#' @description Fits jont models with Normal and non-Normal distributions
#' 
#' @param fixed_long A two-sided formula for fixed effects of longitudinal sub-model
#' @param random_long A one-sided formula for random effects of longitudinal sub-model
#' @param fixed_surv A two-sided formula for survival time, event indicator and 
#'                  baseline covariates to be included in survival sub-model, 
#'                  e.g. \code{cbind(Time, death) ~ x}
#' @param data_long A data frame that includes the longitudinal data
#' @param data_surv A data frame that includes the survival data
#' @param id_long A character string that indicates the column name for the id column in \code{data_long}
#' @param id_surv A character string that indicates the column name for the id column in \code{data_surv}
#' @param model A character string for model identification; options are: 
#'             "nor_nor", "t_t_mod1", "t_t_mod2", "t_t_mod3", "nor_t_mod3",
#'             "t_t_tv", 
#' @param timeVar A character string for the column name of the time variable in \code{data_long}
#' @param bh A character string for baseline hazard specification, "weibull" for Weibull, 
#'           "spline" for b-spline, "piecewise" for piecewise constant
#' @param bh_nknots A numeric value for the number of knots for 
#'        b-spline and piecewie constant specifications for baseline hazard; default is 2; 
#'        not used for \code{bh = "weibull"}
#' @param spline_tv A list with two elements; first element is the name of the time variable, 
#'               and number of knots
#' @param Q A numeric value for the number of quadratures; default is 15              
#' @param priors A list of hyperparameters; theta, Omega, sigma_B, sigma_Z, beta (for tv), zeta, omega, eta. 
#'              See details below.
#' @param ... to be passed to the \code{stan} function 
#' 
#' @details This is a wrapper function for fitting mixed effects models. 
#'          Cauchy distribution is assumed as the prior for theta (QR decomposed alpha), 
#'          half-Cauchy for sigma_B (standard deviations of var-cov of B), 
#'          half-Cauchy for sigma_Z (standard deviation of error),
#'          Cauchy for beta (time-varying degree of freedom parameters),
#'          Cauchy for zeta (baseline hazard spline coefficients),
#'          Cauchy for omega (survival sub-model regression coefficients),
#'          Cauchy for eta (association parameter)
#'
#' @return Returns the output of the \code{stan} function 
#' @examples For examples, see the repository at https://github.com/ozgurasarstat/robjm-run                                              

fit_jm <- function(fixed_long, 
                   random_long, 
                   fixed_surv,
                   data_long, 
                   data_surv, 
                   id_long, 
                   id_surv,
                   model, 
                   timeVar,
                   deriv = NULL,
                   bh = "weibull",
                   bh_nknots = 2, #number of knots for baseline hazard 
                   spline_tv = list("time", 2), #spline for tv dof - same in fit_ld
                   Q = 15, 
                   priors = list(),
                   ...){

  ## be sure that distribution specifications are correct
  if(!(model %in% c("nor_nor", "t_t_mod1", "t_t_mod2", "t_t_mod3", "nor_t_mod3", "t_t_tv"))){
    stop("Model should be one of the followings: nor_nor, t_t_mod1, t_t_mod2, t_t_mod3, 
         nor_t_mod3, t_t_tv")
  }

  ## organise priors
  if(bh %in% c("spline", "piecewise")){
    if(model == "t_t_tv" & length(priors) != 8){
      priors_full <- list(alpha = 5, 
                          Omega = 2, 
                          sigma_B = 5, 
                          sigma_Z = 5,
                          beta = 4.6,
                          zeta = 5,
                          omega = 5,
                          eta = 5)
      for(i in 1:8){
        if(!(names(priors_full)[i] %in% names(priors))){
          priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
        }
      }
    }
    
    if(model != "t_t_tv" & length(priors) != 7){
      priors_full <- list(alpha = 5, 
                          Omega = 2, 
                          sigma_B = 5, 
                          sigma_Z = 5,
                          zeta = 5,
                          omega = 5,
                          eta = 5)
      for(i in 1:7){
        if(!(names(priors_full)[i] %in% names(priors))){
          priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
        }
      }
    }  
  }else if(bh == "weibull"){
      if(model == "t_t_tv" & length(priors) != 9){
        priors_full <- list(alpha = 5, 
                            Omega = 2, 
                            sigma_B = 5, 
                            sigma_Z = 5,
                            beta = 4.6,
                            log_lambda = 5,
                            log_nu = 5,
                            omega = 5,
                            eta = 5)
        for(i in 1:9){
          if(!(names(priors_full)[i] %in% names(priors))){
            priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
          }
        }
      }
      
      if(model != "t_t_tv" & length(priors) != 8){
        priors_full <- list(alpha = 5, 
                            Omega = 2, 
                            sigma_B = 5, 
                            sigma_Z = 5,
                            log_lambda = 5,
                            log_nu = 5,
                            omega = 5,
                            eta = 5)
        for(i in 1:8){
          if(!(names(priors_full)[i] %in% names(priors))){
            priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
          }
        }
      }  
    }
  
  
  ## be sure that id is: 1, 2, 3, ...
  data_long[, id_long] <- rep(1:length(unique(data_long[, id_long])), 
                              as.numeric(table(data_long[, id_long])))
  l_id <- data_long[, id_long]
  
  data_surv[, id_surv] <- rep(1:nrow(data_surv))
  s_id <- data_surv[, id_surv]
  
  ## x and y matrices
  x <- model.matrix(fixed_long, data_long)
  y <- model.frame(fixed_long, data_long)[, 1]
  
  ## create blok-diagonal random effects design matrix
  dmat <- model.matrix(random_long, data_long)
  id_dmat <- data.frame(l_id, dmat)
  id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
  d <- do.call(magic::adiag, id_dmat_list)
  
  ## total number of observations in the longitudinal data
  ## number of covariates in the x and d matrices
  ntot <- nrow(x) 
  p <- ncol(x)
  q <- ncol(id_dmat) - 1
  
  ## number of subjects
  ngroup <- length(unique(l_id))
  
  ## Gauss - Legendre weights and abscissas
  gl_quad <- statmod::gauss.quad(Q)
  wt <- gl_quad$weights
  pt <- gl_quad$nodes
  
  ## total number of observations for quadrature approx.
  ntot_quad <- ngroup * Q
  
  ## extract survival times and event indicator
  mf_surv <- model.frame(fixed_surv, data_surv)
  S <- mf_surv[, 1][, 1]
  E <- mf_surv[, 1][, 2]
  
  ## calculate times for hazard function for quadrature approx
  t_quad <- 0.5 * rep(S, each = Q) * (1 + rep(pt, ngroup))
  
  ## creat e and e_quad matrices
  if(bh %in% c("spline", "piecewise")){

    S_uncensored <- S[E == 1]
    knots <- quantile(S_uncensored, seq(0, 1, bh_nknots)[-c(1, (bh_nknots + 2))])
    
    if(bh == "piecewise"){
      e <- pw_mat(S, knots)
      ncol_e <- ncol(e)
    }else if(bh == "spline"){
      e <- bs(S, knots = knots)
      ncol_e <- ncol(e)
      attributes(e) <- NULL
      e <- matrix(e, ncol = ncol_e) 
    }
    
   if(bh == "piecewise"){
     e_quad <- pw_mat(t_quad, knots)
   }else if(bh == "spline"){
     e_quad <- bs(t_quad, knots = knots)
     attributes(e_quad) <- NULL
     e_quad <- matrix(e_quad, ncol = ncol_e) 
   }  
    
  }

  ## fixed effects for survival sub-model
  c <- model.matrix(fixed_surv, data_surv)[, -1, drop = FALSE]
  ncol_c <- ncol(c)
  c_quad <- apply(c, 2, function(i) rep(i, each = Q))
  
  ## x and d matrices at survival times and quadrature times
  data_long_base <- data_long[!duplicated(l_id), ]
  data_long_quad <- data_long_base[rep(1:ngroup, times = Q), ]
  
  data_long_base[, timeVar] <- S
  data_long_quad[, timeVar] <- t_quad
  
  x_T    <- model.matrix(fixed_long, data_long_base)
  x_quad <- model.matrix(fixed_long, data_long_quad)

  dmat_T <- model.matrix(random_long, data_long_base)
  id_dmat_T <- data.frame(s_id, dmat_T)
  id_dmat_list_T <- lapply(split(id_dmat_T[, -1], id_dmat_T[, 1]), as.matrix)
  d_T <- do.call(magic::adiag, id_dmat_list_T)
  
  dmat_quad <- model.matrix(random_long, data_long_quad)
  id_dmat_quad <- data.frame(rep(s_id, each = Q), dmat_quad)
  id_dmat_list_quad <- lapply(split(id_dmat_quad[, -1], id_dmat_quad[, 1]), as.matrix)
  d_quad <- do.call(magic::adiag, id_dmat_list_quad)
  
  ## prepare x and d matrices for the derivative
  if(!is.null(deriv)){
    
    deriv_fixed_formula  <- deriv$deriv_fixed_formula
    deriv_alpha_ind      <- deriv$deriv_alpha_ind
    deriv_random_formula <- deriv$deriv_random_formula
    deriv_B_ind          <- deriv$deriv_B_ind
    
    x_deriv_T <- model.matrix(deriv_fixed_formula, data_long_base)
    x_deriv_quad <- model.matrix(deriv_fixed_formula, data_long_quad)
    
    dmat_deriv_T <- model.matrix(deriv_random_formula, data_long_base)
    id_dmat_deriv_T <- data.frame(s_id, dmat_deriv_T)
    id_dmat_deriv_list_T <- lapply(split(id_dmat_deriv_T[, -1], id_dmat_deriv_T[, 1]), as.matrix)
    d_deriv_T <- do.call(magic::adiag, id_dmat_deriv_list_T)
    
    dmat_deriv_quad <- model.matrix(deriv_random_formula, data_long_quad)
    id_dmat_deriv_quad <- data.frame(rep(s_id, each = Q), dmat_deriv_quad)
    id_dmat_deriv_list_quad <- lapply(split(id_dmat_deriv_quad[, -1], id_dmat_deriv_quad[, 1]), as.matrix)
    d_deriv_quad <- do.call(magic::adiag, id_dmat_deriv_list_quad)
    
  }
  
  ## extend the weights for quadrature approx.
  wt_quad <- rep(wt, ngroup)
  
  ## prior hyperparameters
  if(bh %in% c("spline", "piecewise")){
    if(model != "t_t_tv"){
      priors_long <- unlist(priors)[1:4]
    }else{
      priors_long <- unlist(priors)[1:5]
    }
    priors_surv <- rev(unlist(priors))[1:3]
  }else if(bh == "weibull"){
    if(model != "t_t_tv"){
      priors_long <- unlist(priors)[1:4]
    }else{
      priors_long <- unlist(priors)[1:5]
    }
    priors_surv <- rev(unlist(priors))[1:4]
  }
  
  ## prepare data as a list to be passed to stan
  ## it is the shared version for 
  ## normal and t (mod1, mod2, mod3, tv), weibull, rate of change etc
  ## extra input will be added for the specific model when necessary
  
  data_stan <- list(ntot = ntot,
                    id = l_id, 
                    y = y, 
                    p = p,
                    q = q,
                    ngroup = ngroup, 
                    x = x, 
                    d = d,
                    priors_long = priors_long,
                    priors_surv = priors_surv,
                    Q = Q,
                    ntot_quad = ntot_quad,
                    S = S,
                    E = E,
                    ncol_c = ncol_c, 
                    c = c, 
                    c_quad = c_quad,
                    x_T = x_T,
                    x_quad = x_quad,
                    d_T = d_T,
                    d_quad = d_quad,
                    wt_quad = wt_quad,
                    t_quad = t_quad)
  
  if(model == "nor_nor"){
    
    if(bh %in% c("spline", "piecewise")){
      
      data_stan$ncol_e <- ncol_e
      data_stan$e <- e
      data_stan$e_quad <- e_quad
      
      res <- stan(model_code = nor_nor_jm, data = data_stan, ...)
      
    }else if(bh == "weibull"){
      
      if(!is.null(deriv)){
        
        data_stan$x_deriv_T <- x_deriv_T
        data_stan$x_deriv_quad <- x_deriv_quad
        data_stan$d_deriv_T <- d_deriv_T
        data_stan$d_deriv_quad <- d_deriv_quad
        data_stan$p_deriv <- length(deriv_alpha_ind)
        data_stan$q_deriv <- length(deriv_B_ind)
        data_stan$deriv_alpha_ind <- as.array(deriv_alpha_ind)
        data_stan$deriv_B_ind <- as.array(deriv_B_ind)
        
        res <- stan(model_code = nor_nor_jm_weibull_deriv, data = data_stan, ...)
        
      }else{
        res <- stan(model_code = nor_nor_jm_weibull, data = data_stan, ...)
      }
      
    }
    
  }
  
  if(model %in% c("t_t_mod1", "t_t_mod2", "t_t_mod3", "nor_t_mod3")){
    
  if(bh %in% c("spline", "piecewise")){
    
    data_stan$ncol_e <- ncol_e
    data_stan$e <- e
    data_stan$e_quad <- e_quad

    if(model == "t_t_mod1"){
      res <- stan(model_code = t_t_jm_mod1, data = data_stan, ...)
    } 
    if(model == "t_t_mod2"){
      res <- stan(model_code =t_t_jm_mod2, data = data_stan, ...)
    }
    if(model == "t_t_mod3"){
      res <- stan(model_code = t_t_jm_mod3, data = data_stan, ...)
    }
  }else if(bh == "weibull"){
    
    if(model == "t_t_mod1"){
      res <- stan(model_code = t_t_jm_mod1_weibull, data = data_stan, ...)
    } 
    if(model == "t_t_mod2"){
      res <- stan(model_code =t_t_jm_mod2_weibull, data = data_stan, ...)
    }
    if(model == "t_t_mod3"){
      res <- stan(model_code = t_t_jm_mod3_weibull, data = data_stan, ...)
    }
    if(model == "nor_t_mod3"){
      res <- stan(model_code = nor_t_jm_mod3_weibull, data = data_stan, ...)
    }
  }
    
 }

  if(model == "t_t_tv"){
    a <- splines::ns(data_long[, spline_tv[[1]]], df = (spline_tv[[2]] + 1))
    ncol_a <- ncol(a)
    attributes(a) <- NULL
    a <- matrix(a, ncol = ncol_a)

  if(bh %in% c("spline", "piecewise")){
    
    data_stan$ncol_e <- ncol_e
    data_stan$e <- e
    data_stan$e_quad <- e_quad
    data_stan$ncol_a <- ncol_a
    data_stan$a <- a
    
    res <- stan(model_code = t_t_tv_jm, data = data_stan, ...)
    
  }else if(bh == "weibull"){
    data_stan$ncol_a <- ncol_a
    data_stan$a <- a
    res <- stan(model_code = t_t_tv_jm_weibull, data = data_stan, ...)
  }
    
  }
  
  out <- list(fixed_long = fixed_long, 
              random_long = random_long, 
              fixed_surv = fixed_surv,
              data_long = data_long, 
              data_surv = data_surv, 
              id_long = id_long, 
              id_surv = id_surv,
              model = model, 
              timeVar = timeVar,
              bh = bh,
              bh_nknots = bh_nknots, 
              spline_tv = spline_tv, 
              Q = Q, 
              priors = priors,
              deriv = deriv,
              res = res
              )p
              
  return(out)  
  
}