#' @title Predict survival probabilities
#' 
#' @description A function to predict survival probabilities
#' 
#' @param object A fitted object from \code{fit_jm}
#' @param newdata 
#' 

predSurv_jm <- function(object, newdata, l, h, inc, ...){

  ##
  ## first predict the random effects for the new subjects
  ##
  
  ###### normal-normal model with Weibull baseline hazard
  
  ## extract the chains
  alpha      <- extract(object$res)$alpha
  Sigma_long <- extract(object$res)$Sigma
  M          <- nrow(alpha)
  Sigma      <- lapply(1:M, function(i) Sigma_long[i, ,])
  sigma_Z    <- matrix(extract(object$res)$sigma_Z)
  log_lambda <- matrix(extract(object$res)$log_lambda)
  log_nu     <- matrix(extract(object$res)$log_nu)
  omega      <- extract(object$res)$omega
  eta        <- matrix(extract(object$res)$eta)
  
  ## create covariate matrices
  id_long <- object$id_long
  id_surv <- object$id_surv
  
  ngroup <- length(unique(newdata[, id_long]))
  nobs   <- as.numeric(table(newdata[, id_long]))
  #newdata[, id_long] <- rep(1:ngroup, nobs)
  
  data_surv <- newdata[!duplicated(newdata[, id_long]), ]
  
  l_id <- newdata[, id_long]
  s_id <- data_surv[, id_surv]
  
  x <- model.matrix(object$fixed_long, newdata)
  y <- model.frame(object$fixed_long, newdata)[, 1]
  
  dmat <- model.matrix(object$random_long, newdata)
  id_dmat <- data.frame(l_id, dmat)
  id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
  d <- do.call(magic::adiag, id_dmat_list)
  
  ## total number of observations in the longitudinal data
  ## number of covariates in the x and d matrices
  ntot <- nrow(x) 
  p <- ncol(x)
  q <- ncol(id_dmat) - 1
  
  ## extract survival times and event indicator
  S <- model.frame(object$fixed_surv, data_surv)[, 1][, 1]
  
  ## Gauss - Legendre weights and abscissas
  Q <- object$Q
  gl_quad <- statmod::gauss.quad(Q)
  wt <- gl_quad$weights
  pt <- gl_quad$nodes
  
  ## calculate times for hazard function for quadrature approx
  t_quad <- 0.5 * rep(S, each = Q) * (1 + rep(pt, ngroup))
  
  ## total number of observations for quadrature approx.
  ntot_quad <- ngroup * Q
  
  ## fixed effects for survival sub-model
  c <- model.matrix(object$fixed_surv, data_surv)[, -1, drop = FALSE]
  ncol_c <- ncol(c)
  c_quad <- apply(c, 2, function(i) rep(i, each = Q))
  
  ## x matrix for log survival density
  timeVar <- object$timeVar
  x_T <- x[!duplicated(l_id), ]
  x_T[, timeVar] <- S
  
  x_quad <- x_T[rep(1:ngroup, times = Q), ]
  x_quad[, timeVar] <- t_quad
  
  ## d matrix for log survival density
  dmat_T <- dmat[!duplicated(l_id), ]
  dmat_T[, timeVar] <- S
  id_dmat_T <- data.frame(s_id, dmat_T)
  id_dmat_list_T <- lapply(split(id_dmat_T[, -1], id_dmat_T[, 1]), as.matrix)
  d_T <- do.call(magic::adiag, id_dmat_list_T)
  
  dmat_quad <- dmat_T[rep(1:ngroup, times = Q), ]
  dmat_quad[, timeVar] <- t_quad
  id_dmat_quad <- data.frame(rep(s_id, each = Q), dmat_quad)
  id_dmat_list_quad <- lapply(split(id_dmat_quad[, -1], id_dmat_quad[, 1]), as.matrix)
  d_quad <- do.call(magic::adiag, id_dmat_list_quad)
  
  ## extend the weights for quadrature approx.
  wt_quad <- rep(wt, ngroup)
  
  data_nor_nor <- list(ntot = ntot,
                       id = l_id, 
                       y = y, 
                       p = p,
                       q = q,
                       ngroup = ngroup, 
                       x = x, 
                       d = d,
                       Q = Q,
                       ntot_quad = ntot_quad,
                       S = S,
                       ncol_c = ncol_c, 
                       c = c, 
                       c_quad = c_quad,
                       x_T = x_T,
                       x_quad = x_quad,
                       d_T = d_T,
                       d_quad = d_quad,
                       wt_quad = wt_quad
                       )
  B_sampled <- list()
  
  for(i in 1:M){
    data_nor_nor$alpha      <- as.array(alpha[i, ])
    data_nor_nor$Sigma      <- Sigma[[i]]
    data_nor_nor$sigma_Z    <- sigma_Z[i, ]
    data_nor_nor$log_lambda <- log_lambda[i, ]
    data_nor_nor$log_nu     <- log_nu[i, ]
    data_nor_nor$omega      <- as.array(omega[i, ])
    data_nor_nor$eta        <- eta[i, ]
    
    res <- stan(model_code = new_rand_eff_nor_nor_jm_weibull, 
                data = data_nor_nor, 
                iter = 1, 
                chains = 1,
                warmup = 0,
                control = list(adapt_delta = 0.9999, max_treedepth = 15)
                )
    
    B_sampled[[i]] <- matrix(extract(res)$B, ncol = q, byrow = T)
    
  }
  
  ##
  ## the calculate the survival probabilities by pluggin in the ratio 
  ##
  
  ft <- seq(l, h, inc)
  n_ft <- length(ft)

  x_base <- x[!duplicated(l_id), , drop = FALSE]
  d_base <- dmat[!duplicated(l_id), , drop = FALSE]

  ft_probs <- list()
  
  for(i in 1:ngroup){

    x_i <- x[nobs_cumsum[i]:nobs_cumsum[i+1], , drop = FALSE]
    d_i <- dmat[nobs_cumsum[i]:nobs_cumsum[i+1], , drop = FALSE]
    c_i <- c[i, , drop = FALSE]
    
    ft_probs_i <- list()
    ft_probs_i[[1]] <- rep(1, M)
    for(j in 2:n_ft){
      ft_probs_i_k <- c()
      for(k in 1:M){
        prob_upper <- surv_prob_calc(t = ft[j], 
                                     x = x_base[i, , drop = FALSE], 
                                     d = d_base[i, , drop = FALSE], 
                                     c = c[i, , drop = FALSE],
                                     timeVar = timeVar, 
                                     log_lambda = log_lambda[k, ],
                                     log_nu = log_nu[k, ],
                                     omega = omega[k, ],
                                     eta = eta[k, ],
                                     alpha = alpha[k, ],
                                     B = B_sampled[[k]][i, ],
                                     wt = wt, 
                                     pt = pt,
                                     Q = Q,
                                     bh = "weibull")
        prob_lower <- surv_prob_calc(t = ft[1], 
                                     x = x_base[i, , drop = FALSE], 
                                     d = d_base[i, , drop = FALSE], 
                                     c = c[i, , drop = FALSE],
                                     timeVar = timeVar, 
                                     log_lambda = log_lambda[k, ],
                                     log_nu = log_nu[k, ],
                                     omega = omega[k, ],
                                     eta = eta[k, ],
                                     alpha = alpha[k, ],
                                     B = B_sampled[[k]][i, ],
                                     wt = wt, 
                                     pt = pt,
                                     Q = Q,
                                     bh = "weibull")
        ft_probs_i_k <- c(ft_probs_i_k, prob_upper/prob_lower)
      }
      ft_probs_i[[j]] <- ft_probs_i_k
    }
        
    ft_probs[[i]] <- ft_probs_i 
    
  }
  
  out <- list()
  for(i in 1:ngroup){
    out_i <- cbind(rep(s_id[i], n_ft), ft, do.call(rbind, lapply(ft_probs[[i]], prob_summary)))
    out[[i]] <- out_i
  }
  out <- do.call(rbind, out)
  colnames(out) <- c("id", "time", "2.5%", "mean", "mean", "97.5%")
  return(out)
  
}
