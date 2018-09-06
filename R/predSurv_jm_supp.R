predSurv_jm_supp <- function(batch_data, 
                             forecast, 
                             B_control, 
                             model,
                             bh,
                             id_long,
                             id_surv,
                             fixed_long,
                             fixed_surv,
                             random_long,
                             timeVar,
                             Q,
                             wt,
                             pt){
  
  ## extract the chains
  alpha      <- rstan::extract(object$res)$alpha
  Sigma_long <- rstan::extract(object$res)$Sigma
  M          <- nrow(alpha)
  Sigma      <- lapply(1:M, function(i) Sigma_long[i, ,])
  sigma_Z    <- matrix(rstan::extract(object$res)$sigma_Z)
  log_lambda <- matrix(rstan::extract(object$res)$log_lambda)
  log_nu     <- matrix(rstan::extract(object$res)$log_nu)
  omega      <- rstan::extract(object$res)$omega
  eta        <- matrix(rstan::extract(object$res)$eta)
  
  if(object$model == "t_t_mod3"){
    phi   <- matrix(rstan::extract(object$res)$phi)
    delta <- matrix(rstan::extract(object$res)$delta)
  }
  
  ##
  ## first predict the random effects for the new subjects
  ##

  data_surv <- batch_data[!duplicated(batch_data[, id_long]), ]
  
  ngroup <- length(unique(batch_data[, id_long]))
  nobs   <- as.numeric(table(batch_data[, id_long]))
  batch_data[, id_long] <- rep(1:ngroup, nobs)
  
  s_id_orig <- data_surv[, id_surv]
  data_surv[, id_surv] <- 1:ngroup
  
  l_id <- batch_data[, id_long]
  s_id <- data_surv[, id_surv]
  
  x <- model.matrix(fixed_long, batch_data)
  y <- model.frame(fixed_long, batch_data)[, 1]
  
  dmat <- model.matrix(random_long, batch_data)
  id_dmat <- data.frame(l_id, dmat)
  id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
  d <- do.call(magic::adiag, id_dmat_list)
  
  ## total number of observations in the longitudinal data
  ## number of covariates in the x and d matrices
  ntot <- nrow(x) 
  p <- ncol(x)
  q <- ncol(id_dmat) - 1
  
  ## extract survival times and event indicator
  mf_surv <- model.frame(fixed_surv, data_surv)
  S <- mf_surv[, 1][, 1]
  E <- mf_surv[, 1][, 2]
  
  ## calculate times for hazard function for quadrature approx
  t_quad <- 0.5 * rep(S, each = Q) * (1 + rep(pt, ngroup))
  
  ## total number of observations for quadrature approx.
  ntot_quad <- ngroup * Q
  
  ## fixed effects for survival sub-model
  c <- model.matrix(fixed_surv, data_surv)[, -1, drop = FALSE]
  ncol_c <- ncol(c)
  c_quad <- apply(c, 2, function(i) rep(i, each = Q))
  
  ## x matrix for log survival density
  x_T <- x[!duplicated(l_id), , drop = FALSE]
  x_T[, timeVar] <- S
  
  x_quad <- x_T[rep(1:ngroup, times = Q), ]
  x_quad[, timeVar] <- t_quad
  
  ## d matrix for log survival density
  dmat_T <- dmat[!duplicated(l_id), , drop = FALSE]
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
  
  if(model == "nor_nor" & bh == "weibull"){
    
    mod <- stan_model(model_code = new_rand_eff_nor_nor_jm_weibull, auto_write = TRUE)
    
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
                         S = as.array(S),
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
      
      B_res <- sampling(mod, 
                        data = data_nor_nor, 
                        iter = B_control$iter, 
                        warmup = B_control$warmup,
                        chains = B_control$chains,
                        control = list(adapt_delta = B_control$adapt_delta, 
                                       max_treedepth = B_control$max_treedepth)
                        )
      B_sampled[[i]] <- matrix(rstan::summary(B_res)$summary[1:(ngroup*q), "50%"], ncol = q, byrow = T)
    }  
  }
  
  if(model == "t_t_mod3" & bh == "weibull"){
    
    mod <- stan_model(model_code = new_rand_eff_t_t_mod3_jm_weibull, auto_write = TRUE)
    
    data_t_t_mod3 <- list(ntot = ntot,
                          id = l_id, 
                          y = y, 
                          p = p,
                          q = q,
                          ngroup = ngroup, 
                          x = x, 
                          d = d,
                          Q = Q,
                          ntot_quad = ntot_quad,
                          S = as.array(S),
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
      data_t_t_mod3$alpha      <- as.array(alpha[i, ])
      data_t_t_mod3$Sigma      <- Sigma[[i]]
      data_t_t_mod3$sigma_Z    <- sigma_Z[i, ]
      data_t_t_mod3$log_lambda <- log_lambda[i, ]
      data_t_t_mod3$log_nu     <- log_nu[i, ]
      data_t_t_mod3$omega      <- as.array(omega[i, ])
      data_t_t_mod3$eta        <- eta[i, ]
      data_t_t_mod3$phi        <- phi[i, ]
      data_t_t_mod3$delta      <- delta[i, ]
      
      B_res <- sampling(mod, 
                        data = data_t_t_mod3, 
                        iter = B_control$iter, 
                        warmup = B_control$warmup,
                        chains = B_control$chains,
                        control = list(adapt_delta = B_control$adapt_delta, 
                                       max_treedepth = B_control$max_treedepth)
                        )
      
      B_sampled[[i]] <- matrix(rstan::summary(B_res)$summary[1:(ngroup*q), "50%"], ncol = q, byrow = T)
    }
  }
  
  ##
  ## the calculate the survival probabilities by pluggin in the ratio 
  ##
  
  ft <- list()
  for(i in 1:ngroup){
    ft[[i]] <- seq(S[i], (S[i] + forecast$h), length.out = forecast$n)
  }
  
  x_base <- x[!duplicated(l_id), , drop = FALSE]
  d_base <- dmat[!duplicated(l_id), , drop = FALSE]
  
  ft_probs <- list()
  
  for(i in 1:ngroup){
    
    ft_i <- ft[[i]]
    
    ft_probs_i <- list()
    ft_probs_i[[1]] <- rep(1, M)
    
    for(j in 2:forecast$n){
      
      ft_probs_i_k <- c()
      
      for(k in 1:M){
        prob_upper <- surv_prob_calc(t = ft_i[j], 
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
        prob_lower <- surv_prob_calc(t = ft_i[1], 
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
    ft_i <- ft[[i]]
    out_i <- data.frame(id = rep(s_id_orig[i], forecast$n),
                        stime = rep(S[i], forecast$n),
                        event = rep(E[i], forecast$n),
                        time = ft_i, 
                        do.call(rbind, lapply(ft_probs[[i]], prob_summary)))
    names(out_i)[5:8] <- c("2.5%", "mean", "median", "97.5%")
    out[[i]] <- out_i
  }
  out <- do.call(rbind, out)
  return(out)
  
}