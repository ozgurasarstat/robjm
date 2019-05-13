predSurv_jm_supp <- function(object, 
                             batch_data, 
                             forecast, 
                             B_control, 
                             Q,
                             wt,
                             pt,
                             last_time,
                             lm_time,
                             probs,
                             return_bsamples){
  
  ## make sure that first row for everyone is prob of 1 at stime
  forecast$n <- forecast$n + 1
  
  ## extract info to be used
  model <- object$model
  bh <- object$bh
  id_long <- object$id_long
  id_surv <- object$id_surv
  fixed_long <- object$fixed_long
  fixed_surv <- object$fixed_surv
  random_long <- object$random_long
  timeVar <- object$timeVar
  deriv <- object$deriv
  
  ## extract the chains
  alpha      <- rstan::extract(object$res)$alpha
  Sigma_long <- rstan::extract(object$res)$Sigma
  M          <- nrow(alpha)
  Sigma      <- lapply(1:M, function(i) Sigma_long[i, ,])
  sigma_Z    <- matrix(rstan::extract(object$res)$sigma_Z)
  log_lambda <- matrix(rstan::extract(object$res)$log_lambda)
  log_nu     <- matrix(rstan::extract(object$res)$log_nu)
  omega      <- rstan::extract(object$res)$omega
  if(is.null(deriv)){
    eta <- matrix(rstan::extract(object$res)$eta)
  }else{
    eta1 <- matrix(rstan::extract(object$res)$eta1)
    eta2 <- matrix(rstan::extract(object$res)$eta2)
  }
  
  if(object$model != "nor_nor"){
    if(object$model == "nor_t_mod3"){
      delta <- matrix(rstan::extract(object$res)$delta)     
    }else if(object$model == "t_nor_mod3"){
      phi   <- matrix(rstan::extract(object$res)$phi)
    }else if(object$model == "t_t_mod3"){
      phi   <- matrix(rstan::extract(object$res)$phi)
      delta <- matrix(rstan::extract(object$res)$delta)     
    }
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
  
  # dmat <- model.matrix(random_long, batch_data)
  # id_dmat <- data.frame(l_id, dmat)
  # id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
  # d <- do.call(magic::adiag, id_dmat_list)

  d <- model.matrix(random_long, batch_data)
  
  ## total number of observations in the longitudinal data
  ## number of covariates in the x and d matrices
  ntot <- nrow(x) 
  p <- ncol(x)
  q <- ncol(d)
  
  ## extract survival times and event indicator
  if(last_time == "surv_time"){
    mf_surv <- model.frame(fixed_surv, data_surv)
    S <- mf_surv[, 1][, 1]
    #E <- mf_surv[, 1][, 2]    
  }else if(last_time == "landmark"){
    
    if(is.null(lm_time)) stop("Define landmark time")
    #if(S < lm_time) stop("Survival time(s) is(are) less than lm_time")
    
    if(length(lm_time) == 1){
      S <- rep(lm_time, ngroup)
    }else{
      S <- lm_time
    }
    
    if(forecast$n > 2){
      forecast$n <- 2
      warning("n in forecast is set to 1")
    }
    
  }else{
    S <- batch_data[!duplicated(batch_data[, id_long], fromLast = TRUE), last_time]
  }
  
  ## calculate times for hazard function for quadrature approx
  t_quad <- 0.5 * rep(S, each = Q) * (1 + rep(pt, ngroup))
  
  ## total number of observations for quadrature approx.
  ntot_quad <- ngroup * Q
  
  ## fixed effects for survival sub-model
  c <- model.matrix(fixed_surv, data_surv)[, -1, drop = FALSE]
  ncol_c <- ncol(c)
  c_quad <- apply(c, 2, function(i) rep(i, each = Q))
  
  ## x matrix for log survival density
  
  batch_data_base <- batch_data[!duplicated(l_id), ]
  batch_data_quad <- batch_data_base[rep(1:ngroup, times = Q), ]
  
  batch_data_base[, timeVar] <- S
  batch_data_quad[, timeVar] <- t_quad
  
  x_T    <- model.matrix(fixed_long, batch_data_base)
  x_quad <- model.matrix(fixed_long, batch_data_quad)

  ## d matrix for log survival density
  # dmat_T <- model.matrix(random_long, batch_data_base)
  # id_dmat_T <- data.frame(s_id, dmat_T)
  # id_dmat_list_T <- lapply(split(id_dmat_T[, -1], id_dmat_T[, 1]), as.matrix)
  # d_T <- do.call(magic::adiag, id_dmat_list_T)
  
  d_T <- model.matrix(random_long, batch_data_base)
  
  
  # dmat_quad <- model.matrix(random_long, batch_data_quad)
  # id_dmat_quad <- data.frame(rep(s_id, each = Q), dmat_quad)
  # id_dmat_list_quad <- lapply(split(id_dmat_quad[, -1], id_dmat_quad[, 1]), as.matrix)
  # d_quad <- do.call(magic::adiag, id_dmat_list_quad)
  
  d_quad <- model.matrix(random_long, batch_data_quad)
  
  ## prepare x and d matrices for the derivative
  if(!is.null(deriv)){
    
    deriv_fixed_formula  <- deriv$deriv_fixed_formula
    deriv_alpha_ind      <- deriv$deriv_alpha_ind
    deriv_random_formula <- deriv$deriv_random_formula
    deriv_B_ind          <- deriv$deriv_B_ind
    
    x_deriv_T <- model.matrix(deriv_fixed_formula, batch_data_base)
    x_deriv_quad <- model.matrix(deriv_fixed_formula, batch_data_quad)
    
    # dmat_deriv_T <- model.matrix(deriv_random_formula, batch_data_base)
    # id_dmat_deriv_T <- data.frame(s_id, dmat_deriv_T)
    # id_dmat_deriv_list_T <- lapply(split(id_dmat_deriv_T[, -1], id_dmat_deriv_T[, 1]), as.matrix)
    # d_deriv_T <- do.call(magic::adiag, id_dmat_deriv_list_T)
    
    d_deriv_T <- model.matrix(deriv_random_formula, batch_data_base)
    
    # dmat_deriv_quad <- model.matrix(deriv_random_formula, batch_data_quad)
    # id_dmat_deriv_quad <- data.frame(rep(s_id, each = Q), dmat_deriv_quad)
    # id_dmat_deriv_list_quad <- lapply(split(id_dmat_deriv_quad[, -1], id_dmat_deriv_quad[, 1]), as.matrix)
    # d_deriv_quad <- do.call(magic::adiag, id_dmat_deriv_list_quad)
    
    d_deriv_quad <- model.matrix(deriv_random_formula, batch_data_quad)
    
  }
  
  ## extend the weights for quadrature approx.
  wt_quad <- rep(wt, ngroup)
  
  # prepare a matrix of indices to select rows of d in for loop in stan
  cumsum_nrepeat <- cumsum(nobs)
  d_ind <- cbind(c(1, (cumsum_nrepeat[-ngroup] + 1)), cumsum_nrepeat)
  
  # prepare a matrix of indices for matrices for quadratures
  Q_ind <- cbind((0:(ngroup-1))*Q+1, (1:ngroup)*Q)
  
  ## nor nor 
  if(model == "nor_nor" & bh == "weibull"){
    
    if(is.null(deriv)){
      mod <- rstan::stan_model(model_code = new_rand_eff_nor_nor_jm_weibull, auto_write = TRUE)
    }else{
      mod <- rstan::stan_model(model_code = new_rand_eff_nor_nor_jm_weibull_deriv, auto_write = TRUE)
    }
    
    data_nor_nor <- list(ntot = ntot,
                         id = as.array(l_id), 
                         y = as.array(y), 
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
                         wt_quad = wt_quad,
                         d_ind = d_ind,
                         Q_ind = Q_ind)
    
    if(!is.null(deriv)){
      data_nor_nor$x_deriv_T <- x_deriv_T
      data_nor_nor$x_deriv_quad <- x_deriv_quad
      data_nor_nor$d_deriv_T <- d_deriv_T
      data_nor_nor$d_deriv_quad <- d_deriv_quad
      data_nor_nor$p_deriv <- length(deriv_alpha_ind)
      data_nor_nor$q_deriv <- length(deriv_B_ind)
      data_nor_nor$deriv_alpha_ind <- as.array(deriv_alpha_ind)
      data_nor_nor$deriv_B_ind <- as.array(deriv_B_ind)
    }
    
    B_sampled <- list()
    
    for(i in 1:M){
      data_nor_nor$alpha      <- as.array(alpha[i, ])
      data_nor_nor$Sigma      <- Sigma[[i]]
      data_nor_nor$sigma_Z    <- sigma_Z[i, ]
      data_nor_nor$log_lambda <- log_lambda[i, ]
      data_nor_nor$log_nu     <- log_nu[i, ]
      data_nor_nor$omega      <- as.array(omega[i, ])
      
      if(is.null(deriv)){
        data_nor_nor$eta      <- eta[i, ]
      }else{
        data_nor_nor$eta1     <- eta1[i, ]
        data_nor_nor$eta2     <- eta2[i, ]
      }

      B_res <- rstan::sampling(mod, 
                        data = data_nor_nor, 
                        iter = B_control$iter, 
                        warmup = B_control$warmup,
                        chains = B_control$chains,
                        cores = B_control$cores,
                        init = B_control$init,
                        control = list(adapt_delta = B_control$adapt_delta, 
                                       max_treedepth = B_control$max_treedepth),
                        pars = c("B")
                        )
      #B_sampled[[i]] <- matrix(rstan::summary(B_res)$summary[1:(ngroup*q), "50%"], ncol = q, byrow = T)
      B_sampled[[i]] <- rstan::extract(B_res)[["B"]] %>% subsample_B(nsel_b = B_control$nsel_b)
    }  
  }

  ## nor-t mod3
  if(model == "nor_t_mod3" & bh == "weibull"){
    
    if(is.null(deriv)){
      mod <- rstan::stan_model(model_code = new_rand_eff_nor_t_mod3_jm_weibull, auto_write = TRUE)
    }else{
      mod <- rstan::stan_model(model_code = new_rand_eff_nor_t_mod3_jm_weibull_deriv, auto_write = TRUE)
    }
    
    data_nor_t_mod3 <- list(ntot = ntot,
                          id = as.array(l_id), 
                          y = as.array(y), 
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
                          wt_quad = wt_quad,
                          d_ind = d_ind,
                          Q_ind = Q_ind)
    
    if(!is.null(deriv)){
      data_nor_t_mod3$x_deriv_T <- x_deriv_T
      data_nor_t_mod3$x_deriv_quad <- x_deriv_quad
      data_nor_t_mod3$d_deriv_T <- d_deriv_T
      data_nor_t_mod3$d_deriv_quad <- d_deriv_quad
      data_nor_t_mod3$p_deriv <- length(deriv_alpha_ind)
      data_nor_t_mod3$q_deriv <- length(deriv_B_ind)
      data_nor_t_mod3$deriv_alpha_ind <- as.array(deriv_alpha_ind)
      data_nor_t_mod3$deriv_B_ind <- as.array(deriv_B_ind)
    }
    
    B_sampled <- list()
    
    for(i in 1:M){
      data_nor_t_mod3$alpha      <- as.array(alpha[i, ])
      data_nor_t_mod3$Sigma      <- Sigma[[i]]
      data_nor_t_mod3$sigma_Z    <- sigma_Z[i, ]
      data_nor_t_mod3$log_lambda <- log_lambda[i, ]
      data_nor_t_mod3$log_nu     <- log_nu[i, ]
      data_nor_t_mod3$omega      <- as.array(omega[i, ])
      
      if(is.null(deriv)){
        data_nor_t_mod3$eta <- eta[i, ]
      }else{
        data_nor_t_mod3$eta1 <- eta1[i, ]
        data_nor_t_mod3$eta2 <- eta2[i, ]
      }
      
      data_nor_t_mod3$delta      <- delta[i, ]
      
      B_res <- rstan::sampling(mod, 
                               data = data_nor_t_mod3, 
                               iter = B_control$iter, 
                               warmup = B_control$warmup,
                               chains = B_control$chains,
                               cores = B_control$cores,
                               init = B_control$init,
                               control = list(adapt_delta = B_control$adapt_delta, 
                                              max_treedepth = B_control$max_treedepth)
      )
      
      #B_sampled[[i]] <- matrix(rstan::summary(B_res)$summary[1:(ngroup*q), "50%"], ncol = q, byrow = T)
      B_sampled[[i]] <- rstan::extract(B_res)[["B"]] %>% subsample_B(nsel_b = B_control$nsel_b)
    }
  }
    
  ## t-nor mod3
  if(model == "t_nor_mod3" & bh == "weibull"){
    
    if(is.null(deriv)){
      mod <- rstan::stan_model(model_code = new_rand_eff_t_nor_mod3_jm_weibull, auto_write = TRUE)
    }else{
      mod <- rstan::stan_model(model_code = new_rand_eff_t_nor_mod3_jm_weibull_deriv, auto_write = TRUE)
    }
    
    data_t_nor_mod3 <- list(ntot = ntot,
                          id = as.array(l_id), 
                          y = as.array(y), 
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
                          wt_quad = wt_quad,
                          d_ind = d_ind,
                          Q_ind = Q_ind)
    
    if(!is.null(deriv)){
      data_t_nor_mod3$x_deriv_T <- x_deriv_T
      data_t_nor_mod3$x_deriv_quad <- x_deriv_quad
      data_t_nor_mod3$d_deriv_T <- d_deriv_T
      data_t_nor_mod3$d_deriv_quad <- d_deriv_quad
      data_t_nor_mod3$p_deriv <- length(deriv_alpha_ind)
      data_t_nor_mod3$q_deriv <- length(deriv_B_ind)
      data_t_nor_mod3$deriv_alpha_ind <- as.array(deriv_alpha_ind)
      data_t_nor_mod3$deriv_B_ind <- as.array(deriv_B_ind)
    }
    
    B_sampled <- list()
    
    for(i in 1:M){
      data_t_nor_mod3$alpha      <- as.array(alpha[i, ])
      data_t_nor_mod3$Sigma      <- Sigma[[i]]
      data_t_nor_mod3$sigma_Z    <- sigma_Z[i, ]
      data_t_nor_mod3$log_lambda <- log_lambda[i, ]
      data_t_nor_mod3$log_nu     <- log_nu[i, ]
      data_t_nor_mod3$omega      <- as.array(omega[i, ])
      
      if(is.null(deriv)){
        data_t_nor_mod3$eta <- eta[i, ]
      }else{
        data_t_nor_mod3$eta1 <- eta1[i, ]
        data_t_nor_mod3$eta2 <- eta2[i, ]
      }
      
      data_t_nor_mod3$phi        <- phi[i, ]

      B_res <- rstan::sampling(mod, 
                               data = data_t_nor_mod3, 
                               iter = B_control$iter, 
                               warmup = B_control$warmup,
                               chains = B_control$chains,
                               cores = B_control$cores,
                               init = B_control$init,
                               control = list(adapt_delta = B_control$adapt_delta, 
                                              max_treedepth = B_control$max_treedepth)
      )
      
      #B_sampled[[i]] <- matrix(rstan::summary(B_res)$summary[1:(ngroup*q), "50%"], ncol = q, byrow = T)
      B_sampled[[i]] <- rstan::extract(B_res)[["B"]] %>% subsample_B(nsel_b = B_control$nsel_b)
    }
  }
  
  ## t-t mod3
  if(model == "t_t_mod3" & bh == "weibull"){
    
    if(is.null(deriv)){
      mod <- rstan::stan_model(model_code = new_rand_eff_t_t_mod3_jm_weibull, auto_write = TRUE)
    }else{
      mod <- rstan::stan_model(model_code = new_rand_eff_t_t_mod3_jm_weibull_deriv, auto_write = TRUE)
    }
    
    data_t_t_mod3 <- list(ntot = ntot,
                          id = as.array(l_id), 
                          y = as.array(y), 
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
                          wt_quad = wt_quad,
                          d_ind = d_ind,
                          Q_ind = Q_ind)
    
    if(!is.null(deriv)){
      data_t_t_mod3$x_deriv_T <- x_deriv_T
      data_t_t_mod3$x_deriv_quad <- x_deriv_quad
      data_t_t_mod3$d_deriv_T <- d_deriv_T
      data_t_t_mod3$d_deriv_quad <- d_deriv_quad
      data_t_t_mod3$p_deriv <- length(deriv_alpha_ind)
      data_t_t_mod3$q_deriv <- length(deriv_B_ind)
      data_t_t_mod3$deriv_alpha_ind <- as.array(deriv_alpha_ind)
      data_t_t_mod3$deriv_B_ind <- as.array(deriv_B_ind)
    }
    
    B_sampled <- list()
    
    for(i in 1:M){
      data_t_t_mod3$alpha      <- as.array(alpha[i, ])
      data_t_t_mod3$Sigma      <- Sigma[[i]]
      data_t_t_mod3$sigma_Z    <- sigma_Z[i, ]
      data_t_t_mod3$log_lambda <- log_lambda[i, ]
      data_t_t_mod3$log_nu     <- log_nu[i, ]
      data_t_t_mod3$omega      <- as.array(omega[i, ])
      
      if(is.null(deriv)){
        data_t_t_mod3$eta <- eta[i, ]
      }else{
        data_t_t_mod3$eta1 <- eta1[i, ]
        data_t_t_mod3$eta2 <- eta2[i, ]
      }
      
      data_t_t_mod3$phi        <- phi[i, ]
      data_t_t_mod3$delta      <- delta[i, ]
      
      B_res <- rstan::sampling(mod, 
                        data = data_t_t_mod3, 
                        iter = B_control$iter, 
                        warmup = B_control$warmup,
                        chains = B_control$chains,
                        cores = B_control$cores,
                        init = B_control$init,
                        control = list(adapt_delta = B_control$adapt_delta, 
                                       max_treedepth = B_control$max_treedepth)
                        )
      
      #B_sampled[[i]] <- matrix(rstan::summary(B_res)$summary[1:(ngroup*q), "50%"], ncol = q, byrow = T)
      B_sampled[[i]] <- rstan::extract(B_res)[["B"]] %>% subsample_B(nsel_b = B_control$nsel_b)
    }
  }
  
  ##
  ## the calculate the survival probabilities by pluggin in the ratio 
  ##
  
  ft <- list()
  for(i in 1:ngroup){
    ft[[i]] <- seq(S[i], (S[i] + forecast$h), length.out = forecast$n)
  }
  
  ft_batch_data_base <- batch_data[!duplicated(l_id), ]
  ft_batch_data_quad <- ft_batch_data_base[rep(1:ngroup, each = Q), ]
  
  ft_probs <- list()
  
  ## number of B samples
  if(B_control$nsel_b == "all"){
    B_length <- B_control$iter - B_control$warmup
  }else if(B_control$nsel_b %in% c("mean", "median")){
    B_length <- 1
  }else{
    B_length <- B_control$nsel_b
  }
  
  probs_length <- B_length * M

  for(i in 1:ngroup){
    
    ft_i <- ft[[i]]
    
    ft_probs_i <- list()
    ft_probs_i[[1]] <- rep(1, probs_length)#rep(1, M)
    
    ft_batch_data_quad_i <- ft_batch_data_quad[((i-1)*Q+1):(i*Q), , drop = FALSE]
    c_quad_i <- c_quad[((i-1)*Q+1):(i*Q), , drop = FALSE]
    
    for(j in 2:forecast$n){
      
      ft_probs_i_k <- c()
      
      for(k in 1:M){
        
        for(kk in 1:B_length){
          
          if(is.null(deriv)){
            prob_upper <- surv_prob_calc(t = ft_i[j],
                                         ft_batch_data_quad_i = ft_batch_data_quad_i,
                                         c_quad_i = c_quad_i,
                                         timeVar = timeVar, 
                                         log_lambda = log_lambda[k, ],
                                         log_nu = log_nu[k, ],
                                         omega = omega[k, ],
                                         eta = eta[k, ],
                                         alpha = alpha[k, ],
                                         B = switch(ifelse(B_length == 1, 1, 2), B_sampled[[k]][i, ], B_sampled[[k]][[i]][kk,]),
                                         #B = switch(B_length == 1, B_sampled[[k]][,i,], B_sampled[[k]][,i,][kk,]),
                                         #B = B_sampled[[k]][i, ],
                                         wt = wt, 
                                         pt = pt,
                                         Q = Q,
                                         bh = "weibull",
                                         deriv = deriv,
                                         fixed_long = object$fixed_long,
                                         random_long = object$random_long)
            prob_lower <- surv_prob_calc(t = ft_i[1], 
                                         ft_batch_data_quad_i = ft_batch_data_quad_i,
                                         c_quad_i = c_quad_i,
                                         timeVar = timeVar, 
                                         log_lambda = log_lambda[k, ],
                                         log_nu = log_nu[k, ],
                                         omega = omega[k, ],
                                         eta = eta[k, ],
                                         alpha = alpha[k, ],
                                         B = switch(ifelse(B_length == 1, 1, 2), B_sampled[[k]][i, ], B_sampled[[k]][[i]][kk,]),
                                         #B = switch(B_length == 1, B_sampled[[k]][,i,], B_sampled[[k]][,i,][kk,]),
                                         #B = B_sampled[[k]][i, ],
                                         wt = wt, 
                                         pt = pt,
                                         Q = Q,
                                         bh = "weibull",
                                         deriv = deriv,
                                         fixed_long = object$fixed_long,
                                         random_long = object$random_long)
          }else{
            prob_upper <- surv_prob_calc(t = ft_i[j],
                                         ft_batch_data_quad_i = ft_batch_data_quad_i,
                                         c_quad_i = c_quad_i,
                                         timeVar = timeVar, 
                                         log_lambda = log_lambda[k, ],
                                         log_nu = log_nu[k, ],
                                         omega = omega[k, ],
                                         eta1 = eta1[k, ],
                                         eta2 = eta2[k, ],
                                         alpha = alpha[k, ],
                                         B = switch(ifelse(B_length == 1, 1, 2), B_sampled[[k]][i, ], B_sampled[[k]][[i]][kk,]),
                                         #B = switch(B_length == 1, B_sampled[[k]][,i,], B_sampled[[k]][,i,][kk,]),
                                         #B = B_sampled[[k]][i, ],
                                         wt = wt, 
                                         pt = pt,
                                         Q = Q,
                                         bh = "weibull",
                                         deriv = deriv,
                                         fixed_long = object$fixed_long,
                                         random_long = object$random_long)
            prob_lower <- surv_prob_calc(t = ft_i[1], 
                                         ft_batch_data_quad_i = ft_batch_data_quad_i,
                                         c_quad_i = c_quad_i,
                                         timeVar = timeVar, 
                                         log_lambda = log_lambda[k, ],
                                         log_nu = log_nu[k, ],
                                         omega = omega[k, ],
                                         eta1 = eta1[k, ],
                                         eta2 = eta2[k, ],
                                         alpha = alpha[k, ],
                                         B = switch(ifelse(B_length == 1, 1, 2), B_sampled[[k]][i, ], B_sampled[[k]][[i]][kk,]),
                                         #B = switch(B_length == 1, B_sampled[[k]][,i,], B_sampled[[k]][,i,][kk,]),
                                         #B = B_sampled[[k]][i, ],
                                         wt = wt, 
                                         pt = pt,
                                         Q = Q,
                                         bh = "weibull",
                                         deriv = deriv,
                                         fixed_long = object$fixed_long,
                                         random_long = object$random_long)
          }
          
          #ft_probs_i_k <- c(ft_probs_i_k, prob_upper/prob_lower)
          ft_probs_i_k <- c(ft_probs_i_k, exp(prob_upper - prob_lower))
        }
          
        }
        

      ft_probs_i[[j]] <- ft_probs_i_k
    }
    
    ft_probs[[i]] <- ft_probs_i 
    names(ft_probs)[[i]] <- as.character(s_id_orig[i])
    
  }
  
  ft_table <- list()
  
  for(i in 1:ngroup){
    ft_i <- ft[[i]]
    ft_table_i <- data.frame(id = rep(s_id_orig[i], forecast$n),
                        #stime = rep(S[i], forecast$n),
                        #event = rep(E[i], forecast$n),
                        time = ft_i, 
                        do.call(rbind, lapply(ft_probs[[i]], prob_summary, probs = probs)))
    names(ft_table_i)[3:ncol(ft_table_i)] <- c("mean", paste0((probs*100), "%"))
    ft_table[[i]] <- ft_table_i
  }
  
  ft_table <- do.call(rbind, ft_table)
  
  if(return_bsamples){
    out <- list(ft_probs = ft_probs, ft_table = ft_table, B_sampled = B_sampled)
  }else{
    out <- list(ft_probs = ft_probs, ft_table = ft_table)
  }
  return(out)
  
}