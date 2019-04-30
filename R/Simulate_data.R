
simulate_data <- function(nsubj = 100, 
                          av_n_i = 5, 
                          eta = c(0.3, 0),
                          phi = 3, 
                          delta = 5,
                          delta0 = 5,
                          beta = c(0, 0, 0),
                          knots_tv = 2,
                          lambda = 0.04,
                          nu = 1.2,
                          model = "tv",
                          omega = c(0.5),
                          controls_time = list(t_min = 0, t_max = 5, incr = 0.01),
                          alpha = c(1, 0.6, 0.4, 0.2),
                          Sigma = matrix(c(0.6, 0.25, 0.25, 0.3), ncol = 2),
                          sigmasq = 0.25,
                          returns = c("repeat_data", "base_data")){

  # convert vectors to matrices
  eta <- matrix(eta, ncol = 1)
  beta <- matrix(beta, ncol = 1)
  omega <- matrix(omega, ncol = 1)
  alpha <- matrix(alpha, ncol = 1)
  
  # discretise time at fine intervals
  t_min <- controls_time$t_min
  t_max <- controls_time$t_max
  incr  <- controls_time$incr
  t     <- seq(t_min, t_max, incr)
  m     <- length(t)
  
  # total number of observations before censoring and selection
  ntotal <- m * nsubj
  
  # longitudinal fixed-effects covariate matrix
  
  x1 <- rbinom(nsubj, 1, 0.5)
  x1_ext <- rep(x1, each = m)
  t_ext <- rep(t, nsubj)
  
  x <- cbind(1, x1_ext, t_ext, x1_ext * t_ext) 
  
  # random effects covariate matrix
  d <- Matrix::bdiag(lapply(1:nsubj, function(i) cbind(1, t)))
  
  # create Bstar
  Bstar_mat <- mvtnorm::rmvnorm(nsubj, mean = rep(0, ncol(Sigma)), sigma = Sigma)
  Bstar_vec <- as.numeric(t(Bstar_mat))
  
  # create gamma distributed Vi  
  if(model == "nor_nor"){
    V_indv <- rep(1, nsubj)
  }else{
    V_indv <- rgamma(nsubj, shape = phi/2, rate = phi/2)              
  }
  
  # create B
  B_mat <- Bstar_mat
  for(i in 1:ncol(B_mat)){
    B_mat[, i] <- B_mat[, i]/sqrt(V_indv)
  }
  B_vec  <- as.numeric(t(B_mat))
  
  # measurement error
  Z_star <- rnorm(ntotal, 0, sqrt(sigmasq))
  
  # create W
  if(model == "nor_nor"){
    W <- rep(1, ntotal)
  }else if(model == "mod1"){
    W <- rep(V_indv, each = m)
  }else if(model == "mod2"){
    W_indv <- rgamma(nsubj, shape = delta/2, rate = delta/2)              
    W <- rep(W_indv, each = m)  
  }else if(model == "mod3"){
    W <- rgamma(ntotal, shape = delta/2, rate = delta/2)
  }else if(model == "tv"){
    t_ext <- rep(t, nsubj)                                       
    a_spline  <- ns(t_ext, df = (knots_tv + 1))                              
    delta_ext <- delta0 * exp(a_spline %*% beta)
    W <- rep(NA, ntotal)
    for (i in 1 : ntotal){
      W[i] <- rgamma(1, shape = delta_ext[i]/2, rate = delta_ext[i]/2)
    }    
  }
  
  Z <- Z_star/sqrt(W)
  
  # longitudinal measurements w/o error
  Y_star <- as.matrix(x %*% alpha + d %*% B_vec)  

  # first derivative of Y_star
  Y_star_deriv <- rep(alpha[2, 1], ntotal) + rep(B_mat[, 2], each = m)
    
  # longitudinal measurements with error
  Y <- Y_star + Z                                  
  
  # data-set
  data <- data.frame(
    id   = rep(1:nsubj, each = m),
    time = t_ext,
    Y = Y,
    Y_star = Y_star,
    Y_star_deriv = Y_star_deriv,
    B1 = rep(B_mat[, 1], each = m),
    B2 = rep(B_mat[, 2], each = m)
  )

  # survival model covariates
  c <- matrix(x1_ext, ncol = 1)
  data$c <- x1_ext

  # hazard and survival probabilities
  data$hazard <- as.numeric(with(data, lambda * nu * time^(nu - 1) * exp(c %*% omega + cbind(Y_star, Y_star_deriv) %*% eta)))
  a_vec <- c(1, rep(incr, (m - 1)))
  data$surv_prob <- unlist(with(data, tapply(hazard, id, function(x) exp(- unlist(lapply(1:m, function(i) sum(a_vec[1:i] * x[1:i])))))))
  
  # select observations randomly by being sure that everyone has data at baseline, 0: non-selection, 1:selection
  data$sel <- unlist(lapply(1:nsubj, function(i) c(1, rbinom((m - 1), 1, (av_n_i - 1)/(m - 1)))))
  
  # uniform random variables for survival 
  data$unif_rand_survival  <- rep(runif(nsubj), each = m)
  
  # censor at the event
  data_censored_event <- data[data$surv_prob > data$unif_rand_survival, ]
  
  # create stime and event indicator
  stime_event <- 
    do.call(rbind, with(data_censored_event, tapply(time, id, function(x) stime_event_fun(x, t_max = t_max, a = incr))))
  data_censored_event$stime <- stime_event[, 1]
  data_censored_event$event <- stime_event[, 2]
  
  # further select visits informatively
  data_censored_event_sel <- dplyr::filter(data_censored_event, sel == 1)
  
  # final data-set: repeat + base 
  repeat_data <- data_censored_event_sel
  base_data   <- dplyr::filter(repeat_data, time == 0)
  
  # return some data
  mget(returns)
  
}
