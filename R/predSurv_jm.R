#' @title Predict survival probabilities
#' 
#' @description A function to predict survival probabilities
#' 
#' @param object A fitted object from \code{fit_jm}
#' @param newdata A data frame for new patients in the longitudinal format
#' @param forecast A numeric vector of three elements, 
#'       first element is landmark time, second is horizon, third is increment
#' 

predSurv_jm <- function(object, 
                        newdata, 
                        forecast = list(h = 5, n = 5), 
                        B_control = list(iter = 500, 
                                         warmup = 250, 
                                         chains = 1, 
                                         adapt_delta = 0.8, 
                                         max_treedepth = 10),
                        batch_control = list(size = 100, 
                                             cores = 1),
                        ...){

  ## be sure that B_control has 3 elements
  if(length(B_control) < 5){
    B_control_f <- list(iter = 500, warmup = 250, chains = 1, adapt_delta = 0.8, max_treedepth = 10)
    for(i in 1:5){
      if(!(names(B_control_f)[i] %in% names(B_control))){
        B_control[names(B_control_f)[i]] <- B_control_f[i]
      }
    }
  }
  
  ## extract the chains
  chains <- list()
  chains$alpha      <- rstan::extract(object$res)$alpha
  chains$Sigma_long <- rstan::extract(object$res)$Sigma
  M                 <- nrow(chains$alpha)
  chains$Sigma      <- lapply(1:M, function(i) Sigma_long[i, ,])
  chains$sigma_Z    <- matrix(rstan::extract(object$res)$sigma_Z)
  chains$log_lambda <- matrix(rstan::extract(object$res)$log_lambda)
  chains$log_nu     <- matrix(rstan::extract(object$res)$log_nu)
  chains$omega      <- rstan::extract(object$res)$omega
  chains$eta        <- matrix(rstan::extract(object$res)$eta)
  
  if(model == "t_t_mod3"){
    chains$phi   <- matrix(rstan::extract(object$res)$phi)
    chains$delta <- matrix(rstan::extract(object$res)$delta)
  }
  
  ## Gauss-Legendre stuff
  gl_quad <- statmod::gauss.quad(object$Q)
  wt <- gl_quad$weights
  pt <- gl_quad$nodes
  
  ## create batches
  
  newdata <- newdata[order(newdata[, object$id_long], newdata[, object$timeVar]), ]
  
  idlist <- newdata[, object$id_long] 
  nobs <- idlist %>% table %>% as.numeric
  nsubj <- idlist %>% unique %>% nrow
  index <- rep(1:nsubj, nobs)
  chunk <- ceiling(index/batch_control$size)
  
  newdata_batch <- split(newdata, chunk)

  
}
