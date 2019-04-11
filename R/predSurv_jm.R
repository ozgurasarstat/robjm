#' @title Predict survival probabilities
#' 
#' @description A function to predict survival probabilities
#' 
#' @param object A fitted object from \code{fit_jm}
#' @param newdata A data frame for new patients in the longitudinal format
#' @param last_time A character for indicating the survival time, surv_time, column name, or landmark
#' @param forecast A numeric vector of three elements, 
#'       first element is landmark time, second is horizon, third is increment
#' 

predSurv_jm <- function(object, 
                        newdata, 
                        last_time = "surv_time", 
                        lm_time = NULL,
                        forecast = list(h = 5, n = 5), 
                        B_control = list(iter = 400, 
                                         warmup = 200, 
                                         chains = 1, 
                                         cores = 1,
                                         init = "random",
                                         adapt_delta = 0.8, 
                                         max_treedepth = 10),
                        batch_control = list(size = 100, 
                                             cores = 1),
                        probs = c(0.025, 0.5, 0.975),
                        return_bsamples = FALSE,
                        ...){

  ## be sure that B_control has 5 elements
  if(length(B_control) < 7){
    B_control_f <- list(iter = 500, warmup = 250, chains = 1, cores = 1,
                        init = "random",
                        adapt_delta = 0.8, max_treedepth = 10)
    for(i in 1:7){
      if(!(names(B_control_f)[i] %in% names(B_control))){
        B_control[names(B_control_f)[i]] <- B_control_f[i]
      }
    }
  }
  
  ## be sure that batch_control has 2 elements
  if(length(batch_control) < 2){
    batch_control_f <- list(size = 100, cores = 1)
    for(i in 1:2){
      if(!(names(batch_control_f)[i] %in% names(batch_control))){
        batch_control[names(batch_control_f)[i]] <- batch_control_f[i]
      }
    }
  }
  
  ## Gauss-Legendre stuff
  Q <- object$Q
  gl_quad <- statmod::gauss.quad(Q)
  wt <- gl_quad$weights
  pt <- gl_quad$nodes
  
  ## create batches
  
  newdata <- newdata[order(newdata[, object$id_long], newdata[, object$timeVar]), ]
  
  idlist <- newdata[, object$id_long] 
  nobs <- idlist %>% table %>% as.numeric
  nsubj <- idlist %>% unique %>% length
  index <- rep(1:nsubj, nobs)
  chunk <- ceiling(index/batch_control$size)
  iterations <- max(chunk)
  
  newdata_batch <- split(newdata, chunk)
  
  chunk_sizes <- lapply(newdata_batch, function(x) unique_length(x[, object$id_long])) %>% unlist

  if(batch_control$cores == 1){

    pred_out <- list()
    
    for(i in 1:iterations){
      
      pred_out[[i]] <- predSurv_jm_supp(batch_data = newdata_batch[[i]], 
                                        object = object,
                                        forecast = forecast, 
                                        B_control = B_control, 
                                        Q = Q,
                                        wt = wt,
                                        pt = pt,
                                        last_time = last_time,
                                        lm_time = lm_time,
                                        probs = probs,
                                        return_bsamples = return_bsamples)
      
    }#for(i in 1:iterations){
    
  }else{
    
    cl <- makeCluster(batch_control$cores)
    
    registerDoSNOW(cl)
    
    pb <- txtProgressBar(max = iterations, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    
    clusterExport(cl, list = c('newdata_batch', 'forecast', 'B_control',
                               'object', 'Q', "wt", "pt"), envir = environment())
    
    pred_out <- foreach(i = 1:iterations, .options.snow = opts, .packages = c("rstan", "magic")) %dopar%
    {
      
      pred_i <- predSurv_jm_supp(batch_data = newdata_batch[[i]], 
                                 object = object,
                                 forecast = forecast, 
                                 B_control = B_control, 
                                 Q = Q,
                                 wt = wt,
                                 pt = pt,
                                 last_time = last_time,
                                 lm_time = lm_time,
                                 probs = probs,
                                 return_bsamples = return_bsamples)
      return(pred_i)
      
    }
    close(pb)  
    stopCluster(cl)
      
  }

  ## combine results  
  out <- combine_pred(x = pred_out, 
                      iterations = iterations, 
                      nsubj = nsubj,
                      chunk_sizes = chunk_sizes)
  
  if(return_bsamples){
    out$bsamples <- combine_bsamples(x  = pred_out, 
                                     iterations = iterations, 
                                     nsubj = nsubj,
                                     chunk_sizes = chunk_sizes)
    
    names(out$bsamples) <- names(out$samples)    
  }

  return(out)
  
}
