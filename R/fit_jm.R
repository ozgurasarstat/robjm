
fit_jm <- function(fixed_long, 
                   random_long, 
                   fixed_surv,
                   data_long, 
                   data_surv, 
                   id_long, 
                   id_surv,
                   model, 
                   bh, #baseline hazard 
                   spline_tv, #spline for tv dof - same in fit_ld
                   Q = 15, 
                   ...){

  ## be sure that distribution specifications are correct
  if(!(model %in% c("nor_nor", "t_t_mod1", "t_t_mod2", "t_t_mod3", "t_t_tv"))){
    stop("Model should be one of the followings: nor_nor, t_t_mod1, t_t_mod2, t_t_mod3, t_t_tv")
  }
  
  ## be sure that id is: 1, 2, 3, ...
  data_long[, id_long] <- rep(1:length(unique(data_long[, id_long])), as.numeric(table(data_long[, id_long])))
  id <- data_long[, id_long]
  
  ## x and y matrices
  x <- model.matrix(fixed_long, data_long)
  y <- model.frame(fixed_long, data_long)[, 1]
  
  ## create blok-diagonal random effects design matrix
  id_dmat <- data.frame(data_long[, id_long], model.matrix(random_long, data_long))
  id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
  d <- do.call(magic::adiag, id_dmat_list)
  
  ## total number of observations in the longitudinal data
  ## number of covariates in the x and d matrices
  ntot <- nrow(x) 
  p <- ncol(x)
  q <- ncol(id_dmat) - 1
  
  ## number of subjects
  ngroup <- length(unique(data_long[, id_long]))
  
  ## Gauss - Legendre weights and abscissas
  gl_quad <- statmod::gauss.quad(Q)
  wt <- gl_quad$weights
  pt <- gl_quad$nodes
  
  ## total number of observations for quadrature approx.
  ntot_quad <- ngroup * Q
  
  ## extract survival times and event indicator
  mf_surv <- model.frame(fixed_surv, data_surv)
  S <- mf_surv[, 1]
  E <- mf_surv[, 2]
  
  ## calculate times for hazard function for quadrature approx
  t_quad <- 0.5 * rep(S, each = Q) * (1 + rep(pt, ngroup))
  
  ## baseline hazard - piecewise constant with one know a.t.m. 
  ncol_e <- 2
  median_S <- median(S)
  e <- cbind(ifelse(S < median_S, 1, 0), ifelse(S >= median_S, 1, 0))#bs(S, df = 3)
  #attributes(e) <- NULL
  #e <- matrix(e, ncol = 3)
  e_quad <- cbind(ifelse(t_quad < median_S, 1, 0), ifelse(t_quad >= median_S, 1, 0))#bs(t_quad, df = 3)
  #attributes(e_quad) <- NULL
  #e_quad <- matrix(e_quad, ncol = 3)
  
  ## fixed effects for survival sub-model
  c <- model.matrix(fixed_surv, data_surv)[, -1, drop = FALSE]
  ncol_c <- ncol(c)
  c_quad <- apply(c, 2, function(i) rep(i, each = Q))
  
  ## x matrix for log survival density
  x_T <- cbind(1, S) 
  x_quad <- cbind(1, t_quad)
  
  ## d matrix for log survival density
  id_dmat_T <- data.frame(data_surv[, id_surv], cbind(1, S))
  id_dmat_list_T <- lapply(split(id_dmat_T[, -1], id_dmat_T[, 1]), as.matrix)
  d_T <- do.call(magic::adiag, id_dmat_list_T)
  
  id_dmat_quad <- data.frame(rep(data_surv[, id_surv], each = Q), cbind(1, t_quad))
  id_dmat_list_quad <- lapply(split(id_dmat_quad[, -1], id_dmat_quad[, 1]), as.matrix)
  d_quad <- do.call(magic::adiag, id_dmat_list_quad)
  
  ## extend the weights for quadrature approx.
  wt_quad <- rep(wt, ngroup)
  
  ## prior hyperparameters
  priors_long <- c(5, 2, 5, 5, 5)
  priors_surv <- c(5, 5, 5)
  
  if(model == "nor_nor"){
    data_nor_nor <- list(ntot = ntot,
                        id = id, 
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
                        ncol_e = ncol_e, 
                        e = e, 
                        e_quad = e_quad,
                        ncol_c = ncol_c, 
                        c = c, 
                        c_quad = c_quad,
                        x_T = x_T,
                        x_quad = x_quad,
                        d_T = d_T,
                        d_quad = d_quad,
                        wt_quad = wt_quad
                        )
    
    res <- stan(model_code = nor_nor_jm, data = data_nor_nor_jm, ...)
    
  }
  
  if(model %in% c("t_t_mod1", "t_t_mod2", "t_t_mod3")){
    data_t_t <- list(ntot = ntot,
                         id = id, 
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
                         ncol_e = ncol_e, 
                         e = e, 
                         e_quad = e_quad,
                         ncol_c = ncol_c, 
                         c = c, 
                         c_quad = c_quad,
                         x_T = x_T,
                         x_quad = x_quad,
                         d_T = d_T,
                         d_quad = d_quad,
                         wt_quad = wt_quad
                     )

    if(model == "t_t_mod1"){
      res <- stan(model_code = t_t_jm_mod1, data = data_t_t, ...)
    } 
    if(model == "t_t_mod2"){
      res <- stan(model_code =t_t_jm_mod2, data = data_t_t, ...)
    }
    if(model == "t_t_mod3"){
      res <- stan(model_code = t_t_jm_mod3, data = data_t_t, ...)
    }
    
 }

  
  if(model == "t_t_tv"){
    a <- splines::ns(data_long[, spline[[1]]], df = spline[[2]])
    ncol_a <- ncol(a)
    attributes(a) <- NULL
    a <- matrix(a, ncol = ncol_a)
    s <- ncol_a
    
    data_t_t_tv <- list(ntot = ntot,
                        id = id, 
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
                        ncol_e = ncol_e, 
                        e = e, 
                        e_quad = e_quad,
                        ncol_c = ncol_c, 
                        c = c, 
                        c_quad = c_quad,
                        x_T = x_T,
                        x_quad = x_quad,
                        d_T = d_T,
                        d_quad = d_quad,
                        wt_quad = wt_quad,
                        s = s,
                        a = a
    )
    
    res <- stan(model_code = t_t_tv_jm, data = data_nor_nor_jm, ...)
    
  }

  return(res)  
  
}