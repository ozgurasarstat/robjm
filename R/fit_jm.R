
fit_jm <- function(fixed_long, 
                   random_long, 
                   fixed_surv,
                   data_long, 
                   data_surv, 
                   id_long, 
                   id_surv,
                   model, 
                   bh, #baseline hazard 
                   spline_tv, 
                   Q = 15, 
                   ...){

  ## be sure that distribution specifications are correct
  if(!(model %in% c("nor_nor", "t_t_mod1", "t_t_mod2", "t_t_mod3", "t_t_tv"))){
    stop("Model should be one of the followings: nor_nor, t_t_mod1, t_t_mod2, t_t_mod3, t_t_tv")
  }
  
  ## be sure that id is: 1, 2, 3, ...
  data_long[, id_long] <- rep(1:length(unique(data_long[, id_long])), as.numeric(table(data_long[, id_long])))
  
  ## x and y matrices
  x <- model.matrix(fixed_long, data_long)
  y <- model.frame(fixed_long, data_long)[, 1]
  
  ## create blok-diagonal random effects design matrix
  id_dmat <- data.frame(data_long[, id_long], model.matrix(random_long, data_long))
  id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
  d <- do.call(magic::adiag, id_dmat_list)
  
  ntot <- nrow(x) 
  p <- ncol(x)
  q <- ncol(id_dmat) - 1
  
  ngroup <- length(unique(id_long))
  
  x <- cbind(1, long_data$obstime)
  
  gl_quad <- statmod::gauss.quad(Q)
  wt <- gl_quad$weights
  pt <- gl_quad$nodes
  
  ntot_quad <- ngroup * Q
  
  mf_surv <- model.frame(fixed_surv, data)
  S <- mf_surv[, 1]
  E <- mf_surv[, 2]
  
  t_quad <- 0.5 * rep(S, each = Q) * (1 + rep(pt, ngroup))
  
  ncol_e <- 2
  median_S <- median(S)
  e <- cbind(ifelse(S < median_S, 1, 0), ifelse(S >= median_S, 1, 0))#bs(S, df = 3)
  #attributes(e) <- NULL
  #e <- matrix(e, ncol = 3)
  e_quad <- cbind(ifelse(t_quad < median_S, 1, 0), ifelse(t_quad >= median_S, 1, 0))#bs(t_quad, df = 3)
  #attributes(e_quad) <- NULL
  #e_quad <- matrix(e_quad, ncol = 3)
  
  c <- model.matrix(cbind(a, b) ~ c, data)[, -1, drop = FALSE]
  ncol_c <- ncol(c)
  c_quad <- apply(c, 2, function(i) rep(i, each = Q))
  
  x_T <- cbind(1, S)
  x_quad <- cbind(1, t_quad)
  
  id_dmat_T <- data.frame(surv_data[, id_surv], cbind(1, S))
  id_dmat_list_T <- lapply(split(id_dmat_T[, -1], id_dmat_T[, 1]), as.matrix)
  d_T <- do.call(magic::adiag, id_dmat_list_T)
  
  id_dmat_quad <- data.frame(rep(surv_data[, id_surv], each = Q), cbind(1, t_quad))
  id_dmat_list_quad <- lapply(split(id_dmat_quad[, -1], id_dmat_quad[, 1]), as.matrix)
  d_quad <- do.call(magic::adiag, id_dmat_list_quad)
  
  wt_quad <- rep(wt, ngroup)
  priors_long <- c(5, 2, 5, 5, 5)
  priors_surv <- c(5, 5, 5)
  
  if(model == "t_t_tv"){
    a <- splines::ns(long_data[, "obstime"], df = 3)
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
    
    res <- stan(model_code = t_t_tv_jm,
                data = data_nor_nor_jm,
                ...
    )
  }

  return(res)  
  
}

