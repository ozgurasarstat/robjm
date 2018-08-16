
 fit_ld <- function(fixed, 
                    random, 
                    data, 
                    id, 
                    model, 
                    spline, 
                    priors = list(), 
                    ...){

   ## fixed: two-sided formula for fixed effects
   ## random: one-sided formula for random effects
   ## data: data frame
   ## id: name of the id column
   ## model: model identifier, "nor_nor", "t_t_mod1", "t_t_mod2", "t_t_mod3", "t_t_tv"
   ## spline: a list: name of the time variable, and number of knots plus 1
   ## priors: prior hyperparameters, order = theta, omega, sigma_Bstar, sigma_B, beta (for tv)
   ## ... to be passed to stan() function

   ## be sure that distribution specifications are correct
   if(!(model %in% c("nor_nor", "t_t_mod1", "t_t_mod2", "t_t_mod3", "t_t_tv"))){
     stop("Model should be one of the followings: nor_nor, t_t_mod1, t_t_mod2, t_t_mod3, t_t_tv")
   }

   ## re-organise priors
   if(model == "t_t_tv" & length(priors) != 5){
     priors_full <- list(theta = 5, 
                         Omega = 2, 
                         sigma_B = 5, 
                         sigma_Z = 5,
                         beta = 4.6)
     for(i in 1:5){
       if(!(names(priors_full)[i] %in% names(priors))){
         priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
       }
     }
   }
   if(model != "t_t_tv" % length(priors) != 4){
     priors_full <- list(theta = 5, 
                         Omega = 2, 
                         sigma_B = 5, 
                         sigma_Z = 5)
     for(i in 1:4){
       if(!(names(priors_full)[i] %in% names(priors))){
         priors[names(priors_full)[i]] <- priors_full[names(priors_full)[i]]
       }
     } 
   }
   
   ## be sure that id is: 1, 2, 3, ...
   data[, id] <- rep(1:length(unique(data[, id])), as.numeric(table(data[, id])))

   ## x and y matrices
   x <- model.matrix(fixed, data)
   y <- model.frame(fixed, data)[, 1]

   ## create blok-diagonal random effects design matrix
   id_dmat <- data.frame(data[, id], model.matrix(random, data))
   id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
   d <- do.call(magic::adiag, id_dmat_list)

   ## Fit the normal - normal model
   if(model == "nor_nor"){
     priors
     dat_nor_nor <- list(ntot = nrow(data),
                         id = data[, id],
                         y = y,
                         x = x,
                         d = d,
                         p = ncol(x),
                         q = ncol(id_dmat) - 1,
                         ngroup = length(unique(data[, id])),
                         priors = unlist(priors)
                         )
     res <- stan(model_code = nor_nor_ld, data = dat_nor_nor, ...)

   }


   ## t-t models - time invariant d.o.f.

   # prepare data-set first, common to t_t_mod1, t_t_mod2, t_t_mod3
   if(model %in% c("t_t_mod1", "t_t_mod2", "t_t_mod3")){
     dat_t_t <- list(ntot = nrow(data),
                     id = data[, id],
                     y = y,
                     x = x,
                     d = d,
                     p = ncol(x),
                     q = ncol(id_dmat) - 1,
                     ngroup = length(unique(data[, id])),
                     priors = unlist(priors)
                     )
   }

   if(model == "t_t_mod1"){
     res <- stan(model_code = t_t_ld_mod1, data = dat_t_t, ...)
   }


   if(model == "t_t_mod2"){
     res <- stan(model_code = t_t_ld_mod2, data = dat_t_t, ...)
   }

   if(model == "t_t_mod3"){
     res <- stan(model_code = t_t_ld_mod3, data = dat_t_t, ...)
   }

   ## time-varying d.o.f for Z

   if(model == "t_t_tv"){
     a <- splines::ns(data[, spline[[1]]], df = spline[[2]])
     ncol_a <- ncol(a)
     attributes(a) <- NULL
     a <- matrix(a, ncol = ncol_a)

     dat_tv <- list(ntot = nrow(data),
                    id = data[, id],
                    y = y,
                    x = x,
                    d = d,
                    p = ncol(x),
                    q = ncol(id_dmat) - 1,
                    ngroup = length(unique(data[, id])),
                    s = ncol(a),
                    a = a,
                    priors = unlist(priors)
                    )

     res <- stan(model_code = t_t_tv_ld, data = dat_tv, ...)
   }

 return(res)

 }
