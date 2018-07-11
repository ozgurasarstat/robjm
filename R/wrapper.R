 
 wrapper <- function(fixed, random, data, id, model, spline, ...){
   
   ## fixed: two-sided formula for fixed effects
   ## random: one-sided formula for random effects
   ## data: data frame
   ## id: name of the id column
   ## model: model identifier: "nor_nor", "t_t_mod1", "t_t_mod2", "t_t_mod3", "t_t_tv"
   ## spline: a list: name of the time variable, and number of knots 
   
   ## set the working directory to the folder where the codes are located
   wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
   setwd(wd)
   
   ## packages to be used
   library(rstan)
   library(magic)
   library(splines)

   ## be sure that id is: 1, 2, 3, ...
   data[, id] <- rep(1:length(unique(data[, id])), as.numeric(table(data[, id])))
      
   ## x matrix
   x <- model.matrix(fixed, data)
   y <- as.matrix(model.frame(fixed, data)[, 1])
        
   ## create blok-diagonal random effects design matrix
   id_dmat <- data.frame(data[, id], model.matrix(random, data))
   id_dmat_list <- lapply(split(id_dmat[, -1], id_dmat[, 1]), as.matrix)
   d <- do.call(adiag, id_dmat_list)
   

   ## Fit the normal - normal model
   id(model = "nor_nor"){
     dat_nor_nor <- list(ntot = nrow(data),
                         id = data[, "id"],
                         y = y, 
                         x = x,
                         d = d,
                         p = ncol(x),
                         q = ncol(id_dmat) - 1,
                         ngroup = length(unique(data[, "id"])),
                         priors = c(5, 2, 5, 5)
                         ) 
     
     source("normal_normal.R")
     res_nor_nor <- stan(model_code = normal_normal, 
                         data = dat_nor_nor,
                         ...
                         )
     res_nor_nor     
   }

   
   ## t-t models
   
   dat_t_t <- list(ntot = nrow(Orthodont),
                   id = Orthodont$id,
                   y = Orthodont$distance, 
                   x = x,
                   d = d,
                   p = ncol(x),
                   q = 2,
                   ngroup = length(unique(Orthodont$id)),
                   priors = c(5, 2, 5, 5)
   ) 
   
   source("t_t_mod1.R")
   res_t_t_mod1 <- stan(model_code = t_t_mod1, 
                        data = dat_t_t,
                        chains = 2, 
                        iter = 2000, 
                        warmup = 500,
                        control = list(adapt_delta = 0.99)
   )
   print(res_t_t_mod1, digits = 4, pars = c("alpha", "Sigma", "sigmasq", "phi"))
   
   source("t_t_mod2.R")
   res_t_t_mod2 <- stan(model_code = t_t_mod2, 
                        data = dat_t_t,
                        chains = 2, 
                        iter = 2000, 
                        warmup = 500,
                        control = list(adapt_delta = 0.99)
   )
   print(res_t_t_mod2, digits = 4, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta"))
   
   source("t_t_mod3.R")
   res_t_t_mod3 <- stan(model_code = t_t_mod3, 
                        data = dat_t_t,
                        chains = 2, 
                        iter = 2000, 
                        warmup = 500,
                        control = list(adapt_delta = 0.99)
   )

   ### TIME-VARYING
   
   a <- ns(Orthodont$age2, df = 1)
   attributes(a) <- NULL
   a <- matrix(a, ncol = 1)
   
   dat_tv <- list(ntot = nrow(Orthodont),
                  id = Orthodont$id,
                  y = Orthodont$distance, 
                  x = x,
                  d = d,
                  p = ncol(x),
                  q = 2,
                  ngroup = length(unique(Orthodont$id)),
                  s = ncol(a),
                  a = a,
                  priors = c(5, 2, 5, 5, 4.6)
   )
   
   source("t_t_tv.R")
   res_t_t_tv <- stan(model_code = t_t_tv, 
                      data = dat_tv,
                      ...)
   res_t_t_tv

 }