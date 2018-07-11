
## set the working directory to the folder
## where the codes are located
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

## packages to be used
library(rstan)
library(nlme)
library(magic)
library(splines)

## load the Orthodont data-set
data(Orthodont) 
head(Orthodont) 

## create id column taking values: 1, 2, 3, ...
## and center age
Orthodont$id <- rep(1:27, each = 4)
Orthodont$age2 <- Orthodont$age - 8

## creat blok-diagonal random effects design matrix
dmat <- cbind(1, Orthodont$age2)
d <- dmat[1:4, ]

for(i in 2:27) d <- adiag(d, dmat[((i-1)*4+1): (i*4),])  

formula_fixed <- distance ~ age2
x <- model.matrix(formula_fixed, data = Orthodont)

## Fit the normal - normal model

dat_nor_nor <- list(ntot = nrow(Orthodont),
            id = Orthodont$id,
            y = Orthodont$distance, 
            x = x,
            d = d,
            p = ncol(x),
            q = 2,
            ngroup = length(unique(Orthodont$id)),
            priors = c(5, 2, 5, 5)
            ) 

source("normal_normal.R")
res_nor_nor <- stan(model_code = normal_normal, 
            data = dat_nor_nor,
            chains = 2, 
            iter = 2000, 
            warmup = 500,
            control = list(adapt_delta = 0.99)
            )
print(res_nor_nor, digits = 4, pars = c("alpha", "Sigma", "sigmasq"))
traceplot(res_nor_nor, pars = "Sigma")

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
print(res_t_t_mod3, digits = 4, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta"))

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

res_t_t_tv <- stan(model_code = t_t_tv, 
                data = dat_tv,
                chains = 2, 
                iter = 2000, 
                warmup = 1000,
                control = list(adapt_delta = 0.99)
                )

print(res_t_t_tv, digits = 4, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta0", "beta"))
traceplot(res_t_t_tv, pars = "beta")

lmefit <- lme(distance ~ age2, data = Orthodont, random = ~ age2|id)
summary(lmefit)
