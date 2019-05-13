
#library(devtools)
#install_github("ozgurasarstat/robustjoint")
library("robjm")

data(Orthodont)
Orthodont$age2 <- Orthodont$age - 8

## fit the normal model
fit_nor_nor <- fit_ld(fixed = distance ~ age2, 
                          random = ~ age2, 
                          data = Orthodont,
                          id = "Subject",
                          model = "nor_nor",
                          chains = 2,
                          cores = 2,
                          iter = 2000,
                          warmup = 1000,
                          control = list(adapt_delta = 0.999)
                          )
print(fit_nor_nor, pars = c("alpha", "Sigma", "sigmasq"))
traceplot(fit_nor_nor, pars = c("alpha", "Sigma", "sigmasq"))

## fit the t model 1
fit_t_t_mod1 <- fit_ld(fixed = distance ~ age2, 
                          random = ~ age2, 
                          data = Orthodont,
                          id = "Subject",
                          model = "t_t_mod1",
                          chains = 2,
                          cores = 2,
                          iter = 2000,
                          warmup = 1000,
                          control = list(adapt_delta = 0.999)
                       )
print(fit_t_t_mod1, pars = c("alpha", "Sigma", "sigmasq", "phi"))
#traceplot(fit_t_t_mod1, pars = c("alpha", "sigma_Bstar", "Omega", "sigmasq", "phi"))
traceplot(fit_t_t_mod1, pars = c("alpha", "Sigma", "sigmasq", "phi"))

## fit the t model 2
fit_t_t_mod2 <- fit_ld(fixed = distance ~ age2, 
                           random = ~ age2, 
                           data = Orthodont,
                           id = "Subject",
                           model = "t_t_mod2",
                           chains = 2,
                           cores = 2,
                           iter = 2000,
                           warmup = 1000,
                           control = list(adapt_delta = 0.999)
                           )
print(fit_t_t_mod2, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta"))
traceplot(fit_t_t_mod2, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta"))

## fit t model 3
fit_t_t_mod3 <- fit_ld(fixed = distance ~ age2, 
                           random = ~ age2, 
                           data = Orthodont,
                           id = "Subject",
                           model = "t_t_mod3",
                           chains = 2,
                           cores = 2,
                           iter = 10000,
                           warmup = 5000,
                           control = list(adapt_delta = 0.999, 
                                          max_treedepth = 15)
                           )
print(fit_t_t_mod3, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta"))
traceplot(fit_t_t_mod3, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta"))

## fit nor t model 3
fit_nor_t_mod3 <- fit_ld(fixed = distance ~ age2, 
                       random = ~ age2, 
                       data = Orthodont,
                       id = "Subject",
                       model = "nor_t_mod3",
                       chains = 2,
                       cores = 2,
                       iter = 2000,
                       warmup = 1000,
                       control = list(adapt_delta = 0.999, 
                                      max_treedepth = 15)
                       )
print(fit_nor_t_mod3, pars = c("alpha", "Sigma", "sigmasq", "delta"))
traceplot(fit_nor_t_mod3, pars = c("alpha", "Sigma", "sigmasq", "delta"))

## fit t with tv dof
fit_t_t_tv <- fit_ld(fixed = distance ~ age2, 
                           random = ~ age2, 
                           data = Orthodont,
                           id = "Subject",
                           model = "t_t_tv",
                           spline = list("age2", 2),
                           chains = 2,
                           cores = 2,
                           iter = 2000,
                           warmup = 1000,
                           control = list(adapt_delta = 0.999, 
                                          max_treedepth = 15)
                         )
print(fit_t_t_tv, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta0", "beta"))
traceplot(fit_t_t_tv, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta0", "beta"))

## fit nor t with tv dof
fit_nor_t_tv <- fit_ld(fixed = distance ~ age2, 
                     random = ~ age2, 
                     data = Orthodont,
                     id = "Subject",
                     model = "nor_t_tv",
                     spline = list("age2", 1),
                     chains = 2,
                     cores = 2,
                     iter = 2000,
                     warmup = 1000,
                     control = list(adapt_delta = 0.999, 
                                    max_treedepth = 15)
)
print(fit_nor_t_tv, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta"))
traceplot(fit_nor_t_tv, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta"))


## fit using nlme package
lmefit <- nlme::lme(distance ~ age2, data = Orthodont, random = ~ age2|Subject)
summary(lmefit)
