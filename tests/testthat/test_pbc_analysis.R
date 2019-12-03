
#devtools::install_github("ozgurasarstat/robjm")
library("robjm")

data(pbc2)

pbc2 <- pbc2[complete.cases(pbc2$alkaline), ]
pbc2_base <- pbc2[!duplicated(pbc2$id), ]

#setwd("C:\\Users\\ozgur.asar\\Documents\\PBC analysis")
#setwd("C:\\Users\\Ozgur Asar\\Documents\\Dropbox\\PBC analysis")

## normal normal model
fit_nor_nor <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                      random_long = ~ year, 
                      fixed_surv = cbind(years, status2) ~ age + drug2, 
                      data_long = pbc2,
                      data_surv = pbc2_base,
                      id_long = "id",
                      id_surv = "id",
                      model = "nor_nor",
                      timeVar = "year",
                      chains = 4,
                      cores = 4,
                      iter = 6000,
                      warmup = 3000,
                      bh = "weibull",
                      pars = c("alpha", "Sigma", "sigma_Z", "sigmasq", 
                               "log_lambda", "lambda", "nu", "log_nu", "omega", "eta", "B"),
                      control = list(adapt_delta = 0.6))

print(fit_nor_nor$res, pars = c("alpha", "Sigma", "sigmasq", "log_lambda", 
                                "log_nu", "omega", "eta"))

traceplot(fit_nor_nor$res, pars = c("alpha", "Sigma", "sigmasq", "log_lambda", 
                                    "log_nu", "omega", "eta"))

#saveRDS(fit_nor_nor, "fit_nor_nor.rds")

## t t model 1
fit_t_t_mod1 <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                       random_long = ~ year, 
                       fixed_surv = cbind(years, status2) ~ age + drug2, 
                       data_long = pbc2,
                       data_surv = pbc2_base,
                       id_long = "id",
                       id_surv = "id",
                       model = "t_t_mod1",
                       timeVar = "year",
                       chains = 4,
                       cores = 4,
                       iter = 6000,
                       warmup = 3000,
                       bh = "weibull",
                       pars = c("alpha", "Sigma", "sigma_Z", "sigmasq", "phi",
                                "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                       control = list(adapt_delta = 0.6, max_treedepth = 10))

print(fit_t_t_mod1$res, pars = c("alpha", "Sigma", "sigmasq", "phi", 
                                 "log_lambda", "log_nu", "omega", "eta"))
traceplot(fit_t_t_mod1$res, pars = c("alpha", "Sigma", "sigmasq", "phi", 
                                     "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_t_t_mod1, "fit_t_t_mod1.rds")

# ## t - t, mod2
# 
# fit_t_t_mod2 <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
#                        random_long = ~ year, 
#                        fixed_surv = cbind(years, status2) ~ age + drug2, 
#                        data_long = pbc2,
#                        data_surv = pbc2_base,
#                        id_long = "id",
#                        id_surv = "id",
#                        model = "t_t_mod2",
#                        timeVar = "year",
#                        chains = 4,
#                        cores = 4,
#                        iter = 2000,
#                        warmup = 1000,
#                        bh = "weibull",
#                        #pars = c("alpha", "Sigma", "sigmasq", "phi", "delta",
#                        #        "log_lambda", "log_nu", "omega", "eta"),
#                        control = list(adapt_delta = 0.8))
# 
# print(fit_t_t_mod2$res, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta",
#                              "log_lambda", "log_nu", "omega", "eta"))
# 
# traceplot(fit_t_t_mod2$res, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta",
#                                  "log_lambda", "log_nu", "omega", "eta"))
# 
# saveRDS(fit_t_t_mod2, "fit_t_t_mod2.rds")

## nor - t, mod2

fit_nor_t_mod2 <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                         random_long = ~ year, 
                         fixed_surv = cbind(years, status2) ~ age + drug2, 
                         data_long = pbc2,
                         data_surv = pbc2_base,
                         id_long = "id",
                         id_surv = "id",
                         model = "nor_t_mod2",
                         timeVar = "year",
                         chains = 4,
                         cores = 4,
                         iter = 6000,
                         warmup = 3000,
                         bh = "weibull",
                         pars = c("alpha", "Sigma", "sigma_Z", "sigmasq", "delta",
                                  "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                         control = list(adapt_delta = 0.6))

print(fit_nor_t_mod2$res, pars = c("alpha", "Sigma", "sigmasq",  "delta",
                                   "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_mod2$res, pars = c("alpha", "Sigma", "sigmasq", "delta",
                                       "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_mod2, "fit_nor_t_mod2.rds")

# ## t - t, mod3
# fit_t_t_mod3 <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
#                        random_long = ~ year, 
#                        fixed_surv = cbind(years, status2) ~ age + drug2, 
#                        data_long = pbc2,
#                        data_surv = pbc2_base,
#                        id_long = "id",
#                        id_surv = "id",
#                        model = "t_t_mod3",
#                        timeVar = "year",
#                        chains = 4,
#                        cores = 4,
#                        iter = 6000,
#                        warmup = 3000,
#                        bh = "weibull",
#                        pars = c("alpha", "Sigma", "sigmasq", "phi", "delta",
#                                 "log_lambda", "lambda", "log_nu", "nu", 
#                                 "omega", "eta", "B"),
#                        control = list(adapt_delta = 0.6, max_treedepth = 10))
# 
# print(fit_t_t_mod3$res, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta",
#                              "log_lambda", "log_nu", "omega", "eta"))
# 
# traceplot(fit_t_t_mod3$res, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta",
#                              "log_lambda", "log_nu", "omega", "eta"))
# 
# saveRDS(fit_t_t_mod3, "fit_t_t_mod3.rds")

## nor - t, mod3
fit_nor_t_mod3 <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                         random_long = ~ year, 
                         fixed_surv = cbind(years, status2) ~ age + drug2, 
                         data_long = pbc2,
                         data_surv = pbc2_base,
                         id_long = "id",
                         id_surv = "id",
                         model = "nor_t_mod3",
                         timeVar = "year",
                         chains = 4,
                         cores = 4,
                         iter = 6000,
                         warmup = 3000,
                         bh = "weibull",
                         pars = c("alpha", "Sigma", "sigmasq", "delta",
                                  "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                         control = list(adapt_delta = 0.6, max_treedepth = 10))

print(fit_nor_t_mod3$res, pars = c("alpha", "Sigma", "sigmasq", "delta",
                                   "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_mod3$res, pars = c("alpha", "Sigma", "sigmasq", "delta",
                                       "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_mod3, "fit_nor_t_mod3.rds")

# ## t - t, time-varying
# fit_t_t_tv_2knots <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2,
#                             random_long = ~ year,
#                             fixed_surv = cbind(years, status2) ~ age + drug2,
#                        data_long = pbc2,
#                        data_surv = pbc2.id,
#                        id_long = "id",
#                        id_surv = "id",
#                        model = "t_t_tv",
#                        timeVar = "year",
#                        spline_tv = list("year", 2),
#                        chains = 4,
#                        cores = 4,
#                        iter = 6000,
#                        warmup = 3000,
#                        bh = "weibull",
#                        #pars = c("alpha", "Sigma", "sigma_Z", "sigmasq", "phi", "delta0", "beta",
#                        #         "log_lambda", "lambda", "log_nu", "nu", "omega", "eta"),
#                        control = list(adapt_delta = 0.6))
# 
# print(fit_t_t_tv_2knots$res, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta0", "beta",
#                                    "log_lambda", "log_nu", "omega", "eta"))
# 
# traceplot(fit_t_t_tv_2knots$res, pars = c("alpha", "Sigma", "sigmasq", "phi", "delta0", "beta",
#                                        "log_lambda", "log_nu", "omega", "eta"))
# 
# saveRDS(fit_t_t_tv_2knots, "fit_t_t_tv_2knots.rds")

## nor - t, time-varying
fit_nor_t_tv_1knot <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                             random_long = ~ year, 
                             fixed_surv = cbind(years, status2) ~ age + drug2,
                             data_long = pbc2,
                             data_surv = pbc2.id,
                             id_long = "id",
                             id_surv = "id",
                             model = "nor_t_tv",
                             timeVar = "year",
                             spline_tv = list("year", 1), 
                             chains = 4,
                             cores = 4,
                             iter = 6000,
                             warmup = 3000,
                             bh = "weibull",
                             pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta", 
                                      "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                             control = list(adapt_delta = 0.6, max_treedepth = 15))

print(fit_nor_t_tv_1knot$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                       "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_tv_1knot$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                           "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_tv_1knot, "fit_nor_t_tv_1knot.rds")

## nor - t, time-varying
fit_nor_t_tv_2knots <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                              random_long = ~ year, 
                              fixed_surv = cbind(years, status2) ~ age + drug2,
                              data_long = pbc2,
                              data_surv = pbc2.id,
                              id_long = "id",
                              id_surv = "id",
                              model = "nor_t_tv",
                              timeVar = "year",
                              spline_tv = list("year", 2), 
                              chains = 4,
                              cores = 4,
                              iter = 6000,
                              warmup = 3000,
                              bh = "weibull",
                              pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta", 
                                       "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                              control = list(adapt_delta = 0.6, max_treedepth = 15))

print(fit_nor_t_tv_2knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                        "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_tv_2knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                            "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_tv_2knots, "fit_nor_t_tv_2knots.rds")

## nor - t, time-varying - 3 knots
fit_nor_t_tv_3knots <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                              random_long = ~ year, 
                              fixed_surv = cbind(years, status2) ~ age + drug2,
                              data_long = pbc2,
                              data_surv = pbc2.id,
                              id_long = "id",
                              id_surv = "id",
                              model = "nor_t_tv",
                              timeVar = "year",
                              spline_tv = list("year", 3), 
                              chains = 4,
                              cores = 4,
                              iter = 6000,
                              warmup = 3000,
                              bh = "weibull",
                              pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta", 
                                       "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                              control = list(adapt_delta = 0.6, max_treedepth = 15))

print(fit_nor_t_tv_3knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                        "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_tv_3knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                            "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_tv_3knots, "fit_nor_t_tv_3knots.rds")

## nor - t, time-varying - 4 knots

fit_nor_t_tv_4knots <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                              random_long = ~ year, 
                              fixed_surv = cbind(years, status2) ~ age + drug2, 
                              data_long = pbc2,
                              data_surv = pbc2.id,
                              id_long = "id",
                              id_surv = "id",
                              model = "nor_t_tv",
                              timeVar = "year",
                              spline_tv = list("year", 4), 
                              chains = 4,
                              cores = 4,
                              iter = 6000,
                              warmup = 3000,
                              bh = "weibull",
                              pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta", 
                                       "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                              control = list(adapt_delta = 0.6, max_treedepth = 15))

print(fit_nor_t_tv_4knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                        "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_tv_4knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                            "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_tv_4knots, "fit_nor_t_tv_4knots.rds")

## nor - t, time-varying - 5 knots

fit_nor_t_tv_5knots <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                              random_long = ~ year, 
                              fixed_surv = cbind(years, status2) ~ age + drug2,
                              data_long = pbc2,
                              data_surv = pbc2.id,
                              id_long = "id",
                              id_surv = "id",
                              model = "nor_t_tv",
                              timeVar = "year",
                              spline_tv = list("year", 5), 
                              chains = 4,
                              cores = 4,
                              iter = 6000,
                              warmup = 3000,
                              bh = "weibull",
                              pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta", 
                                       "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                              control = list(adapt_delta = 0.6, max_treedepth = 10))

print(fit_nor_t_tv_5knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                        "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_tv_5knots$res, pars = c("alpha", "Sigma", "sigmasq", "delta0", "beta",
                                            "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_tv_5knots, "fit_nor_t_tv_5knots.rds")

## nor - t, time-varying dof and scale
fit_nor_t_tv_dof_scale_1knot <- fit_jm(fixed_long = log(alkaline) ~ age + year + drug2, 
                                       random_long = ~ year, 
                                       fixed_surv = cbind(years, status2) ~ age + drug2,
                                       data_long = pbc2,
                                       data_surv = pbc2.id,
                                       id_long = "id",
                                       id_surv = "id",
                                       model = "nor_t_tv_dof_scale",
                                       timeVar = "year",
                                       spline_tv = list("year", 1), 
                                       chains = 4,
                                       cores = 4,
                                       iter = 6000,
                                       warmup = 3000,
                                       bh = "weibull",
                                       pars = c("alpha", "Sigma", "sigma0_Z", "psi", "delta0", "beta", 
                                                "log_lambda", "lambda", "log_nu", "nu", "omega", "eta", "B"),
                                       control = list(adapt_delta = 0.6, max_treedepth = 15))

print(fit_nor_t_tv_dof_scale_1knot$res, pars = c("alpha", "Sigma", "sigma0_Z", "psi", "delta0", "beta",
                                                 "log_lambda", "log_nu", "omega", "eta"))

traceplot(fit_nor_t_tv_dof_scale_1knot$res, pars = c("alpha", "Sigma", "sigma0_Z", "psi", "delta0", "beta",
                                                     "log_lambda", "log_nu", "omega", "eta"))

#saveRDS(fit_nor_t_tv_dof_scale_1knot, "fit_nor_t_tv_dof_scale_1knot.rds")
