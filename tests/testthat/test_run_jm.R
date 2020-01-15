
#library(devtools)
#install_github("ozgurasarstat/robjm")
library("robjm")

data(aids)
data(aids.id)

aids$drug2 <- ifelse(aids$drug == "ddC", 0, 1)
aids.id$drug2 <- ifelse(aids.id$drug == "ddC", 0, 1)

id_first_100 <- aids.id$patient[1:100]

long_data <- aids[aids$patient %in% id_first_100, ]
surv_data <- aids.id[aids.id$patient %in% id_first_100, ]

## normal normal model
fit_nor_nor <- fit_jm(fixed_long = CD4 ~ obstime, 
                      random_long = ~ obstime, 
                      fixed_surv = cbind(Time, death) ~ drug2, 
                      data_long = long_data,
                      data_surv = surv_data,
                      id_long = "patient",
                      id_surv = "patient",
                      model = "nor_t_tv_dof_scale",
                      timeVar = "obstime",
                      bh = "weibull",
                      spline_tv = list("obstime", 3),
                      chains = 2,
                      cores = 2,
                      iter = 2000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.9, max_treedepth = 15) 
                      )
print(fit_nor_nor, pars = c("alpha", "Sigma", "sigmasq", "lambda", "nu", "omega", "eta"))
traceplot(fit_nor_nor, 
          pars = c("alpha", "Sigma", "sigmasq", "zeta", "omega", "eta"),
          inc_warmup = FALSE)
#pairs(fit_nor_nor, pars = c("alpha", "Sigma", "sigmasq", "zeta", "omega", "eta"))

## t t model 1
fit_t_t_mod1 <- fit_jm(fixed_long = CD4 ~ obstime, 
                      random_long = ~ obstime, 
                      fixed_surv = cbind(Time, death) ~ drug2, 
                      data_long = long_data,
                      data_surv = surv_data,
                      id_long = "patient",
                      id_surv = "patient",
                      model = "t_t_mod1",
                      timeVar = "obstime",
                      chains = 2,
                      cores = 2,
                      iter = 2000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99)
                      )
print(fit_t_t_mod1, pars = c("alpha", "Sigma", "sigmasq", "phi", "zeta", "omega", "eta"))
traceplot(fit_t_t_mod1, 
          pars = c("alpha", "Sigma", "sigmasq", "phi", "zeta", "omega", "eta"),
          inc_warmup = FALSE)
pairs(fit_t_t_mod1, pars = c("alpha", "Sigma", "sigmasq", "phi", "zeta", "omega", "eta"))

fit_t_t_mod2 <- fit_jm(fixed_long = CD4 ~ obstime, 
                       random_long = ~ obstime, 
                       fixed_surv = cbind(Time, death)~drug2, 
                       data_long = long_data,
                       data_surv = surv_data,
                       id_long = "patient",
                       id_surv = "patient",
                       model = "t_t_mod2",
                       chains = 4,
                       cores = 4,
                       iter = 2000,
                       warmup = 1000,
                       control = list(adapt_delta = 0.99)
                       )
print(fit_t_t_mod2, pars = c("alpha", "Sigma", "phi", "sigmasq", "delta", "zeta", "omega", "eta"))

fit_t_t_mod3 <- fit_jm(fixed_long = CD4 ~ obstime, 
                       random_long = ~ obstime, 
                       fixed_surv = cbind(Time, death)~drug2, 
                       data_long = long_data,
                       data_surv = surv_data,
                       id_long = "patient",
                       id_surv = "patient",
                       model = "t_t_mod3",
                       chains = 4,
                       cores = 4,
                       iter = 2000,
                       warmup = 1000,
                       control = list(adapt_delta = 0.99)
                       )
print(fit_t_t_mod2, pars = c("alpha", "Sigma", "phi", "sigmasq", "delta", "zeta", "omega", "eta"))

fit_t_t_tv <- fit_jm(fixed_long = CD4 ~ obstime, 
                       random_long = ~ obstime, 
                       fixed_surv = cbind(Time, death)~drug2, 
                       data_long = long_data,
                       data_surv = surv_data,
                       id_long = "patient",
                       id_surv = "patient",
                       model = "t_t_tv",
                       spline_tv = list("obstime", 2), 
                       chains = 4,
                       cores = 4,
                       iter = 2000,
                       warmup = 1000,
                       control = list(adapt_delta = 0.99)
                     )
print(fit_t_t_tv, pars = c("alpha", "Sigma", "phi", "sigmasq", "delta0", "beta", "zeta", "omega", "eta"))

fit_nor_t_fixed_dof_mod3 <- fit_jm(fixed_long = CD4 ~ obstime, 
                       random_long = ~ obstime, 
                       fixed_surv = cbind(Time, death)~drug2, 
                       data_long = long_data,
                       data_surv = surv_data,
                       id_long = "patient",
                       id_surv = "patient",
                       model = "nor_t_fixed_dof_mod3",
                       deriv = list(deriv_fixed_formula = ~ 1,
                                    deriv_alpha_ind = 2, 
                                    deriv_random_formula = ~ 1, 
                                    deriv_B_ind = 2),
                       timeVar = "obstime",
                       chains = 4,
                       cores = 4,
                       iter = 2000,
                       warmup = 1000,
                       control = list(adapt_delta = 0.99))
print(fit_nor_t_fixed_dof_mod3$res, 
      pars = c("alpha", "Sigma", "sigmasq",  "omega", "eta"))
