
 # loading the packages to be used
 #library(rstudioapi) 
 #library(ggplot2)
 library(JM)
 #library(survminer)
 library(robjm)

 # set the wd to where the codes are located
 #setwd(dirname(getActiveDocumentContext()$path))
 
 # source the supp. functions
 #source("Supplementary_functions.R")
 #source("Simulate_data.R")
 
 # check the data
 # g <- ggplot(sim_nor_nor$repeat_data, aes(time, Y, group = id))
 # g + geom_line() + facet_grid(. ~ event)
 # 
 # summary(unlist(with(sim_data$repeat_data, tapply(time, id, diff))))
 # 
 # summary(as.numeric(table(sim_data$repeat_data$id)))
 # 
 # km_fit <- survfit(Surv(stime, event) ~ 1, data = sim_data$base_data)
 # ggsurvplot(km_fit, data = sim_data$base_data, risk.table = TRUE)
 
 # nor nor
 sim_nor_nor <- simulate_data(model = "tv", 
                              nsubj = 100, 
                              av_n_i = 5,
                              eta = c(0.2, 0.5),
                              beta = c(0, 0, 0)) 
                              
 fit_nor_nor <- fit_jm(fixed_long = Y ~ time, 
                       random_long = ~ time, 
                       fixed_surv = cbind(stime, event) ~ c, 
                       data_long = sim_nor_nor$repeat_data,
                       data_surv = sim_nor_nor$base_data,
                       id_long = "id",
                       id_surv = "id",
                       model = "nor_nor",
                       timeVar = "time",
                       chains = 2,
                       cores = 2,
                       iter = 2000,
                       warmup = 1000,
                       bh = "weibull",
                       #pars = c("alpha", "Sigma", "sigma_Z", "sigmasq", "log_lambda", "lambda", "nu", "log_nu", "omega", "eta1", "eta2", "log_lik"),
                         deriv = list(deriv_fixed_formula = ~ 1,
                                      deriv_alpha_ind = 2, 
                                      deriv_random_formula = ~ 1, 
                                      deriv_B_ind = 2),
                       control = list(adapt_delta = 0.8)
                       )
 print(fit_nor_nor$res, pars = c("alpha", "Sigma", "sigmasq", "log_lambda", "log_nu", "omega", "eta1", "eta2"))
 traceplot(fit_nor_nor$res, pars = c("alpha", "Sigma", "sigmasq", "lambda", "nu", "omega", "eta1", "eta2"))

 # log_lik1 <- extract_log_lik(fit_nor_nor$res, merge_chains = FALSE)
 # rel_n_eff <- relative_eff(exp(log_lik1))
 # loo(log_lik1, r_eff = rel_n_eff, cores = 2)
 
 lme_fit <- lme(fixed = Y ~ time, random = ~ time|id, data = sim_nor_nor$repeat_data)
 summary(lme_fit)
 
 cox_fit <- coxph(Surv(stime, event) ~ c, data = sim_nor_nor$base_data, x = T)
 summary(cox_fit)
 
 dform <- list(fixed = ~ 1, indFixed = 2, random = ~ 1, indRandom = 2)
 
 joint_fit <- jointModel(lme_fit, cox_fit, timeVar = "time",
                         method = "weibull-PH-aGH",verbose = T, 
                         #derivForm = dform, 
                         #parameterization = "both",
                         iter.EM = 1000)
 summary(joint_fit)$"CoefTable-Long"
 summary(joint_fit)$"CoefTable-Event"
 summary(joint_fit)$"D"
 summary(joint_fit)$sigma^2
 
 sim_nor_nor_test <- simulate_data(model = "nor_nor",  
                                   nsubj = 4,
                                   av_n_i = 5,
                                   eta = c(0.2, 0.5))
 
 fore_nor_nor <- predSurv_jm(object = fit_nor_nor, 
                             newdata = sim_nor_nor_test$repeat_data, 
                             forecast = list(h = 5, n = 1),
                             B_control = list(iter = 30, warmup = 15, init = 0,
                                              nsel_b = 1),
                             batch_control = list(size = 2, cores = 2),
                             return_bsamples = TRUE)
 fore_nor_nor$output 
 sim_nor_nor_test$base_data
 
 sim_mod1 <- simulate_data(model = "mod1")
 fit_t_t_mod1 <- fit_jm(fixed_long = Y ~ time, 
                      random_long = ~ time,  
                      fixed_surv = cbind(stime, event) ~ c, 
                      data_long = sim_mod1$repeat_data,
                      data_surv = sim_mod1$base_data,
                      id_long = "id",
                      id_surv = "id",
                      bh = "spline",
                      model = "t_t_mod1",
                      timeVar = "time",
                      chains = 2,
                      cores = 2,
                      iter = 2000,
                      warmup = 1000,
                      control = list(adapt_delta = 0.99)
                      )
 print(fit_t_t_mod1, pars = c("alpha", "Sigma", "sigmasq", "phi", "lambda", "nu", "omega", "eta")) 
 
 sim_mod2 <- simulate_data(model = "mod2")
 fit_t_t_mod2 <- fit_jm(fixed_long = Y ~ time, 
                        random_long = ~ time,  
                        fixed_surv = cbind(stime, event) ~ c, 
                        data_long = sim_mod2$repeat_data,
                        data_surv = sim_mod2$base_data,
                        id_long = "id",
                        id_surv = "id",
                        model = "t_t_mod2",
                        bh = "spline",
                        chains = 2,
                        cores = 2,
                        iter = 2000,
                        warmup = 1000,
                        control = list(adapt_delta = 0.99))
 print(fit_t_t_mod2, pars = c("alpha", "Sigma", "phi", "sigmasq", "phi", "delta","lambda", "nu", "omega", "eta")) 
 traceplot(fit_t_t_mod2, pars = c("alpha", "Sigma", "phi", "sigmasq", "phi", "delta","lambda", "nu", "omega", "eta")) 
 
 sim_mod3 <- simulate_data(model = "mod3", 
                           nsubj = 100, 
                           av_n_i = 5,
                           eta = c(0.2, 0))
 fit_t_t_mod3 <- fit_jm(fixed_long = Y ~ time, 
                        random_long = ~ time,  
                        fixed_surv = cbind(stime, event) ~ c, 
                        data_long = sim_mod3$repeat_data,
                        data_surv = sim_mod3$base_data,
                        id_long = "id",
                        id_surv = "id",
                        model = "t_t_mod3",
                        timeVar = "time",
                        bh = "spline",
                        # deriv = list(deriv_fixed_formula = ~ 1,
                        #              deriv_alpha_ind = 2, 
                        #              deriv_random_formula = ~ 1, 
                        #              deriv_B_ind = 2),
                        chains = 2,
                        cores = 2,
                        iter = 1000,
                        warmup = 500,
                        control = list(adapt_delta = 0.9)
                        )
 print(fit_t_t_mod3$res, pars = c("alpha", "Sigma", "phi", "sigmasq", "phi", "delta","lambda", "nu", "omega", "eta")) 
 traceplot(fit_t_t_mod3$res, pars = c("alpha", "Sigma", "phi", "sigmasq", "phi", "delta","lambda", "nu", "omega", "eta")) 
 
 sim_mod3_test <- simulate_data(model = "mod3", 
                                nsubj = 4, 
                                av_n_i = 5,
                                eta = c(0.2, 0.6))
 fore_mod3 <- predSurv_jm(object = fit_t_t_mod3, 
                          newdata = sim_mod3_test$repeat_data, 
                          forecast = list(h = 5, n = 5),
                          B_control = list(adapt_delta = 0.9, iter = 500, warmup = 250),
                          batch_control = list(size = 2, cores = 2))
 
 sim_tv <- simulate_data(nsubj = 200, model = "tv", beta = c(-2, -1, 1))
 fit_t_t_tv <- fit_jm(fixed_long = Y ~ time, 
                        random_long = ~ time,  
                        fixed_surv = cbind(stime, event) ~ c, 
                        data_long = sim_tv$repeat_data,
                        data_surv = sim_tv$base_data,
                        id_long = "id",
                        id_surv = "id",
                      timeVar = "time",
                        model = "t_t_tv",
                        bh = "spline",
                        spline_tv = list("time", 2), 
                        chains = 2,
                        cores = 2,
                        iter = 2000,
                        warmup = 1000,
                        control = list(adapt_delta = 0.8, max_treedepth = 10)
                      )
 print(fit_t_t_tv$res, pars = c("alpha", "Sigma", "phi", "sigmasq", "phi", "delta0", "beta","lambda", "nu", "omega", "eta")) 
 traceplot(fit_t_t_tv, pars = c("alpha", "Sigma", "phi", "sigmasq", "phi", "delta0", "beta","lambda", "nu", "omega", "eta")) 
 
 a <- cbind(log(extract(fit_t_t_tv)$delta0), extract(fit_t_t_tv)$beta)
 colnames(a) <- c("log.delta0", "beta1", "beta2", "beta3") 
 pairs(a) 
