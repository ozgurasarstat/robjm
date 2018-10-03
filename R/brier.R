
#'
#'@param pred A numericvector for failure probabilities
#'
 brier <- function(pred, T, E, lm_time, h){

   km_fit <- survival::survfit(survival::Surv(T, E) ~ 1)
   
   G1 <- summary(km_fit, times = T)$surv
   G2 <- summary(km_fit, times = lm_time + h)$surv
   
   T_leq <- ifelse(T <= lm_time + h, 1, 0)
   T_gr  <- ifelse(T > lm_time + h, 1, 0)
   
   expr1 <- (pred - T_leq) %>% '^'(2)
   expr2 <- (T_leq * E)/G1 + T_gr/G2
 
   out <- mean(expr1 * expr2)
   return(out)
   
 }