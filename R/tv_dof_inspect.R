tv_dof_inspect <- function(time_range = c(0, 5),
                           incr = 0.01,
                           nknots = 2,
                           delta0 = 5,
                           beta = c(-2, 0, 2),
                           plot = TRUE,
                           return_d = TRUE){
  
  beta <- as.matrix(beta, ncol = 1)
  
  d <- data.frame(time = seq(time_range[1], time_range[2], incr))
  
  d$a <- splines::ns(d$time, df = (nknots + 1))
  
  d$delta_ij <- delta0 * exp(d[, -1] %*% beta)
  
  if(plot){
    g <- ggplot(d, aes(x = time, y = delta_ij)) + geom_line() + labs(x = "t", y = bquote(delta(t)))
    print(g)   
  }
  
  if(return_d){
    return(d)
  }
  
}