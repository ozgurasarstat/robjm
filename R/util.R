#' @title Creates a matrix of indicator values
#' @description Creates a matrix of zeros and ones, ones meaning being included in the interval
#' @param x A numeric vector for the time variable
#' @param knots A numeric vector that includes the knot locations
#'
pw_mat <- function(x, knots){
  
  matrix(unlist(lapply(x, in_interval, knots = knots)), nrow = length(x), byrow = TRUE)
  
}

#' @title CreateS a vector indicator values
#' @description Creates a vector of zeros and one, one meaning the value being included in the interval
#' @param x A numeric value for the time variable
#' @param knots A numeric vector that includes the knot locations
in_interval <- function(x, knots){
  
  intervals <- c(min(knots) - 10, knots, max(knots) + 10)#10 is arbitrary
  
  out <- rep(0, (length(knots) + 1))
    
  for(i in 1:(length(intervals) - 1)){
    if(x >= intervals[i] & x < intervals[i+1]){
      out[i] <- 1
    }
  }

  out
  
}

#' Function to summarise the probability estimates
#' @param x a numeric vector
prob_summary <- function(x){
  quants <- quantile(x, c(0.025, 0.975))
  out <- c(quants[1], mean(x), median(x), quants[2])
  return(out)
}

#' Function to combine prediction results
#' @param x a list
combine_pred <- function(x, nsubj = nsubj){

   if(nsubj == 1){
     samples <- fore_nor_nor[[1]]$ft_probs
     output  <- fore_nor_nor[[1]]$ft_table
   }else{
    iterations <- length(x)
    
    idlist  <- c()
    samples <- c()
    output  <- c()
    
    for(i in 1:iterations){
      
      iterations_nsubj <- length(x[[i]])
      
      for(j in 1:iterations_nsubj){
        
        if(length(x[[i]][[1]]) == 1){
          samples <- cbind(samples, x[[i]][[1]])
        }else{
          samples <- cbind(samples, x[[i]][[1]][[j]])
        }
        
      }
      
      idlist <- c(idlist, names(x[[i]][[1]]))
      output  <- rbind(output, x[[i]][[2]])
      
    }
    
    samples <- lapply(seq_len(ncol(samples)), function(i) samples[,i])
    names(samples) <- idlist  
  }
  
  return(list(samples = samples, output = output))

}