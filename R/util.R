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

