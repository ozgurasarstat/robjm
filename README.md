# robjm: An R Package for Robust Longitudinal and Joint Models

Install and load using the following lines:

library(devtools)  
install_github("ozgurasarstat/robjm")  
library(robjm)  

TO DO: 

add priors and Q into control list

x_T, x_quad, d_T, d_quad only work for random intercept and slope -- ask timeVar from the user?

baseline hazard is divided into 2 as piecewise constants

prediction functions

no need for 1/phi

but delta0 to exp

make delta = exp(2+exp(beta0 + a beta))