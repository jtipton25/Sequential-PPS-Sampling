##
## adaptive sampling using a directional linear neighborhood
## 

N <- 100 #population size
W <- matrix(0, nrow = N, ncol = N)
W[cbind(1:(N - 1), 2:N)] <- 1
W


