library(permute)
library(MASS)

##
## Exact Simulation of sampling design
##

ctrl <- how(maxperm = 10000000)
make.prob.perm <- function(N){
  perm.mat <- rbind(1:N, allPerms(N, control = ctrl))
  pi <- matrix(nrow = dim(perm.mat)[1], ncol = dim(perm.mat)[2])
  for(i in 1:dim(perm.mat)[1]){
	  for(k in 1:dim(perm.mat)[2]){
		  x <- perm.mat[i, 1:k]
  		fn <- ecdf(x)
	  	pi[i, perm.mat[i,k]] <- fn(x)[k]
  	}
  #	pi[i,] <- pi[i, ][perm.mat[i,]]
  }
  apply(pi, 2, mean)
}

make.HN <- function(N){ # calculate the harmonic number H_n <- note this is an approximation
	-digamma(1) + digamma(N+1)
}

make.pi <- function(N){
	if(floor(N) < N){
		"N must be an integer"
	} else if(N < 2){
		"N must be greater than 2"
	} else if(N == 2){
		bias <- c(1 / 4, 0)
	} else {
		HN <- make.HN(N)
		bias <- vector(length = N)
		for(j in 1:(N-1)){
			#bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2) #* (N - j)
			#bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
			#bias[j] <- (1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j) #This works well
			bias[j] <- (HN - 1) / (N * (N - 1)) * (N - j)
		}
		bias[N] <- 0
		return(bias)
	}
}

##
## Analytic Formula for sampling weights
##

##
## \hat{\pi}_i - \pi_i  = \frac{H_N - 1} {N (N - 1)} * (N - i)
##



## Works for N=1:8 relatively quick, slows down dramatically for N > 8
N <- 4
## Works fast
make.pi(N) + 1:N / N
 ## Works slow
make.prob.perm(N)
##
make.pi(N) + 1:N / N - make.prob.perm(N)

