##
## First try to estimate pi_kl for known ecdf sampling
##

##
## Need to automate the calculations to try and find patterns...
##

##
## libraries and functions
##

library(permute)
library(MASS)

make.pi <- function(N){
	if(floor(N) < N){
		"N must be an integer"
	} else if(N < 2){
		"N must be greater than 2"
	} else if(N == 2){
		bias <- c(1 / 4, 0)
	} else {
		bias <- vector(length = N)
		for(j in 1:(N-1)){
			#bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2) #* (N - j)
			#bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
			bias[j] <- (1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
		}
		bias[N] <- 0
		return(bias)
	}
}

make.ecdf.sample <- function(s, N){
	x <- 1:N
	fn <- ecdf(x)
	fn(x)
	samp <- which(rbinom(N, 1, prob = fn(x)) == 1)
	x[samp]
}

ctrl <- how(maxperm = 10000000)

make.pi.kl <- function(N){
	perm.mat <- rbind(1:N, allPerms(N, control = ctrl))
	pi <- matrix(nrow = dim(perm.mat)[1], ncol = dim(perm.mat)[2])
	pi.kl <- matrix(nrow = N, ncol = N)
	for(i in 1:dim(perm.mat)[1]){
		for(k in 1:dim(perm.mat)[2]){
			x <- perm.mat[i, 1:k]
			fn <- ecdf(x)
			pi[i, perm.mat[i,k]] <- fn(x)[k]
		}
		#	pi[i,] <- pi[i, ][perm.mat[i,]]
	}
	for(k in 1:N){
		for(l in k:N){
			pi.kl[k, l] <- sum(pi[, k] * pi[, l]) / factorial(N)
		}
	}
	return(pi.kl)
}

##
## Automate design pi.kl weights
##

##
## N = 3
##
N <- 3
pi.kl <- make.pi.kl(N)
pi.kl
pi.kl[1, 2] # pi.kl for design
pi.kl[1, 2] - (1*2) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 2] - (1*2) /(N^2))

##
## N = 4
##
N <- 4
pi.kl <- make.pi.kl(N)
pi.kl

pi.kl[1, 2] # pi.kl for design
pi.kl[1, 2] - (1*2) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 2] - (1*2) /(N^2))

pi.kl[1, 3] # pi.kl for design
pi.kl[1, 3] - (1*3) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 3] - (1*3) /(N^2))

pi.kl[2, 3] # pi.kl for design
pi.kl[2, 3] - (2*3) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[2, 3] - (2*3) /(N^2))

##
## N = 5
##

N <- 5
pi.kl <- make.pi.kl(N)
pi.kl

pi.kl[1, 2] # pi.kl for design
pi.kl[1, 2] - (1*2) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 2] - (1*2) /(N^2))

pi.kl[1, 3] # pi.kl for design
pi.kl[1, 3] - (1*3) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 3] - (1*3) /(N^2))

pi.kl[1, 4] # pi.kl for design
pi.kl[1, 4] - (1*4) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 4] - (1*4) /(N^2))

pi.kl[2, 3] # pi.kl for design
pi.kl[2, 3] - (2*3) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[2, 3] - (2*3) /(N^2))

pi.kl[2, 4] # pi.kl for design
pi.kl[2, 4] - (2*4) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[2, 4] - (2*4) /(N^2))

pi.kl[3, 4] # pi.kl for design
pi.kl[3, 4] - (3*4) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[3, 4] - (3*4) /(N^2))

##
## N = 6
##

N <- 6
pi.kl <- make.pi.kl(N)
pi.kl

pi.kl[1, 2] # pi.kl for design
pi.kl[1, 2] - (1*2) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 2] - (1*2) /(N^2))

pi.kl[1, 3] # pi.kl for design
pi.kl[1, 3] - (1*3) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 3] - (1*3) /(N^2))

pi.kl[1, 4] # pi.kl for design
pi.kl[1, 4] - (1*4) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 4] - (1*4) /(N^2))

pi.kl[1, 5] # pi.kl for design
pi.kl[1, 5] - (1*5) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[1, 5] - (1*5) /(N^2))

pi.kl[2, 3] # pi.kl for design
pi.kl[2, 3] - (2*3) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[2, 3] - (2*3) /(N^2))

pi.kl[2, 4] # pi.kl for design
pi.kl[2, 4] - (2*4) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[2, 4] - (2*4) /(N^2))

pi.kl[2, 5] # pi.kl for design
pi.kl[2, 5] - (2*5) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[2, 5] - (2*5) /(N^2))

pi.kl[3, 4] # pi.kl for design
pi.kl[3, 4] - (3*4) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[3, 4] - (3*4) /(N^2))

pi.kl[3, 5] # pi.kl for design
pi.kl[3, 5] - (3*5) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[3, 5] - (3*5) /(N^2))

pi.kl[4, 5] # pi.kl for design
pi.kl[4, 5] - (4*5) /(N^2) # bias of pi.kl for design relative to known ecdf
as.fractions(pi.kl[4, 5] - (4*5) /(N^2))



##
## Try to come up with an analytic solution
##

make.HN <- function(N){ # calculate the harmonic number H_n
	-digamma(1) + digamma(N+1)
}

make.pi.kl <- function(N){
	if(floor(N) < N){
		"N must be an integer"
	} else if(N < 2){
		"N must be greater than 2"
	} else if(N == 2){
		bias <- c(1 / 4, 0)
	} else {
		bias <- matrix(nrow = N, ncol = N)
		HN <- make.HN(N)
		for(k in 1:N){
			for(l in 1:N){
        #bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2) #* (N - j)
			  #bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
			  #bias[k, l] <- - (1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - k) * (N - l) +
			  #	(1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - k)  + 
			  #	(1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - l)
				
				#bias[k, l] <- - (1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - k) * (N - l) +
				#	(HN - 1) / (N * (N - 1)) * (N - k) + 
	      #  (HN - 1) / (N * (N - 1)) * (N - l)
				
				bias[k, l] <- (N - k) * (HN - N) / (N - 1) + (N - l) * (HN - N) / (N - 1) - (N - k) * (N - l) * (HN - N) / (N - 1)
			}
		}
		return(bias)
	}
}
	
N <- 3
pi.kl.bias <- make.pi.kl(N) ## this has the right behavior along the last column and row...
pi.kl.bias
#as.fractions(pi.kl.bias)
make.pi(N) ## check last row and column of above
#as.fractions(make.pi(N))
















