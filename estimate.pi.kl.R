##
## First try to estimate pi_kl for known ecdf sampling
##

##
## Need to automate the calculations to try and find patterns...
##

##
## libraries and functions
##

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

library(permute)
library(MASS)

##
## Exact Simulation of sampling design
##

make.ecdf.sample <- function(s, N){
  x <- 1:N
  fn <- ecdf(x)
  fn(x)
  samp <- which(rbinom(N, 1, prob = fn(x)) == 1)
  x[samp]
}

iter <- 10000

N <- 4
out.ecdf <- lapply(1:iter, make.ecdf.sample, N)

idx <- c()
out.list <- vector('list', length = iter)
for(s in 1:iter){
	temp <- matrix(nrow = N, ncol = N)
	for(k in 1:N){
		for(l in 1:N){
			temp[k, l] <- sum(c(k, l) %in% out.ecdf[[s]])
			if(s == iter){
				idx <- cbind(idx, c(k, l))
			}
		}
	}
	out.list[[s]] <- which(temp == 2)
}

prob.kl <- vector(length = N^2)

for(k in 1:(N^2)){
	temp <- vector(length = iter)
	for(s in 1:iter){
		temp[s] <- k %in% out.list[[s]]
	}
	prob.kl[k] <- sum(temp) / iter
}

idx
prob.kl

##
## N = 3
##

pi.12 <- (1*1 + 1*2/3 + 1*1/2 + 1/3*1 + 1/2*2/3 + 1/2*1/3) / factorial(3)
pi.12 # exact pi.kl for sampling design
(1*2)/(3^2) # exact pi.kl for ecdf sampling
pi.12 - (1*2)/(3^2) # bias for pi.kl between design and ecdf
as.fractions(pi.12 - (1*2)/(3^2))

##
## N = 4
##

pi.12 <- (1*1 + 1*1 + 1*2/3 + 1*1/2 + 1*2/3 + 1*1/2 + 1/2*1 + 1/2*1 + 1/3*1 + 1/4*1 + 1/4*1 + 1/3*1 + 1/2*2/3 + 1/2*1/2 + 1/3*1/2 + 1/4*1/2 + 1/3*1/2 + 1/4*1/3 + 1/2*2/3 + 1/2*1/2 + 1/3*1/2 + 1/4*1/2 + 1/3*1/2 + 1/4*1/3) / factorial(4)
pi.12 # exact pi.kl for sampling design
(1*2)/(4^2) # exact pi.kl for ecdf sampling
pi.12 - (1*2)/(4^2) # bias for pi.kl between design and ecdf
as.fractions(pi.12 - (1*2)/(4^2))

pi.13 <- (1*1 + 1*3/4 + 1*1 + 1*1 + 1*3/4 + 1*2/3 + 1/2*3/4 + 1/2*1 + 1/3*1 + 1/4*1 + 1/4*2/3 + 1/3*3/4 + 1/2*1 + 1/2*1 + 1/3*1 + 1/4*1 + 1/3*1 + 1/4*1 + 1/2*3/4 + 1/2*2/3 + 1/3*3/4 + 1/4*2/3 + 1/3*1/2 + 1/4*1/2) / factorial(4)
pi.13 # exact pi.kl for sampling design
(1*3)/(4^2) # exact pi.kl for ecdf sampling
pi.13 - (1*3)/(4^2) # bias for pi.kl between design and ecdf
as.fractions(pi.13 - (1*3)/(4^2))

#pi.23 <- (1*1 + 1*3/4 + 2/3*1 + 1/2*1 + 2/3*3/4 + 1/2*2/3 + 1*3/4 + 1*1 + 1*1 + 1*1 + 1*2/3 + 1*3/4 + 2/3*1 + 1/2*1 + 1/2*1 + 1/2*1 + 1/2*1 + 1/3*1 + 2/3*3/4 + 1/2*2/3 + 1/2*3/4 + 1/2*2/3 + 1/2*1/3 + 1/3*1/2) / factorial(4)
pi.23 <- (1*1 + 1*3/4 + 2/3*1 + 1/2*1 + 2/3*3/4 + 1/2*2/3 + 1*3/4 + 1*1 + 1*1 + 1*1 + 1*2/3 + 1*3/4 + 2/3*1 + 1/2*1 + 1/2*1 + 1/2*1 + 1/2*1 + 1/3*1 + 2/3*3/4 + 1/2*2/3 + 1/2*3/4 + 1/2*2/3 + 1/2*1/2 + 1/3*1/2) / factorial(4)
pi.23 # exact pi.kl for sampling design
(2*3)/(4^2) # exact pi.kl for ecdf sampling
pi.23 - (2*3)/(4^2) # bias for pi.kl between design and ecdf
as.fractions(pi.23 - (2*3)/(4^2))

as.fractions(c(pi.12, pi.13, pi.23) - c((1*2)/(4^2), (1*3)/(4^2), (2*3)/(4^2))) # bias for pi.kl between design and ecdf
c(23*6, 49*3, 131)

##
## N = 5 <- very long and complicated, will see if there is any pattern, perhaps a denominator of (N!)^2 in the bias?
##





##
## Automate design pi.kl weights
##

N <- 4

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

make.pi.kl <- function(N){
	if(floor(N) < N){
		"N must be an integer"
	} else if(N < 2){
		"N must be greater than 2"
	} else if(N == 2){
		bias <- c(1 / 4, 0)
	} else {
		bias <- matrix(nrow = N, ncol = N)
		for(k in 1:N){
			for(l in 1:N){
        #bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2) #* (N - j)
			  #bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
			  bias[k, l] <- - (1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - k) * (N - l) +
			  	(1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - k)  + 
			  	(1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - l)
			}
		}
		return(bias)
	}
}
	
N <- 3
pi.kl.bias <- make.pi.kl(N) ## this has the right behavior along the last column and row...
pi.kl.bias
as.fractions(pi.kl.bias)
make.pi(N) ## check last row and column of above
as.fractions(make.pi(N))
