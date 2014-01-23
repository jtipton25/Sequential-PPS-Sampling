##
## libraries and functions
##

source('sequential.functions.R')

##
## Simulate dbh
##

N <- 1000 # finite population size
n <- 100 # expected sample size

##
## Start with a population that is a mixture of Gamma Distributions
##

alpha <- c(2, 4, 6, 8) # gamma mixture parameter
beta <- c(12, 10, 8, 6) # gamma mixture parameter

dbh <- make.gamma.mixture(N, alpha, beta)

##
## Plot DBH
##

layout(matrix(1:4, 2, 2))
for(i in 1:4){ # plot the four distibutions to mix
	curve(dgamma(x, alpha[i], beta[i]), from = 0, to = 4)
}

par(mfrow = c(1, 1))
#pdf(file = 'dbh1.pdf')
hist(dbh, breaks = 20)
#dev.off()

## 
## True mean and variance for DBH
##

dbh.mn <- mean(dbh)
dbh.var <- var(dbh)

##
## Simulate Biomass
##

b0 <- 5
b1 <- 2
s2 <- 1/4

bio <- make.sim.biomass(dbh, b0, b1, s2)
bio.mn <- mean(bio)
bio.var <- var(bio)

##
## Plot Data
##

#make.model.plot(dbh, bio, file = "fullModel.pdf")
make.model.plot(dbh, bio)

##
## empirical estimator results
##

make.est.result <- function(rep, dbh, bio, n = N / 2, method = 'design'){
	out <- make.samp(dbh, bio, n, method)
	n.samp <- length(out$samp)
  if(method == 'srs'){
  	return(mean(out$bio))
  }
	if(method == 'pps'){
		return(1 / N * sum(out$bio / out$p[out$samp]))
	}
	if(method == 'ecdf'){
		return(1 / N * sum(out$bio / out$p[out$samp]))
	}
	if(method == 'design'){
		z <- vector(length = n.samp - 1)
		q <- 1 / (N - 1:(n.samp - 1))
		for(i in 2:n.samp){
			z[i - 1] <- sum(out$bio[1:(i - 1)]) + out$bio[i] / q[i - 1]
		}
		return(1 / (N * n) * (N*out$bio[1] + sum(z)))
	}
}

library(parallel)
library(snowfall)
library(rlecuyer)
cps <- detectCores()
sfInit(parallel = TRUE, cpus = cps)
sfClusterSetupRNG() 
sfExportAll()

temp.srs <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = N / 2, method = 'srs') 
summary(temp.srs)
mean(temp.srs)
var(temp.srs)

temp.pps <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = N / 2, method = 'pps') 
summary(temp.pps)
mean(temp.pps)
var(temp.pps)

temp.ecdf <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = N / 2, method = 'ecdf') 
summary(temp.ecdf)
mean(temp.ecdf)
var(temp.ecdf)

temp.design <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = N / 2, method = 'design') 
summary(temp.design)
mean(temp.design)
var(temp.design)


bio.mn
mean(temp.srs)
mean(temp.ecdf)
mean(temp.pps)
mean(temp.design)

var(temp.srs)
var(temp.ecdf)
var(temp.pps)
var(temp.design)


dbh.mn


# make.two.stage.samp <- function(rep, dbh, bio, n, n0, N0, method = 'design'){
# 	N <- length(dbh)
#	if(method == 'srs'){
#		p <- n * rep(1 / N, N)
#		samp <- sample(1:N, n, prob = p)
#	}
#	if(method == 'pps'){
#		p <- n * dbh / sum(dbh)
#		samp <- sample(1:N, n, prob = p)
#	}
#	if(method == 'ecdf'){
#		p <- 2 * n / N * ecdf(dbh)(dbh)
#		samp <- which(rbinom(1:N, 1, p) == 1)
#	}
# # 	if(method == 'design'){
# 		dbh <- sample(dbh)
# 		idx0 <- sample(1:N0, n0)
# 		idx <- sample(N0:N)
# 		dbh.srs <- dbh[idx0]	
# 		p <- vector(length = N - N0)
# 		for(k in (N0 + 1):N){
# 			x <- dbh[idx[1:(k - N0)]]
# 			fn <- ecdf(x)
# 			p[k - N0] <- fn(x)[k - N0]
# 		}
# 		#p <- p * 2 * n / N # potential sample size adjustment
# 		samp <- which(rbinom(N - N0, 1, p) == 1)
# 		#p <- make.pi(N)
# 	}  
# 	prob <- p[samp]
# 	dbh.samp <- dbh[idx][samp]
# 	n <- length(dbh.srs) + length(dbh.samp)
# 	n.samp <- length(dbh.samp)
# 	z <- vector(length = n.samp)
# 	z[1] <- sum(dbh.srs) + dbh.samp[1] / prob[1]
# 	for(i in 2:n.samp){
# 		z[i] <- sum(dbh.srs) + sum(dbh.samp[1:(i - 1)]) + dbh.samp[i] / prob[i]
# 	}
#   1 / (N * n) * (N * n0 * mean(dbh.srs) + sum(z))
# }

# test <- sapply(1:100, make.two.stage.samp, dbh = dbh, bio = bio, n = n, n0 = 50, N0 = N-4, method = 'design')
# summary(test)
# mean(na.omit(test))
# dbh.mn
