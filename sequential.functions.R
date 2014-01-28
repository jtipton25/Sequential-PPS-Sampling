##
## libraries and functions
##

<<<<<<< HEAD:sequential.functions.R
make.gamma.mixture <- function(N, alpha, beta){
	# makes a mixture of gamma distributions given alpha and beta vectors
	n <- length(alpha)
	samp <- sample(1:n, N, replace = TRUE)    
	dbh <- rgamma(N, alpha[samp], beta[samp])
}

make.sim.biomass <- function(dbh, b0, b1, s2){
	# make a power law distribution function with exponential random normal error and restricts values less than 0
	N <- length(dbh)
	bio <- b0 * dbh ^ b1 * exp(rnorm(N, 0, s2))
	bio[bio < 0] <- min(bio)
	if(length(bio[bio < 0]) == 0){
		return(bio)
	} else{
		return(bio)
		"some values of biomass were less than 0 so were set to min(bio)"
	}
}

make.model.plot <- function(dbh, bio, file = 'pathname'){ # plot model fits and save to file using file = 'pathname'
	if(file != 'pathname'){
		pdf(file = "fullModel")
	}
	N <- length(dbh)
	layout(matrix(1:4, nrow = 2))
	hist(dbh, breaks = floor(N / 10))
	hist(bio, breaks = floor(N / 10))
	# model
	model <- lm(log(bio) ~ log(dbh))
	newdbh <- seq(min(dbh), max(dbh), length.out = N)
	preds <- predict(model, newdata = data.frame(dbh = newdbh), interval = "predict")
	
	plot(bio ~ dbh)
	curve(exp(model$coeff[1]) * x^model$coeff[2], add = TRUE)
	polygon(c(rev(newdbh), newdbh), c(rev(exp(preds[, 3])), exp(preds[, 2])), col = adjustcolor('grey80', alpha.f=0.5), border = NA)
	lines(newdbh, exp(preds[ ,3]), lty = 'dashed', col = 'red')
	lines(newdbh, exp(preds[ ,2]), lty = 'dashed', col = 'red')
	
	plot(log(bio) ~ log(dbh), main = 'log scale')
	abline(model)
	polygon(c(rev(log(newdbh)), log(newdbh)), c(rev(preds[, 3]), preds[, 2]), col = adjustcolor('grey80', alpha.f=0.5), border = NA)
	lines(log(newdbh), preds[ ,3], lty = 'dashed', col = 'red')
	lines(log(newdbh), preds[ ,2], lty = 'dashed', col = 'red')
	if(file != 'pathname'){
		dev.off()
	}
}

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
			bias[j] <- (1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
		}
		bias[N] <- 0
		return(bias + 1:N / N)
	}
}

make.samp <- function(dbh, bio, n, method = 'srs'){
	N <- length(dbh)
	if(method == 'srs'){
		p <- n * rep(1 / N, N)
		samp <- sample(1:N, n, prob = p)
	}
	if(method == 'pps'){
		p <- n * dbh / sum(dbh)
		samp <- sample(1:N, n, prob = p)
	}
	if(method == 'ecdf'){
		p <- 2 * n / N * ecdf(dbh)(dbh)
		samp <- which(rbinom(1:N, 1, p) == 1)
	}
	if(method == 'design'){
		idx <- sample(1:N)
		p <- vector(length = N)
		for(k in 1:N){
			x <- dbh[idx[1:k]]
			fn <- ecdf(x)
			p[k] <- fn(x)[k]
		}
		p <- p * 2 * n / N # potential sample size adjustment
		samp <- which(rbinom(N, 1, p) == 1)
		idx <- order(dbh)
		#p <- make.pi(N)
	}
	if(method == 'strat'){
	  if(n %% 3 == 0){
		  n1 <- n / 3
  		n2 <- n / 3
	   	n3 <- n / 3
  	}
	  if(n %% 3 == 1){
		  n1 <- floor(n / 3)
  		n2 <- floor(n / 3)
	   	n3 <- ceiling(n / 3)
  	}
	  if(n %% 3 == 2){
		  n1 <- floor(n / 3)
		  n2 <- ceiling(n / 3)
		  n3 <- ceiling(n / 3)
  	}	
	  N <- length(dbh)
  	dbh.quant <- make.gamma.mixture(N, alpha = c(2, 4, 6, 8), beta = c(12, 10, 8, 6))
	  quant <- quantile(dbh.quant, probs = c(1/3, 2/3))
  	q.idx.1 <- which(dbh < quant[1])
  	q.idx.2 <- which(quant[1] <= dbh & dbh < quant[2])
  	q.idx.3 <- which(quant[2] <= dbh)
  	samp1 <- sample(q.idx.1, n1)
  	samp2 <- sample(q.idx.2, n2)
  	samp3 <- sample(q.idx.3, n3)
  	samp <- c(samp1, samp2, samp3)
  	dbh.samp1 <- dbh[samp1]
  	dbh.samp2 <- dbh[samp2]
  	dbh.samp3 <- dbh[samp3]
  	bio.samp1 <- bio[samp1]
  	bio.samp2 <- bio[samp2]
  	bio.samp3 <- bio[samp3]
  	nh <- c(length(samp1), length(samp2), length(samp3))
  	Nh <- c(length(q.idx.1), length(q.idx.2), length(q.idx.3))
  	muh.dbh <- c(mean(dbh.samp1), mean(dbh.samp2), mean(dbh.samp3))
  	sh.dbh <- c(var(dbh.samp1), var(dbh.samp2), var(dbh.samp3))
  	muh.bio <- c(mean(bio.samp1), mean(bio.samp2), mean(bio.samp3))
  	dbh.mn <- 1 / N * sum(Nh * muh.dbh)
  	bio.mn <- 1 / N * sum(Nh * muh.bio)
  	dbh.samp <- dbh[samp]
  	bio.samp <- bio[samp]
  	list(dbh = dbh.samp, dbh.mn = dbh.mn, bio = bio.samp, bio.mn = bio.mn, quant = quant, samp = samp)
	} else {
	prob <- p[samp]
	dbh.samp <- dbh[samp]
	dbh.mn <- 1 / N * sum(dbh.samp / prob)
	bio.samp <- bio[samp]
	bio.mn <- 1 / N * sum(bio.samp / prob)	
	list(dbh = dbh.samp, dbh.mn = dbh.mn, bio = bio.samp, bio.mn = bio.mn, p = p, samp = samp)
	}
}

make.bias.est <- function(iter, dbh, bio, n, method = 'srs'){ # estimate bias in regression coefficients from sampling design
	if(method == 'strat'){
		out <- make.samp.strat(dbh, bio, n)	
	}	else {
		out <- make.samp(dbh, bio, n, method)
	}
	model <- lm(log(bio) ~ log(dbh), data = data.frame(dbh = out$dbh, bio = out$bio))
	model.wt <- lm(log(bio) ~ log(dbh), weights = out$p[out$samp], data = data.frame(dbh = out$dbh, bio = out$bio))
	newdata <- data.frame(dbh = dbh[ - out$samp])
	preds <- predict(model, newdata = newdata)
	preds.wt <- predict(model.wt, newdata = newdata)
	pred.mse <- mean((exp(preds) - bio[ - out$samp])^2)
	pred.wt.mse <- mean((exp(preds.wt) - bio[ - out$samp])^2)
	#c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2], pred.mse, pred.wt.mse)
	bias <- c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2], pred.mse, pred.wt.mse)
	names(bias) <- c('EST intercept OLS', 'EST slope OLS', 'SE intercept OLS', 'SE slope OLS', 'EST intercept WLS', 'EST slope WLS', 'SE intercept WLS', 'SE slope WLS', 'MSPE OLS', 'MSPE WLS')
	return(bias)
}

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
	if(method == 'strat'){
		return(out$bio.mn)
	}
}
=======
source("sequential.functions.R")

library(parallel)
library(snowfall)
library(rlecuyer)
cps <- detectCores()
sfInit(parallel = TRUE, cpus = cps)
sfClusterSetupRNG() 

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



sfExportAll()

temp.srs <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = n, method = 'srs') 
summary(temp.srs)
mean(temp.srs)
var(temp.srs)

temp.pps <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = n, method = 'pps') 
summary(temp.pps)
mean(temp.pps)
var(temp.pps)

temp.ecdf <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = n, method = 'ecdf') 
summary(temp.ecdf)
mean(temp.ecdf)
var(temp.ecdf)

temp.design <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = n, method = 'design') 
summary(temp.design)
mean(temp.design)
var(temp.design)

temp.strat <- sfSapply(1:1000, make.est.strat.result, dbh = dbh, bio = bio, n = n) 
summary(temp.strat)
mean(temp.strat)
var(temp.strat)


bio.mn
mean(temp.srs)
mean(temp.ecdf)
mean(temp.pps)
mean(temp.design)
mean(temp.strat)

var(temp.srs)
var(temp.ecdf)
var(temp.pps)
var(temp.design)
var(temp.strat)

out <- make.samp.strat(dbh, bio, n)

make.model.plot(out$dbh, out$bio)

dbh.mn

>>>>>>> f2051084a56dd102b5547ab2037cb6c20fb222d4:stratified.R
