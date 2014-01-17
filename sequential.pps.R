##
## Sequential PPS sampling design simulation
##

##
## libraries and functions
##

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

<<<<<<< HEAD
=======
make.pps.samp <- function(dbh, bio, n){ # sample dbh and bio using pps on dbh
	N <- length(dbh)
	p <- dbh / sum(dbh)
	samp <- sample(1:N, n, prob = p)
	prob <- p[samp]
	#dbh.pps.mn <- 1 / N * sum(1 / n * dbh[samp] / p[samp])
	dbh.pps <- dbh[samp]
	dbh.mn <- 1 / N * sum(dbh.pps / (n * prob))
	dbh.mn
	#delta.kl <- 1 - 
	#dbh.pps.var <- 
	#bio.pps.mn <- 1 / N * sum(1 / n * bio[samp] / p[samp])	
	bio.pps <- bio[samp]
	bio.mn <- 1 / N * sum(bio.pps / (n * prob))
	#bio.pps.var <- 
	list(dbh = dbh.pps, dbh.mn = dbh.mn, bio = bio.pps, bio.mn = bio.mn, p = p, samp = samp)
	#list(dbh.pps = dbh, dbh.mn = dbh.mn, dbh.var = dbh.var, bio = bio, bio.mn = bio.mn, bio.var = bio.var)
}

## add in a MSPE criteria
make.bias.pps.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
	out <- make.pps.samp(dbh, bio, n)
	model <- lm(log(out$bio) ~ log(out$dbh))
	model.wt <- lm(log(out$bio) ~ log(out$dbh), weights = out$p[out$samp])
	#list(coef = model$coef, coef.wt = model.wt$coef)
	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2])
}

>>>>>>> 6823e974f89a36aa39a1342b4a6b6f1caf759ab0
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

<<<<<<< HEAD
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
	p <- make.pi(N)
	}  
	prob <- p[samp]
	dbh.samp <- dbh[samp]
	dbh.mn <- 1 / N * sum(dbh.samp / prob)
	bio.samp <- bio[samp]
  bio.mn <- 1 / N * sum(bio.samp / prob)	
	list(dbh = dbh.samp, dbh.mn = dbh.mn, bio = bio.samp, bio.mn = bio.mn, p = p, samp = samp)
}

make.bias.est <- function(iter, dbh, bio, n, method = 'srs'){ # estimate bias in regression coefficients from sampling design
	out <- make.samp(dbh, bio, n, method)
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
=======
make.ecdf.samp <- function(dbh, bio, n){ # sample dbh and bio using pps on dbh
	N <- length(dbh)
	p <- 2 * n / N * ecdf(dbh)(dbh)
	samples <- rbinom(1:N, 1, p)
	samp <- which(samples == 1)
	prob <- p[samp]
	dbh.ecdf <- dbh[samp]
	dbh.mn <- 1 / N * sum(dbh.ecdf / prob)
	#delta.kl <- 1 - 
	#dbh.pps.var <- 
	#bio.pps.mn <- 1 / N * sum(1 / n * bio[samp] / p[samp])	
	bio.ecdf <- bio[samp]
	bio.mn <- 1 / N * sum(bio.ecdf / prob)
	#bio.pps.var <- 
	list(dbh = dbh.ecdf, dbh.mn = dbh.mn, bio = bio.ecdf, bio.mn = bio.mn, p = p, samp = samp)
	#list(dbh.pps = dbh, dbh.mn = dbh.mn, dbh.var = dbh.var, bio = bio, bio.mn = bio.mn, bio.var = bio.var)
}

## add in a MSPE criteria
make.bias.ecdf.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
	out <- make.pps.samp(dbh, bio, n)
	model <- lm(log(out$bio) ~ log(out$dbh))
	model.wt <- lm(log(out$bio) ~ log(out$dbh), weights = out$p[out$samp])
	#list(coef = model$coef, coef.wt = model.wt$coef)
	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2])
}

make.design.samp <- function(dbh, bio, n){ # sample according to the proposed design
	N <- length(dbh)
	p <- vector(length = N)
	for(k in 1:N){
		x <- dbh[1:k]
		fn <- ecdf(x)
		p[k] <- fn(x)[k]
	}
	p <- p * 2 * n / N # potential sample size adjustment
	samp <- which(rbinom(N, 1, p) == 1)
	pi <- (make.pi(N) + 1:N/N)#[samples == 1]
	pi.samp <- pi[samp]
	dbh.design <- dbh[samp]
	dbh.mn <- 1 / N * sum(dbh.design / pi.samp)
	bio.design <- bio[samp]
	bio.mn <- 1 / N * sum(bio.design / pi.samp)
		list(dbh = dbh.design, dbh.mn = dbh.mn, bio = bio.design, bio.mn = bio.mn, p = pi, samp = samp)
}
## add in a MSPE criteria
make.bias.design.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
	out <- make.design.samp(dbh, bio, n)
	model <- lm(log(out$bio) ~ log(out$dbh))
	model.wt <- lm(log(out$bio) ~ log(out$dbh), weights = out$p[out$samp])
	#list(coef = model$coef, coef.wt = model.wt$coef)
	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2])
}

>>>>>>> 6823e974f89a36aa39a1342b4a6b6f1caf759ab0

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
<<<<<<< HEAD

bio <- make.sim.biomass(dbh, b0, b1, s2)
bio.mn <- mean(bio)
bio.var <- var(bio)

##
## Plot Data
##

#make.model.plot(dbh, bio, file = "fullModel.pdf")
make.model.plot(dbh, bio)
=======
>>>>>>> 6823e974f89a36aa39a1342b4a6b6f1caf759ab0

##
## SRS sampling
##

out.srs <- make.samp(dbh, bio, n, method = 'srs')
dbh.mn - out.srs$dbh.mn
bio.mn - out.srs$bio.mn

##
## Plot relationship for SRS sample
##

<<<<<<< HEAD
make.model.plot(out.srs$dbh, out.srs$bio)

##
## Estimate Bias from SRS sampling
##

bias.srs <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = 100, method = 'srs')
idx.mn <- c(1:2, 5:6)
idx.var <- c(3:4, 7:8)
idx.mspe <- 9:10
apply(bias.srs[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.srs[idx.mn, ], 1, var) - apply(bias.srs[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.srs[idx.mspe, ], 1, mean) 
=======
#make.model.plot(dbh, bio, file = "fullModel.pdf")
make.model.plot(dbh, bio)
>>>>>>> 6823e974f89a36aa39a1342b4a6b6f1caf759ab0

##
## True PPS sampling design
##

#out.pps <- make.pps.samp(dbh, bio, n)
out.pps <- make.samp(dbh, bio, n, method = 'pps')
dbh.mn - out.pps$dbh.mn
bio.mn - out.pps$bio.mn

##
## Plot relationship for PPS sample
##

make.model.plot(out.pps$dbh, out.pps$bio)

##
## Estimate Bias from pps sampling
##

<<<<<<< HEAD
bias.pps <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = 100, method = 'pps')
apply(bias.pps[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.pps[idx.mn, ], 1, var) - apply(bias.pps[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.pps[idx.mspe, ], 1, mean) 
=======
bias.pps <- sapply(1:1000, make.bias.pps.est, dbh = dbh, bio = bio, n = 100)
rownames(bias.pps) <- c('EST intercept OLS', 'EST slope OLS', 'SE intercept OLS', 'SE slope OLS', 'EST intercept WLS', 'EST slope WLS', 'SE intercept WLS', 'SE slope WLS')
idx.mn <- c(1:2, 5:6)
idx.var <- c(3:4, 7:8)
apply(bias.pps[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.pps[idx.mn, ], 1, var) - apply(bias.pps[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased

##
## ecdf sampling design
##

out.ecdf <- make.ecdf.samp(dbh, bio, n)
dbh.mn - out.ecdf$dbh.mn
bio.mn - out.ecdf$bio.mn

##
## Plot relationship for ecdf sample
##

make.model.plot(out.ecdf$dbh, out.ecdf$bio)

##
## Estimate Bias from ecdf sampling
##

bias.ecdf <- sapply(1:1000, make.bias.ecdf.est, dbh = dbh, bio = bio, n = 100)
rownames(bias.ecdf) <- c('EST intercept OLS', 'EST slope OLS', 'SE intercept OLS', 'SE slope OLS', 'EST intercept WLS', 'EST slope WLS', 'SE intercept WLS', 'SE slope WLS')
idx.mn <- c(1:2, 5:6)
idx.var <- c(3:4, 7:8)
apply(bias.ecdf[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.ecdf[idx.mn, ], 1, var) - apply(bias.ecdf[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
>>>>>>> 6823e974f89a36aa39a1342b4a6b6f1caf759ab0

##
## ecdf sampling design
##

#out.ecdf <- make.ecdf.samp(dbh, bio, n)
out.ecdf <- make.samp(dbh, bio, n, method = "ecdf")
dbh.mn - out.ecdf$dbh.mn
bio.mn - out.ecdf$bio.mn

##
## Plot relationship for ecdf sample
##

make.model.plot(out.ecdf$dbh, out.ecdf$bio)

##
## Estimate Bias from ecdf sampling
##

bias.ecdf <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = 100, method = 'ecdf')
apply(bias.ecdf[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.ecdf[idx.mn, ], 1, var) - apply(bias.ecdf[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.ecdf[idx.mspe, ], 1, mean) # seems to be similar to pps (at least for this particular finite population)

##
## Sequential ecdf design
##

##
## Sample the first element (x1) with probability 1
## Sample the second element (x2) with probability 1 if x2 > x1
##                                                 1/2 if x2 <= x1  
##
## Sample the third element (x3) with probability 1 if x3 > x1 & x2
##                                                2/3 if x3 > x1 & x3 < x2 
##                                                2/3 if x3 > x2 & x3 < x1
##                                                1/3 if x3 < x1 & x2
##
## And so on
##

##
## Testing the behavior of the probability of being sampled pi_i for one population
##

<<<<<<< HEAD
#out.design <- make.design.samp(dbh, bio, n)
out.design <- make.samp(dbh, bio, n, method = "design")
make.model.plot(out.design$dbh, out.design$bio)
make.model.plot(dbh, bio)

bias.design <- sapply(1:100, make.bias.est, dbh = dbh, bio = bio, n = 100, method = 'design')
apply(bias.design[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.design[idx.mn, ], 1, var) - apply(bias.design[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.design[idx.mspe, ], 1, mean) # seems to be similar to pps (at least for this particular finite population)
=======
out.design <- make.design.samp(dbh, bio, n)

make.model.plot(out.design$dbh, out.design$bio)
make.model.plot(dbh, bio)

bias.design <- sapply(1:100, make.bias.design.est, dbh = dbh, bio = bio, n = 100)
rownames(bias.design) <- c('EST intercept OLS', 'EST slope OLS', 'SE intercept OLS', 'SE slope OLS', 'EST intercept WLS', 'EST slope WLS', 'SE intercept WLS', 'SE slope WLS')
idx.mn <- c(1:2, 5:6)
idx.var <- c(3:4, 7:8)
apply(bias.design[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.design[idx.mn, ], 1, var) - apply(bias.design[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters - seems to be unbiased








##
## Older code...
##


make.design.samp <- function(dbh, bio, n){
	N <- length(dbh)
	p <- vector(length = N)
	for(k in 1:N){
		x <- dbh[1:k]
		fn <- ecdf(x)
		p[k] <- fn(x)[k]
	}
	#	p <- p * 2 * n / N # potential sample size adjustment
	samples <- rbinom(N, 1, p)
	dbh.design <- dbh[samples == 1]
	bio.design <- bio[samples == 1]
	list(dbh = dbh.design, bio = bio.design)
}

out.design <- make.design.samp(dbh, bio, n)

make.model.plot(out.design$dbh, out.design$bio)
make.model.plot(dbh, bio)
>>>>>>> 6823e974f89a36aa39a1342b4a6b6f1caf759ab0

apply(bias.srs[idx.mspe, ], 1, mean)
apply(bias.pps[idx.mspe, ], 1, mean)
apply(bias.ecdf[idx.mspe, ], 1, mean)
apply(bias.design[idx.mspe, ], 1, mean)



##
<<<<<<< HEAD
## Older code
=======
## Older code...
##



iter <- 10000
est.mn.dbh <- vector(length = iter)
est.var.dbh <- vector(length = iter)
est.mn.bio <- vector(length = iter)
est.var.bio <- vector(length = iter)
bio.fit <- vector(mode = 'list', length = iter)
bio.fit.weight <- vector(mode = 'list', length = iter)
bio.fit.srs <- vector(mode = 'list', length = iter)
pi.save <- rep(0, length = N)
n.save <- vector(length = iter)
s2.core <- 1/8 # Sampling error from coring tree
coefs <- matrix(nrow = iter, ncol = 2)

for(i in 1:iter){
	if(i %% 100 == 0){
		cat(i, ' ')
	}
	p <- vector(length = N)
	samp.srs <- sample(1:N, n) #SRS sampling for comparison
	samp <- sample(1:N)
	dbh.samp <- dbh[samp]
	bio.samp <- bio[samp]
	for(k in 1:N){
		x <- dbh.samp[1:k]
		fn <- ecdf(x)
		p[k] <- fn(x)[k]
	}
#	p <- p * 2 * n / N #potential sample size adjustment
	samples <- rbinom(N, 1, p)
	dbh.pps <- dbh.samp[which(samples == 1)]
	error <- rnorm(N, 0, s2.core)
	bio.pps <- (bio.samp*exp(error))[which(samples == 1)] #exponential error due to measurement error sampling the cores, not sampling design
	
	n.save[i] <- sum(samples)
	probs <- p[which(samples == 1)]
	p[samp] <- p
	pi.save <- pi.save + 1 / iter * p
	est.mn.dbh[i] <- 1 / N * sum(dbh.pps / probs) # check these
	est.var.dbh[i] <- (1 / N) * var(dbh.pps / probs)
	est.mn.bio[i] <- 1 / N * sum(bio.pps / probs) # check these
	est.var.bio[i] <- (1 / N) * var(bio.pps / probs)
	model.pps <- lm(log(bio.pps) ~ log(dbh.pps))
	coefs[i, ] <- model.pps$coef
	bio.fit[[i]] <- exp(model.pps$coef[1]) * dbh^model.pps$coef[2]  
	model.srs <- lm(log((bio*exp(error))[samp.srs]) ~ log(dbh[samp.srs]))
	bio.fit.srs[[i]] <- exp(model.srs$coef[1]) * dbh^model.srs$coef[2]  
	model.pps.weight <- lm(log(bio.pps) ~ log(dbh.pps), weights = 1 / fn(dbh.pps))
	bio.fit.weight[[i]] <- exp(model.pps.weight$coef[1]) * dbh^model.pps.weight$coef[2]  
}

#load("pps.test.RData")

apply(coefs, 2, mean) ## Looks pretty unbiased to me
c(log(a), b) 


mean(est.mn.dbh)
var(est.mn.dbh)
mean(est.var.dbh)
mean(dbh)
var(dbh) / N^2
mean(n.save)


## Plots suggest CLT behavior
layout(matrix(c(1:4), nrow = 2))
hist(est.mn.dbh, freq = FALSE, breaks = 20)
curve(dnorm(x, mean(est.mn.dbh), sqrt(N*mean(est.var.dbh))), add = TRUE)  
## Density function for dbh
#curve((1/4 * dbeta(x, alpha[1], beta[1]) + 1/4 * dbeta(x, alpha[2], beta[2]) + 1/4 * dbeta(x, alpha[3], beta[3]) + 1/4 * dbeta(x, alpha[4], beta[4])), ylab = "")  
hist(est.mn.bio, freq = FALSE, breaks = 20)
curve(dnorm(x, mean(est.mn.bio), sqrt(N*var(est.var.bio))), add = TRUE)
hist(n.save) # mean is about 1/2 number of trees N for this distribution
cdffn <- ecdf(dbh)
#hist(pi.save - 2 * n / N * cdffn(dbh) - 1/N, breaks = 20) # Not the same as pps sampling # controlling sample size
#mean(pi.save - 2 * n / N * cdffn(dbh)- 1/N) # controlling sample size

hist(pi.save - cdffn(dbh)); abline(v = 0, lwd = 8)
mean(pi.save - cdffn(dbh))
MSPE.srs <- vector(length = iter)
MSPE.pps <- vector(length = iter)
MSPE.pps.weight <- vector(length = iter)
for(i in 1:iter){
	MSPE.srs[i] <- mean((bio.fit.srs[[i]] - bio)^2)
	MSPE.pps[i] <- mean((bio.fit[[i]] - bio)^2)
	MSPE.pps.weight[i] <- mean((bio.fit.weight[[i]] - bio)^2)
}

mean(MSPE.srs)
mean(MSPE.pps)
mean(MSPE.pps.weight)


##
##
>>>>>>> 6823e974f89a36aa39a1342b4a6b6f1caf759ab0
##

#make.pps.samp <- function(dbh, bio, n){ # sample dbh and bio using pps on dbh
#	N <- length(dbh)
#	p <- dbh / sum(dbh)
#	samp <- sample(1:N, n, prob = p)
#	prob <- p[samp]
	#dbh.pps.mn <- 1 / N * sum(1 / n * dbh[samp] / p[samp])
#	dbh.pps <- dbh[samp]
#	dbh.mn <- 1 / N * sum(dbh.pps / (n * prob))
#	dbh.mn
	#delta.kl <- 1 - 
	#dbh.pps.var <- 
	#bio.pps.mn <- 1 / N * sum(1 / n * bio[samp] / p[samp])	
#	bio.pps <- bio[samp]
#	bio.mn <- 1 / N * sum(bio.pps / (n * prob))
	#bio.pps.var <- 
#	list(dbh = dbh.pps, dbh.mn = dbh.mn, bio = bio.pps, bio.mn = bio.mn, p = p, samp = samp)
	#list(dbh.pps = dbh, dbh.mn = dbh.mn, dbh.var = dbh.var, bio = bio, bio.mn = bio.mn, bio.var = bio.var)
#}

# add MSPE for total biomass
#make.bias.pps.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
#	out <- make.pps.samp(dbh, bio, n)
#	model <- lm(log(bio) ~ log(dbh), data = out)
#	model.wt <- lm(log(bio) ~ log(dbh), weights = out$p[out$samp], data = out)
#	newdata <- data.frame(dbh = dbh[ - out$samp])
#	preds <- predict(model, newdata = newdata)
#	preds.wt <- predict(model.wt, newdata = newdata)
#	pred.mse <- mean((exp(preds) - bio[ - out$samp])^2)
#	pred.wt.mse <- mean((exp(preds.wt) - bio[ - out$samp])^2)
#	#list(coef = model$coef, coef.wt = model.wt$coef, pred.mse = pred.mse, pred.wt.mse = pred.wt.mse)
#	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2], pred.mse, pred.wt.mse)
#}


#make.ecdf.samp <- function(dbh, bio, n){ # sample dbh and bio using pps on dbh
#	N <- length(dbh)
#	p <- 2 * n / N * ecdf(dbh)(dbh)
#	samples <- rbinom(1:N, 1, p)
#	samp <- which(samples == 1)
#	prob <- p[samp]
#	dbh.ecdf <- dbh[samp]
#	dbh.mn <- 1 / N * sum(dbh.ecdf / prob)
	#delta.kl <- 1 - 
	#dbh.pps.var <- 
	#bio.pps.mn <- 1 / N * sum(1 / n * bio[samp] / p[samp])	
#	bio.ecdf <- bio[samp]
#	bio.mn <- 1 / N * sum(bio.ecdf / prob)
	#bio.pps.var <- 
#	list(dbh = dbh.ecdf, dbh.mn = dbh.mn, bio = bio.ecdf, bio.mn = bio.mn, p = p, samp = samp)
#}

# add MSPE for total biomass
#make.bias.ecdf.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
#	out <- make.pps.samp(dbh, bio, n)
#	model <- lm(log(bio) ~ log(dbh), data = out)
#	model.wt <- lm(log(bio) ~ log(dbh), weights = out$p[out$samp], data = out)
#	newdata <- data.frame(dbh = dbh[ - out$samp])
#	preds <- predict(model, newdata = newdata)
#	preds.wt <- predict(model.wt, newdata = newdata)
#	pred.mse <- mean((exp(preds) - bio[ - out$samp])^2)
#	pred.wt.mse <- mean((exp(preds.wt) - bio[ - out$samp])^2)
	#list(coef = model$coef, coef.wt = model.wt$coef, pred.mse = pred.mse, pred.wt.mse = pred.wt.mse)
	#list(coef = model$coef, coef.wt = model.wt$coef)
#	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2], pred.mse, pred.wt.mse)
#}

#make.design.samp <- function(dbh, bio, n){ # sample according to the proposed design
#	N <- length(dbh)
#	p <- vector(length = N)
#	for(k in 1:N){
#		x <- dbh[1:k]
#		fn <- ecdf(x)
#		p[k] <- fn(x)[k]
#	}
#	p <- p * 2 * n / N # potential sample size adjustment
#	samp <- which(rbinom(N, 1, p) == 1)
#	pi <- (make.pi(N) + 1:N/N)#[samples == 1]
#	pi.samp <- pi[samp]
#	dbh.design <- dbh[samp]
#	dbh.mn <- 1 / N * sum(dbh.design / pi.samp)
#	bio.design <- bio[samp]
#	bio.mn <- 1 / N * sum(bio.design / pi.samp)
#	list(dbh = dbh.design, dbh.mn = dbh.mn, bio = bio.design, bio.mn = bio.mn, p = pi, samp = samp)
#}

# add MSPE for total biomass
#make.bias.design.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
#	out <- make.design.samp(dbh, bio, n)
#	model <- lm(log(bio) ~ log(dbh), data = data.frame(dbh = out$dbh, bio = out$bio))
#	model.wt <- lm(log(bio) ~ log(dbh), weights = out$p[out$samp], data = data.frame(dbh = out$dbh, bio = out$bio))
#	newdata <- data.frame(dbh = dbh[ - out$samp])
#	preds <- predict(model, newdata = newdata)
#	preds.wt <- predict(model.wt, newdata = newdata)
#	pred.mse <- mean((exp(preds) - bio[ - out$samp])^2)
#	pred.wt.mse <- mean((exp(preds.wt) - bio[ - out$samp])^2)
#	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2], pred.mse, pred.wt.mse)
	#list(coef = model$coef, coef.wt = model.wt$coef, pred.mse = pred.mse, pred.wt.mse = pred.wt.mse)
	#list(coef = model$coef, coef.wt = model.wt$coef)
	#c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2], pred.mse, pred.wt.mse)
#}

##
## Informative sampling test idea
##
model.no.informative <- lm(log(bio.pps) ~ log(dbh.pps))
model.informative <- lm(log(bio.pps) ~ log(dbh.pps) + 1 / fn(dbh.pps) + log(dbh.pps) * 1 / fn(dbh.pps))
anova(model.no.informative, model.informative)
