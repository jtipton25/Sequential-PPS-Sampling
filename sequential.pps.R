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

make.pps.samp <- function(dbh, bio, n){ # sample dbh and bio using pps on dbh
	N <- length(dbh)
	p <- dbh / sum(dbh)
	samp <- sample(1:N, n, prob = p)
	prob <- p[samp]
	#dbh.pps.mn <- 1 / N * sum(1 / n * dbh[samp] / p[samp])
	dbh.pps <- dbh[samp]
	dbh.mn <- 1 / N * sum(dbh[samp] / (n * p[samp]))
	#delta.kl <- 1 - 
	#dbh.pps.var <- 
	#bio.pps.mn <- 1 / N * sum(1 / n * bio[samp] / p[samp])	
	bio.pps <- bio[samp]
	bio.mn <- 1 / N * sum(bio[samp] / (n * p[samp]))
	#bio.pps.var <- 
	list(dbh = dbh.pps, dbh.mn = dbh.mn, bio = bio.pps, bio.mn = bio.mn, p = prob)
	#list(dbh.pps = dbh, dbh.mn = dbh.mn, dbh.var = dbh.var, bio = bio, bio.mn = bio.mn, bio.var = bio.var)
}

make.bias.pps.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
	out <- make.pps.samp(dbh, bio, n)
	model <- lm(log(out$bio) ~ log(out$dbh))
	model.wt <- lm(log(out$bio) ~ log(out$dbh), weights = out$p)
	#list(coef = model$coef, coef.wt = model.wt$coef)
	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2])
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
			#bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2) #* (N - j)
			#bias[j] <- (factorial(N - 2) / factorial(N) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
			bias[j] <- (1 / (N * (N - 1)) * (sum(1:(N-1) * 1:(N-1) / (2:N))) + 1 /2 - (N - 1) / N) * (N - j)
		}
		bias[N] <- 0
		return(bias)
	}
}

make.design.samp <- function(dbh, bio, n){ # sample according to the proposed design
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
	pi <- (make.pi(N) + 1:N/N)[samples == 1]
	list(dbh = dbh.design, bio = bio.design, p = pi)
}

make.bias.design.est <- function(iter, dbh, bio, n){ # estimate bias in regression coefficients from pps sampling
	out <- make.design.samp(dbh, bio, n)
	model <- lm(log(out$bio) ~ log(out$dbh))
	model.wt <- lm(log(out$bio) ~ log(out$dbh), weights = out$p)
	#list(coef = model$coef, coef.wt = model.wt$coef)
	c(summary(model)$coef[, 1], summary(model)$coef[, 2], summary(model.wt)$coef[, 1], summary(model.wt)$coef[, 2])
}


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
## True PPS sampling design
##

out.pps <- make.pps.samp(dbh, bio, n)
dbh.mn - out.pps$dbh.mn
bio.mn - out.pps$bio.mn

##
## Plot relationship for PPS sample
##

make.model.plot(out.pps$dbh, out.pps$bio)

##
## Estimate Bias from pps sampling
##

bias.pps <- sapply(1:100, make.bias.pps.est, dbh = dbh, bio = bio, n = 100)
rownames(bias.pps) <- c('EST intercept OLS', 'EST slope OLS', 'SE intercept OLS', 'SE slope OLS', 'EST intercept WLS', 'EST slope WLS', 'SE intercept WLS', 'SE slope WLS')
idx.mn <- c(1:2, 5:6)
idx.var <- c(3:4, 7:8)
apply(bias.pps[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
apply(bias.pps[idx.mn, ], 1, var) * dim(bias.pps)[2] - apply(bias.pps[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters

##
## Sequential PPS design
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


out.design <- make.design.samp(dbh, bio, n)

make.model.plot(out.design$dbh, out.design$bio)
make.model.plot(dbh, bio)

bias.design <- sapply(1:1000, make.bias.design.est, dbh = dbh, bio = bio, n = 100)
rownames(bias.design) <- c('EST intercept OLS', 'EST slope OLS', 'SE intercept OLS', 'SE slope OLS', 'EST intercept WLS', 'EST slope WLS', 'SE intercept WLS', 'SE slope WLS')
idx.mn <- c(1:2, 5:6)
idx.var <- c(3:4, 7:8)
apply(bias.design[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
apply(bias.design[idx.mn, ], 1, var) * (dim(bias.design)[2] - 1) - apply(bias.design[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters



##
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
##


##
## Informative sampling test idea
##
model.no.informative <- lm(log(bio.pps) ~ log(dbh.pps))
model.informative <- lm(log(bio.pps) ~ log(dbh.pps) + 1 / fn(dbh.pps) + log(dbh.pps) * 1 / fn(dbh.pps))
anova(model.no.informative, model.informative)




#save.image(file = "pps.test.RData")

hist(pi.save - 2 * n / N * cdffn(dbh), freq = FALSE)
lines(density(pi.save - 2 * n / N * cdffn(dbh)))

