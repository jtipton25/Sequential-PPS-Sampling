##
## Sequential PPS sampling design simulation
##

setwd('~/Google Drive/PalEON/PalEON Meeting Data/Forest ECDF Sampling/')

##
## Simulate dbh
##

N <- 1000
n <- 100
#alpha <- 2
#beta <- 6
#curve(dgamma(x, alpha, beta))
#dbh <- rgamma(N, alpha, beta)
#dbh <- runif(N)
alpha <- c(2, 4, 6, 8) 
beta <- c(12, 10, 8, 6)
layout(matrix(1:4, 2, 2))
for(i in 1:4){
	curve(dgamma(x, alpha[i], beta[i]), from = 0, to = 4)
}
samp.density <- sample(1:4, N, replace = TRUE)
dbh <- 5 + 10 * rgamma(N, alpha[samp.density], beta[samp.density])
par(mfrow = c(1, 1))
#pdf(file = 'dbh1.pdf')
hist(dbh, breaks = 20)
#dev.off()
truth.dbh <- 1 / N * sum(dbh)
truth.dbh

##
## Simulate Biomass
##

a <- 2
b <- 2
s2 <- 1 / 4
epsilon <- rnorm(N, 0, s2)
bio <- a * dbh ^ b * exp(epsilon)
bio[bio < 0]
plot(log(bio) ~ log(dbh))
model <- lm(log(bio) ~ log(dbh))
abline(model)
plot(resid(model) ~ fitted(model))
abline(h=0)
plot(bio ~ dbh)
curve(exp(model$coeff[1]) * x^model$coeff[2], add = TRUE)

truth.bio <- mean(bio)

## Plot of data
#pdf(file = 'dbhLogModel.pdf')
plot(log(bio) ~ log(dbh))
model <- lm(log(bio) ~ log(dbh))
abline(model)
newdbh <- seq(min(dbh), max(dbh), length.out = 1000)
preds <- predict(model, newdata = data.frame(dbh = newdbh), interval = "confidence")
polygon(c(rev(log(newdbh)), log(newdbh)), c(rev(preds[, 3]), preds[, 2]), col = 'grey80', border = NA)
lines(log(newdbh), preds[ ,3], lty = 'dashed', col = 'red')
lines(log(newdbh), preds[ ,2], lty = 'dashed', col = 'red')
#dev.off()
#pdf(file = 'dbhModelResid.pdf')
plot(resid(model) ~ fitted(model))
abline(h=0)
#dev.off()
#pdf(file = 'dbhModel.pdf')
plot(bio ~ dbh)
curve(exp(model$coeff[1]) * x^model$coeff[2], add = TRUE)
polygon(c(rev(newdbh), newdbh), c(rev(exp(preds[, 3])), exp(preds[, 2])), col = 'grey80', border = NA)
lines(newdbh, exp(preds[ ,3]), lty = 'dashed', col = 'red')
lines(newdbh, exp(preds[ ,2]), lty = 'dashed', col = 'red')
#dev.off()

##
## True PPS sampling design
##

p <- dbh / sum(dbh)
hist(p)

samp <- sample(1:N, n, prob = p)
estimate.dbh <- 1 / N * sum(1 / n * dbh[samp] / p[samp])
estimate.bio <- 1 / N * sum(1 / n * bio[samp] / p[samp])

truth.dbh - estimate.dbh
truth.bio - estimate.bio

##
## Compare True allometric relationship to pps allometric relationship
##

model.pps <- lm(log(bio[samp]) ~ log(dbh[samp]))
model
model.pps

layout(matrix(1:4, 2))
plot(log(bio) ~ log(dbh))
abline(model)
plot(bio ~ dbh)
curve(exp(model$coeff[1]) * x^model$coeff[2], add = TRUE)
plot(log(bio[samp]) ~ log(dbh[samp]))
abline(model.pps)
plot(bio[samp] ~ dbh[samp])
curve(exp(model.pps$coeff[1]) * x^model.pps$coeff[2], add = TRUE)


##
## Estimate Bias from pps sampling
##
iter <- 100
coefs <- matrix(nrow = iter, ncol = 2)
for(j in 1:iter){
	samp <-  sample(1:N, n, prob = p)
	model.pps <- lm(log(bio[samp]) ~ log(dbh[samp]))
	coefs[j, ] <- model.pps$coef
}

apply(coefs, 2, mean)
c(log(a), b) 
## Seems to be unbiasedly estimating the regression parameters

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

#dbh <- runif(N)
#alpha <- 2
#beta <- 6

hist(dbh, freq = FALSE, breaks = 20)
#curve(dgamma(x, alpha, beta))

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

