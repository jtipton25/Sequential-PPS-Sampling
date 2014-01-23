##
## Sequential PPS sampling design simulation
##

##
## libraries and functions
##


source('sequential.functions.R')

##
## Simulate dbh
##

N <- 300 # finite population size
n <- 30 # expected sample size

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
## SRS sampling
##

out.srs <- make.samp(dbh, bio, n, method = 'srs')
dbh.mn - out.srs$dbh.mn
bio.mn - out.srs$bio.mn

##
## Plot relationship for SRS sample
##

make.model.plot(out.srs$dbh, out.srs$bio)

##
## Estimate Bias from SRS sampling
##

bias.srs <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = n, method = 'srs')
idx.mn <- c(1:2, 5:6)
idx.var <- c(3:4, 7:8)
idx.mspe <- 9:10
apply(bias.srs[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.srs[idx.mn, ], 1, var) - apply(bias.srs[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.srs[idx.mspe, ], 1, mean) 

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

bias.pps <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = n, method = 'pps')
apply(bias.pps[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.pps[idx.mn, ], 1, var) - apply(bias.pps[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.pps[idx.mspe, ], 1, mean) 

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

bias.ecdf <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = n, method = 'ecdf')
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

#out.design <- make.design.samp(dbh, bio, n)
out.design <- make.samp(dbh, bio, n, method = "design")
make.model.plot(out.design$dbh, out.design$bio)
make.model.plot(dbh, bio)

bias.design <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = n, method = 'design')
apply(bias.design[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.design[idx.mn, ], 1, var) - apply(bias.design[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.design[idx.mspe, ], 1, mean) # seems to be similar to pps (at least for this particular finite population)



##
## Estimate Bias from stratified sampling
##

bias.strat <- sapply(1:1000, make.bias.est, dbh = dbh, bio = bio, n = n, method = 'strat')
apply(bias.strat[idx.mn,], 1, mean) - rep(c(log(b0), b1), 2) # seems to be unbiasedly estimating the regression parameters
## variance of the means - mean of the variances
apply(bias.strat[idx.mn, ], 1, var) - apply(bias.design[idx.var, ], 1, mean) # variance of estimator vs estimated variance for regression parameters seems to be unbiased
# mean MSPE
apply(bias.strat[idx.mspe, ], 1, mean) # seems to be similar to pps (at least for this particular finite population)

##
## Compare bias estimates
##

apply(bias.srs[idx.mspe, ], 1, mean)
apply(bias.pps[idx.mspe, ], 1, mean)
apply(bias.ecdf[idx.mspe, ], 1, mean)
apply(bias.design[idx.mspe, ], 1, mean)
apply(bias.strat[idx.mspe, ], 1, mean)





