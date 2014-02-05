##
## libraries and functions
##

source("sequential.functions.R")

library(xtable)
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
pdf(file = 'dbhModel.pdf')
make.model.plot(dbh, bio)
dev.off()

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

temp.strat <- sfSapply(1:1000, make.est.result, dbh = dbh, bio = bio, n = n, method = 'strat') 
summary(temp.strat)
mean(temp.strat)
var(temp.strat)


tab <- cbind(c(mean(temp.srs) - bio.mn, mean(temp.ecdf) - bio.mn, mean(temp.pps) - bio.mn, mean(temp.design) - bio.mn, mean(temp.strat) - bio.mn), c(var(temp.srs) / var(temp.srs), var(temp.ecdf) / var(temp.srs), var(temp.pps) / var(temp.srs), var(temp.design) / var(temp.srs), var(temp.strat) / var(temp.srs)))
rownames(tab) = c('SRS', 'ECDF', 'PPS', "AECDF", "STSI")
colnames(tab) <- c('Bias', 'Relative Efficiency')
xtable(tab)


out <- make.samp(dbh, bio, n, method = 'strat')

make.model.plot(out$dbh, out$bio)

dbh.mn

