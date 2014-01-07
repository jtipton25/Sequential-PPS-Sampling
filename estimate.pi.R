library(permute)
library(MASS)
ctrl <- how(maxperm = 10000000)
N <- 8
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
as.fractions(apply(pi, 2, mean))
## N = 5, denominator of 1200
137*4
237*3
437*2

137*4 - 240
237*3 - 240*2
437*2 - 240*3
1037 - 240*4

## N = 6 denominator of 600
49*5
79*4
129*3
229*2

49*5 - 100
79*4 - 200
129*3 - 300
229*2 - 400
529 - 500

## N = 7 denominator of 3039960
363*6*517
559*5*517
## 517 has divisors of 11 and 47
## 5880
300*5880
1343*3*517
2323*2*517
5263*517

363*6*517 - 840*517
559*5*517 - 840*2*517
300*5880 - 840*3*517
1343*3*517 - 840*4*517
2323*2*517 - 840*5*517
5263*517 - 840*6*517
115

## N=8 denominator of 2*2*2*2*2*2*5*7*7 * 5*17*19 * 3*3*3*3*41 * 11*457 
1615
## 1615 = 5*17*19
15680/5027
3321
## 3321 = 3*3*3*3*41
5027/11
## 5027 = 11*457 
15680/(2^6)/5
##15680 = 2*2*2*2*2*2*5*7*7
as.fractions(apply(pi, 2, mean))

