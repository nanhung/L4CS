library(randtoolbox)



  
u <- sobol(128, 7)
X <- matrix(0,128,7)
class()

set.seed(1234)
X[,1] <- u[,1]*(2-0.8)+0.8
min(X[,1]); max(X[,1]); hist(X[,1])
X[,2] <- floor(u[,2]*(10-6+1))+6
min(X[,2]); max(X[,2]); hist(X[,2])
X[,3] <- 2*qnorm(u[,3])+1
min(X[,3]); max(X[,3]); hist(X[,3])

#
a <-0
b<-20
c<-9
Ind1 <- which(u[,4]<(c-a)/(b-a))
Ind2 <- which(u[,4]>=(c-a)/(b-a))
X[Ind1,4]<-a+sqrt(u[Ind1,4]*(b-a)*(c-a))  
X[Ind2,4]<-b-sqrt((1-u[Ind2,4])*(b-a)*(b-c))  
min(X[,4]); max(X[,4]); hist(X[,4]) #

#
X[,5]<-qgamma(u[,5],2,0.5)
min(X[,5]); max(X[,5]); hist(X[,5])

#
Ind1 <- which(u[,7]<0.5)
Ind2 <- which(u[,7]>=0.5)
X[Ind1,7]<-1+log(2*(u[Ind1,7]))/2
X[Ind2,7]<-1-log(2*(1-u[Ind2,7]))/2
hist(X[,7])

############################  
library(EnvStats) # to use rtri
library(GGally)

runifdisc<-function(n, min, max) {
  sample(min:max, n, replace=T)
}

x1 <- runif(128, 0.8, 2)
x2 <- runifdisc(128, 6, 10)
x3 <- rnorm(128, 1, 2)
x4 <- rtri(128, 0.9, 20, 10.45)
x5 <- rgamma(128, 2, 2)
x6 <- rlnorm(128, 1, 2)
x7 <- 


m <- data.frame(x1,x2,x3,x4,x5,x6)


plot(m)





hist(x1, )
plot(density(x1))
