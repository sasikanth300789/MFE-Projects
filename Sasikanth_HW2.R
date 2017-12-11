
# defining fucntions

lgm <- function(n,x0){
  a <- 7^5
  b <- 0
  m <- 2^31 - 1
  rand <- numeric(n)
  rand[1] <- x0
  for (i in 2:n) {
    rand[i] = (rand[i-1]*a)%%m
  }
  return (rand/m) 
}

simnorm <- function(n,m){
  unimatrix <- matrix(lgm(n,m),n/2,2)
  Zone <- sqrt(-2*log(unimatrix[,1]))*cos(2*pi*unimatrix[,2])
  Ztwo <- sqrt(-2*log(unimatrix[,1]))*sin(2*pi*unimatrix[,2])
  Z <- c(Zone, Ztwo)
  return (Z)
}

simmvn <- function(n,m,mu,Sigma){
  rho <- Sigma[1,2]
  Z <- matrix(simnorm(2*n,m),ncol=2)
  X <- mu[1] + Sigma[1,1]*Z[,1]
  Y <- mu[2] + rho*Sigma[2,2]*X + sqrt(1-rho^2)*Sigma[2,2]*Z[,2]
  return (cbind(X,Y))
}

#Q1  Correlation
Sigma <- matrix(c(1,-0.7,-0.7,1),2,2)
mu <- c(0,0)
data <- matrix(simmvn(n = 10000, 0.001,mu, Sigma), ncol =2)
xbar <- data[,1] - mean(data[,1])
ybar <- data[,2] - mean(data[,2])
correlation <- sum(xbar*ybar)/sqrt(sum(xbar*xbar)*sum(ybar*ybar))

#Q2 Monte Carlo Simulation
#library(MASS)
Sigma <- matrix(c(1,0.6,0.6,1),2,2)
mu <- c(0,0)
data<- matrix(simmvn(n = 10000,0.001, mu, Sigma),ncol=2)

func <-function(X,Y){
  return (X^3 + sin(Y) + Y*X^2)
}

temp = 0
for (i in 1:nrow(data)){
  temp = temp + max(0,func(data[i,1],data[i,2]))  
}

E <- temp/i

# Question 3
n <- 100000
rand <- simnorm(n,2)
a1 <- 5*rand^2 + sin(rand*sqrt(5))
Ea1 <- sum(a1)/n
a2 <- exp(0.5/2) * cos(rand*sqrt(0.5))
Ea2 <- sum(a2)/n
a3 <- exp(3.2/2) * cos(rand*sqrt(3.2))
Ea3 <- sum(a3)/n
a4 <- exp(6.5/2) * cos(rand*sqrt(6.5))
Ea4 <- sum(a4)/n

var(a1)
var(a2)
var(a3)
var(a4)

#Q3.c
ran <- -rand
b1 <- (5*ran^2 + sin(ran*sqrt(5)) + a1)/2
Eb1 <- sum(b1)/n
b2 <- (exp(0.5/2) * cos(ran*sqrt(0.5)) + a2)/2
Eb2 <- sum(b2)/n
b3 <- (exp(3.2/2) * cos(ran*sqrt(3.2)) + a3)/2
Eb3 <- sum(b3)/n
b4 <- (exp(6.5/2) * cos(ran*sqrt(6.5)) + a4)/2
Eb4 <- sum(b4)/n

var(b1)
var(b2)
var(b3)
var(b4)

#Question 4
n <- 10000
r <- 0.04
sigma <- 0.2
K <- 100
s0 <- 88
T <- 5
rand <- simnorm(n,0.001)
price <- numeric()
for (i in 1:length(rand)){
  price[i] <- exp(-r*T)*max(0,s0*(exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*rand[i])) - K) 
}
Ca1 <- mean(price)
var(price)

#4.b
library("OptionPricing")
price2 <- numeric()
for (j in 1:length(rand)){
  price2[j] <- 0.5*(price[j] + exp(-r*T)*max(0,s0*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*(-rand[j])) - K))
}
Ca2 <- mean(price2)
var(price2)

blsprice <- BS_EC(5,100,0.04,0.2,88)[1]


# Question 5
S0 <- 88
sigma <- 0.18
r <- 0.04
n <- 10
nsims <- 10000
S <- matrix(NA,nrow=nsims,ncol=10)
for (i in 1:n){
  Z <- simnorm(nsims,0.01)
  S[,i] <- S0 * exp(sigma*sqrt(i)*Z + (r-0.5*sigma^2)*i)   
}
ES <- c(S0,colMeans(S))
plot(0:n,ES,ylim=c(0,300),main="Stock Price Path", xlab= "Time", ylab="Price", type ='o', col = "red", pch=10, lty = 1, lwd =2)
legend(8,300,c("E(S)"),col=("red"), lty=c(1), lwd=c(2.5))

## simulate stock paths
library(sde)
suppressWarnings(suppressMessages(library(sde)))
t<- n/nsims
simS <- matrix(NA,ncol=7, nrow = nsims+1)
simS[,7] <- seq(0,10,10/nsims)
col <- c("lightblue","yellow","green","purple","brown","orange")
for (j in 1:6) {
  simS[,j] <- S0*cumprod(1+ c(0,(r-sigma*sigma*0.5)*t + sigma*sqrt(t)*simnorm(nsims,j/10)))
  lines(simS[,7],simS[,j],type='l', col=col[j])
}

#5.d
sigma <- 0.35
nsims <- 10000
S <- matrix(NA,nrow=nsims,ncol=10)
for (i in 1:n){
  Z <- simnorm(nsims,0.01)
  S[,i] <- S0 * exp(sigma*sqrt(i)*Z + (r-0.5*sigma^2)*i)   
}
ES <- c(S0,colMeans(S))
plot(0:n,ES,ylim=c(0,300),main="Stock Price Path", xlab= "Time", ylab="Price", type ='o', col = "red", pch=10, lty = 1, lwd =2)
legend(8,300,c("E(S)"),col=("red"), lty=c(1), lwd=c(2.5))

## simulate stock paths
library(sde)
t<- n/nsims
simS <- matrix(NA,ncol=7, nrow = nsims+1)
simS[,7] <- seq(0,10,10/nsims)
col <- c("lightblue","yellow","green","purple","brown","orange")
for (j in 1:6) {
  simS[,j] <- S0*cumprod(1+ c(0,(r-sigma*sigma*0.5)*t + sigma*sqrt(t)*simnorm(nsims,j/10)))
  lines(simS[,7],simS[,j],type='l', col=col[j])
}
  


#Question 6
# Euler Method
n <- 10000
tmax <- 1
h <- tmax/n
time <- seq(0,1,1/n)
F <- numeric()
F[1]<-0
for (j in 1:(n+1)){
  F[j+1] = F[j] + 4*h*sqrt(1-time[j]^2)
}
Ia <- F[j]

#6.b Monte Carlo
uniform <- lgm(n,0.01)
M <- numeric()
for (i in 1:length(uniform)){
  M[i] <- 4*sqrt(1-uniform[i]^2)
}
Ib <- mean(M)

#6.c Importance Sampling
n <- 10000
x<- lgm(n,1)
a<- 0.74
# used for g(x) function
max <- 1/(1-a/3)
#since the range is (0,1) then h(x) is uniform
h <- lgm(n,2)
sample <- numeric()
for(i in 1:n){
  if( x[i] <= (1-a*h[i]^2)){
    sample <- c(sample,h[i])
  }
}
# after the sample for t(x) is obtained, use them to calculate the expectation value
t <- numeric()
for (i in 1:length(sample)) {
  t[i] <- 4*(1-a/3)*sqrt(1-sample[i]^2)/(1-a*sample[i]^2)
}
Ic <- mean(t)