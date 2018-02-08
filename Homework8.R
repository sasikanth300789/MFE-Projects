
# Project 8
# Question1
setwd("/Users/sruthi/Desktop/237G")
rm(list = ls())

r0 <- 0.05
kappa <- 0.82
rbar <- 0.05
sigma <- 0.1
t <- 0
T <- 0.5

vasicek <- function(r0,kappa,sigma,rbar,start,end,paths){
  
days <- floor((end-start)*366)
dt <- 1/366

set.seed(0)
Z <- matrix(rnorm(paths*days,0,1),nrow=paths,ncol=days)

Rsim <- cbind(rep(r0,paths),matrix(NA,paths,days))

for (i in 1:paths){
  for (j in 2:(days+1)){
  
    Rsim[i,j] <- Rsim[i,j-1] + kappa*(rbar - Rsim[i,j-1] )*dt + sigma*sqrt(dt)*Z[i,j-1] 
    
  }
}

R <- -rowSums(Rsim[,-1])*dt
# Discount Value
Discount <- mean(exp(R))
return (Discount)
}

#1.a 
Price1a <- 1000*vasicek(r0,kappa,sigma,rbar,t,T,1000)

#1.b
# Coupon Paying Bond
C <- c(rep(30,7),1030)
T <- c(seq(0.5,4,0.5))
CPB <- sum(sapply(1:length(C), function(k)  C[k]*vasicek(r0,kappa,sigma,rbar,0,T[k],1000) ))

# Option Value
r0 <- 0.05
kappa <- 0.82
rbar <- 0.05
sigma <- 0.1
T<- 0.5
K <- 980
maturity <- 0.25
dt <- 1/366

Rsimulate <- function(rbar,kappa,sigma,T,t,r0,paths){

start <- t
end <- T
days <- floor((end-start)*366)
dt <- 1/366

set.seed(0)
Z <- matrix(rnorm(paths*days,0,1),nrow=paths,ncol=days)
Rsim <- cbind(rep(r0,paths),matrix(NA,paths,days))

for (i in 1:paths){
  for (j in 2:(days+1)){
    
    Rsim[i,j] <- Rsim[i,j-1] + kappa*(rbar - Rsim[i,j-1])*dt + sigma*sqrt(dt)*Z[i,j-1] 
    
  }
}
return (Rsim)
}

explicit <- function(rbar,sigma,kappa,T,t,r){
  B <- (1/kappa)*(1-exp(-kappa*(T-t)))
  A <- exp( ((rbar - 0.5*(sigma/kappa)^2)*(B-(T-t))) - (1/kappa)*((0.5*sigma*B)^2) )
  price <- A*exp(-B*r)
  return (price)
}

Rsim <- Rsimulate(rbar,kappa,sigma,maturity,0,r0,100)
R <- -rowSums(Rsim[,-1])*dt
optionprice <- mean(exp(R)*sapply(1:length(R),function(k) max(1000*explicit(rbar,sigma,kappa,T,maturity,Rsim[k,ncol(Rsim)])-K,0)))

# 1.d Monte Carlo Simulation
# the coupon paying bond will be given have different C and T vector at option maturity
T <- c(seq(0.5,4.0,0.5))
C <- c(rep(30,7),1030)
maturity <- 0.25

Rsim <- Rsimulate(rbar,kappa,sigma,maturity,0,r0,100)
R <- -rowSums(Rsim[,-1])*dt
price1d <- mean(exp(R) * (sapply(1:length(R), function(k) max(sum(sapply(1:length(C), function(i)  C[i]*vasicek( Rsim[k,ncol(Rsim)],kappa,sigma,rbar,maturity,T[i],100)))-K,0))) )

# 1.e
# explicit pricing

T <- c(seq(0.5,4.0,0.5))
C <- c(rep(30,7),1030)
maturity <- 0.25

Rsim <- Rsimulate(rbar,kappa,sigma,maturity,0,r0,100)
R <- -rowSums(Rsim[,-1])*dt
price1e <- mean(exp(R) * (sapply(1:length(R), function(k) max(sum(sapply(1:length(C), function(i)  C[i]*explicit(rbar,sigma,kappa,T[i],maturity, Rsim[k,ncol(Rsim)] )))-K,0))) )

# Question 2
# CIR model

r0 <- 0.05
kappa <- 0.92
rbar <- 0.055
sigma <- 0.12
t <- 0
T <- 0.5
S <- 1
dt <- 1/366
C <- 1000
K <- 980
dt <- 1/366

CIR <- function(r0,kappa,sigma,rbar,t,T){
  
  days <- floor((T-t)*366)
  paths <- 1000
  dt <- 1/366
  
  set.seed(0)
  Z <- matrix(rnorm(paths*days,0,1),nrow=paths,ncol=days)
  Rsim <- cbind(rep(r0,paths),matrix(NA,paths,days))
  
  for (i in 1:paths){
    for (j in 2:(days+1)){
        Rsim[i,j] <- Rsim[i,j-1] + kappa*(rbar - Rsim[i,j-1] )*dt + sqrt(Rsim[i,j-1])*sigma*sqrt(dt)*Z[i,j-1] 
    }
  }
  R <- -rowSums(Rsim[,-1])*dt
  # Discount Value
  Discount <- list(price = mean(exp(R)), Rsim = Rsim)
  return (Discount)
}

# Simulation till 0.5 years
Rsim <- CIR(r0,kappa,sigma,rbar,t,T)$Rsim
R <- -rowSums(Rsim[,-1])*dt
rT <- Rsim[,ncol(Rsim)]

price2a <- mean(exp(R)*sapply(1:length(R), function(k) max( (C*CIR(rT[k],kappa,sigma,rbar,T,S)$price-K),0)))

# 2.b
# Implicit Difference Method

explicitCIR <- function(rt,rbar,sigma,kappa,t,T){
  
  h1 <- sqrt( kappa^2 + 2*sigma^2 )
  h2 <- 0.5*(kappa+h1)
  h3 <- (2*kappa*rbar)/(sigma^2)  
  B <- (exp(h1*(T-t)) -1)/(h2*(exp(h1*(T-t)) -1) + h1)
  A <-  (h1*exp(h2*(T-t))/(h2*(exp(h1*(T-t)) -1) + h1))^h3 
  price <- A*exp(-B*rt)
  return (list(price=price,A=A,B=B))
}

require(matrixcalc)

IFD <- function(r0,sigma,kappa,rbar,K,T,S,type){

  # taking minimum value as 0
  dt <- 1/366
  dr <- r0/50
  R <- seq((r0+20*dr),0,-dr)
  N <- length(R)-1
  M <- T/dt
  
  # A matrix  
  A <- matrix(NA,nrow=N+1, ncol = N+1)
  A[1,] <- c(1,-1,rep(0,N-1))
  A[N+1,] <- c( rep(0,N-1),1,-1)
  
  for (i in 2:N){
    dummy <- c(rep(0,N+1))
    pu <- -(dt/(2*dr))*(kappa*(rbar - dr*(N+1-i)) + (sigma^2)*(N+1-i) )
    pm <- 1+(dt/dr)*(dr*dt*(N+1-i) +(sigma^2)*(N+1-i))
    pd <- -(dt/(2*dr))*(-kappa*(rbar - dr*(N+1-i)) + (sigma^2)*(N+1-i))
    dummy[(i-1):(i+1)] <- c(pu,pm,pd)
    A[i,] <- dummy
  }
  
  B <- vector("list", length = M+1)
  F <- vector("list", length = M+1)
  
  if (type == "call"){
    F[[M+1]] <- c(sapply(1:length(R), function(k) max(0,(1000*explicitCIR(R[k],rbar,sigma,kappa,T,S)$price)-K)))
    B[[M+1]] <- c(0,F[[M+1]][2:(length(R)-1)],1000*explicitCIR(R[1],rbar,sigma,kappa,T,S)$price-1000*explicitCIR(R[2],rbar,sigma,kappa,T,S)$price)
  } else {
    F[[M+1]] <- c(sapply(1:length(R), function(k) max(0,K-1000*explicitCIR(R[k],rbar,sigma,kappa,T,S)$price)))
    B[[M+1]] <- c(-1000*explicitCIR(R[length(R)-1],rbar,sigma,kappa,T,S)$price+1000*explicitCIR(R[length(R)],rbar,sigma,kappa,T,S)$price,F[[M+1]][2:(length(R)-1)],0)
  }
  
  for (j in M:1){
    
    if(type == "call") B[[j]] <- c(0,F[[j+1]][2:N],1000*explicitCIR(R[1],rbar,sigma,kappa,T,S)$price-1000*explicitCIR(R[2],rbar,sigma,kappa,T,S)$price) else B[[j]] <- c(-1000*explicitCIR(R[length(R)-1],rbar,sigma,kappa,T,S)$price+1000*explicitCIR(R[length(R)],rbar,sigma,kappa,T,S)$price,F[[j+1]][2:N],0)
    temp <- matrix.inverse(A) %*% B[[j]]
    F[[j]] <- pmax(temp, 0, na.rm = TRUE)
  }
  
  imputprice <- cbind(R,F[[1]])
  return(imputprice[which(R==r0),2])
}

price2b <- IFD(r0,sigma,kappa,rbar,K,T,S,"call")

# 2.c
# Explicit function
K <- K/1000

A <- explicitCIR(r0,rbar,sigma,kappa,T,S)$A
B <- explicitCIR(r0,rbar,sigma,kappa,T,S)$B

theta <- sqrt( kappa^2 + 2*sigma^2 )
phi <- (kappa+theta)/(sigma^2)
rstar <- log(A/K)/B
chi <- (2*theta)/((sigma^2)*(exp(theta*(T-t))-1))

price2c <- explicitCIR(r0,rbar,sigma,kappa,0,S)$price*pchisq(2*r0*(chi+phi+B), df = (4*kappa*rbar)/(sigma^2) , ncp = 2*(chi^2)*r0*exp(theta*(T-t))/(phi+chi+B) ) - K*explicitCIR(r0,rbar,sigma,kappa,0,T)$price*pchisq(2*r0*(phi+chi) ,df = 4*kappa*rbar/(sigma^2), ncp = 2*(chi^2)*r0*exp(theta*(T))/(phi+chi) )


# G2++ Model

a <- 0.1
b <- 0.3
sigma <- 0.03
x0 <- y0 <- 0
phi0 <- r0 <- 0.03
rho <- 0.7
eta <- 0.08
T <- 0.5
S <- 1
t <- 0
dt <- 1/366
K <- 950/1000

G2 <- function(a,b,sigma,x0,y0,phi0,r0,rho,eta,T,t,K){
  
  days <- floor((T-t)*366)
  paths <- 1000
  dt <- 1/366
  
  Z1 <- matrix(rnorm(paths*days,0,1),nrow=paths,ncol=days)
  Z2 <- matrix(rnorm(paths*days,0,1),nrow=paths,ncol=days)
  Z3 <- rho*Z1 + sqrt(1-(rho^2))*Z2
  
  Xsim <- cbind(rep(x0,paths),matrix(NA,paths,days))
  Ysim <- cbind(rep(y0,paths),matrix(NA,paths,days))
  Rsim <- cbind(rep(r0,paths),matrix(NA,paths,days))
  
  for (i in 1:paths){
    for (j in 2:(days+1)){
      Xsim[i,j] <- Xsim[i,j-1] - a*Xsim[i,j-1]*dt + sigma*sqrt(dt)*Z1[i,j-1] 
      Ysim[i,j] <- Ysim[i,j-1] - b*Ysim[i,j-1]*dt + eta*sqrt(dt)*Z3[i,j-1]
      Rsim[i,j] <- Xsim[i,j] + Ysim[i,j] + phi0
    }
  }
  R <- -rowSums(Rsim[,-1])*dt
  # Discount Value
  Discount <- list(price = mean(exp(R)), R = R, Rsim = Rsim, Xsim = Xsim, Ysim = Ysim)
  return (Discount)
}

#simulation till T = 0.5
G2sim <- G2(a,b,sigma,x0,y0,phi0,r0,rho,eta,T,t,K)
Rsim <- G2sim$Rsim
Xsim <- G2sim$Xsim
Ysim <- G2sim$Ysim
R <- G2sim$R

xT <- Xsim[,ncol(Xsim)]
yT <- Ysim[,ncol(Ysim)]
rT <- Rsim[,ncol(Rsim)]

price3a <- 1000*mean( exp(R)*(sapply(1:length(R), function(k) max(K-(G2(a,b,sigma,xT[k],yT[k],phi0,rT[k],rho,eta,S,T,K)$price),0))) )

# 3.b
# Explicit Formula
Sigma <- sqrt(((sigma^2)/(2*a^3))*((1-exp(-a*(S-T)))^2)*(1-exp(-2*a*(T-t))) + ((eta^2)/(2*b^3))*((1-exp(-b*(S-T)))^2)*(1-exp(-2*b*(T-t))) + (2*rho*(sigma*eta/(a*b*(a+b))))*(1-exp(-b*(S-T)))*(1-exp(-a*(S-T)))*(1-exp(-(a+b)*(T-t))))
  
priceT <- G2(a,b,sigma,x0,y0,phi0,r0,rho,eta,T,t,K)$price
priceS <- G2(a,b,sigma,x0,y0,phi0,r0,rho,eta,S,t,K)$price
price3b <- 1000*(-priceS*pnorm((log(K*priceT/priceS)/Sigma) - Sigma*0.5)) + 1000*(priceT*K*pnorm((log(K*priceT/priceS)/Sigma) + Sigma*0.5))
  
  

  

