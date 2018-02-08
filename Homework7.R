
# Question 1
setwd("/Users/sruthi/Desktop/237G")
rm(list = ls())

# 1.a explicit method function
explicit <- function(Smax,Smin,sigma,r,dt,dx,K,T) {

# probabilities  
pu <- dt*((sigma^2)/(2*dx^2) + (r-0.5*sigma^2)/(2*dx))
pm <- 1-dt*(sigma^2)/(dx^2) - r*dt
pd <- dt*((sigma^2)/(2*dx^2) - (r-0.5*sigma^2)/(2*dx))
X <- seq(log(Smax),log(Smin),-dx)
S <- exp(seq(log(Smax),log(Smin),-dx))

N <- 0.5*(length(S) - 1)
M <- T/dt

# A matrix 
A <- matrix(NA,nrow=(2*N+1), ncol = (2*N+1))
A[1,] <- c(pu,pm,pd,rep(0,2*(N-1)))
A[(2*N+1),] <- c(rep(0,2*(N-1)),pu,pm,pd)

for (i in 2:(2*N)){
  dummy <- c(rep(0,length(S)))
  dummy[(i-1):(i+1)] <- c(pu,pm,pd)
  A[i,] <- dummy
}

B <- vector("list", length = M+1)
F <- vector("list", length = M+1)
B[[M+1]] <- c(rep(0,2*N),(S[length(S)-1]-S[length(S)]))
F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,K-S[k])))

# matrix solution
for (j in M:1){
  B[[j]] <- c(rep(0,2*N),(S[length(S)-1]-S[length(S)]))
  F[[j]] <- A %*% F[[(j+1)]] + B[[j]]
}

explicitprices <- cbind(S,F[[1]])
return (explicitprices)
}

#Implicit Method function
implicit <- function(Smax,Smin,sigma,r,dt,dx,K,T){

# probabilities  
ipu <- -dt*((sigma^2)/(2*dx^2) + (r-0.5*sigma^2)/(2*dx))
ipm <- 1 + dt*(sigma^2)/(dx^2) + r*dt
ipd <- -dt*((sigma^2)/(2*dx^2) - (r-0.5*sigma^2)/(2*dx))
X <- seq(log(Smax),log(Smin),-dx)
S <- exp(seq(log(Smax),log(Smin),-dx))

N <- 0.5*(length(S) - 1)
M <- T/dt

# A matrix
A <- matrix(NA,nrow=(2*N+1), ncol = (2*N+1))
A[1,] <- c(1,-1,rep(0,(2*N-1)))
A[(2*N+1),] <- c(rep(0,(2*N-1)),1,-1)

for (i in 2:(2*N)){
  dummy <- c(rep(0,length(S)))
  dummy[(i-1):(i+1)] <- c(ipu,ipm,ipd)
  A[i,] <- dummy
}

B <- vector("list", length = M+1)
F <- vector("list", length = M+1)

F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,K-S[k])))
B[[M+1]] <- c(0,F[[M+1]][2:(length(S)-1)] ,-(S[length(S)-1]-S[length(S)]))

# matrix solutions
require(matrixcalc)

for (j in M:1){
  B[[j]] <- c(0,F[[j+1]][2:(2*N)],-(S[length(S)-1]-S[length(S)]))
  F[[j]] <- matrix.inverse(A) %*% B[[j]]
}

implicitprice <- cbind(S,F[[1]])

return (implicitprice)
}

# Crank-Nicolson Finite Difference 

cranknicolson <- function(Smax,Smin,sigma,r,dt,dx,K,T){

# probabilities  
cpu <- -0.25*dt*((sigma^2)/(dx^2) + (r-0.5*sigma^2)/(dx))
cpm <- 1 + 0.5*dt*(sigma^2)/(dx^2) + 0.5*r*dt
cpd <- -0.25*dt*((sigma^2)/(dx^2) - (r-0.5*sigma^2)/(dx))
X <- seq(log(Smax),log(Smin),-dx)
S <- exp(seq(log(Smax),log(Smin),-dx))

N <- 0.5*(length(S) - 1)
M <- T/dt

# A matrix
A <- matrix(NA, nrow=(2*N+1), ncol=(2*N+1))
A[1,] <- c(1,-1,rep(0,(2*N-1)))
A[(2*N+1),] <- c(rep(0,(2*N-1)),1,-1)

for (i in 2:(2*N)){
  dummy <- c(rep(0,length(S)))
  dummy[(i-1):(i+1)] <- c(cpu,cpm,cpd)
  A[i,] <- dummy
}

F <- vector("list", length = M+1)

F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,K-S[k])))

# matrix solutions
require(matrixcalc)

for (i in M:1){
  
  z <- c(rep(NA,length(S)))
  z[1] <- 0
  z[length(S)] <- -(S[length(S)-1] - S[length(S)])
  z[2:(length(S)-1)] <- c(sapply(2:(length(S)-1), function(k) -cpu*F[[i+1]][k-1] - (cpm-2)*F[[i+1]][[k]] - cpd*F[[i+1]][[k+1]]))
  F[[i]] <- matrix.inverse(A) %*% z

}

cnputprice <- cbind(S,F[[1]])
return(cnputprice)
}

#1.a
#main code

S0 <- 10
Smax <- 16
Smin <- 4
sigma <- 0.2
r <- 0.04
dt <- 0.002
dx <- c(sigma*sqrt(dt),sigma*sqrt(3*dt),sigma*sqrt(4*dt))
K <- 10
T <- 0.5

explicitprices <- vector("list", length = 3)
implicitprices <- vector("list", length = 3)
cnprices <- vector("list", length = 3)

for (j in 1:length(dx)){
explicitprices[[j]] <- explicit(16,4,0.2,0.04,0.002,dx[j],K,T)
implicitprices[[j]] <- implicit(16,4,0.2,0.04,0.002,dx[j],K,T)
cnprices[[j]] <- cranknicolson(16,4,0.2,0.04,0.002,dx[j],K,T)
}

Pa <- sapply(1:3, function(k) explicitprices[[k]][ which.min(abs(explicitprices[[k]][,1]-S0)),2] )
Pb <- sapply(1:3, function(k) implicitprices[[k]][ which.min(abs(implicitprices[[k]][,1]-S0)),2] )
Pc <- sapply(1:3, function(k) cnprices[[k]][ which.min(abs(cnprices[[k]][,1]-S0)),2] )

P <- rbind(Pa,Pb,Pc)
colnames(P) <- c("dx1","dx2","dx3")
rownames(P) <- c("Pa","Pb","Pc")

# output european put prices
#1.a
P

#1.b
require(OptionPricing)
current <- seq(4,16,1)

# EFD method erros
EFDerrors <- cbind(current,NA,NA,NA)
colnames(EFDerrors) <- c("Price","dx1_error","dx2_error", "dx3-error")
for (j in 1:nrow(EFDerrors)){

  blsprice <- BS_EP(T,K,r,sigma,current[j])[1]
  EFDerrors[j,2] <- explicitprices[[1]][ which.min(abs(explicitprices[[1]][,1]-current[j])),2] - blsprice
  EFDerrors[j,3] <- explicitprices[[2]][ which.min(abs(explicitprices[[2]][,1]-current[j])),2] - blsprice
  EFDerrors[j,4] <- explicitprices[[3]][ which.min(abs(explicitprices[[3]][,1]-current[j])),2] - blsprice
}

# IFD method erros
IFDerrors <- cbind(current,NA,NA,NA)
colnames(IFDerrors) <- c("Price","dx1_error","dx2_error", "dx3-error")
for (j in 1:nrow(IFDerrors)){
  
  blsprice <- BS_EP(T,K,r,sigma,current[j])[1]
  IFDerrors[j,2] <- implicitprices[[1]][ which.min(abs(implicitprices[[1]][,1]-current[j])),2] - blsprice
  IFDerrors[j,3] <- implicitprices[[2]][ which.min(abs(implicitprices[[2]][,1]-current[j])),2] - blsprice
  IFDerrors[j,4] <- implicitprices[[3]][ which.min(abs(implicitprices[[3]][,1]-current[j])),2] - blsprice
}

# CN method erros
CNerrors <- cbind(current,NA,NA,NA)
colnames(CNerrors) <- c("Price","dx1_error","dx2_error", "dx3-error")
for (j in 1:nrow(CNerrors)){
  
  blsprice <- BS_EP(T,K,r,sigma,current[j])[1]
  CNerrors[j,2] <- cnprices[[1]][ which.min(abs(cnprices[[1]][,1]-current[j])),2] - blsprice
  CNerrors[j,3] <- cnprices[[2]][ which.min(abs(cnprices[[2]][,1]-current[j])),2] - blsprice
  CNerrors[j,4] <- cnprices[[3]][ which.min(abs(cnprices[[3]][,1]-current[j])),2] - blsprice
}


# Question 2
rm(list = ls())

#American Put and Call

# 1.a explicit method function

EFD <- function(S0,sigma,r,dt,ds,K,T,type){

# minimum value is 0
S <- seq((S0+10*ds),0,-ds)
N <- length(S)-1
M <- T/dt

# A matrix
A <- matrix(NA,nrow=N+1, ncol = N+1)
A[1,] <- c(dt*(0.5*r*(N-1) + 0.5*(sigma*(N-1))^2 ), 1-dt*(r+(sigma*(N-1))^2), dt*(-0.5*r*(N-1) + 0.5*(sigma*(N-1))^2 ),rep(0,N-2))
A[N+1,] <- c( rep(0,N-2), dt*(0.5*r + 0.5*(sigma)^2), 1-dt*(r+(sigma)^2), dt*(-0.5*r + 0.5*(sigma)^2))

for (i in 2:N){
  dummy <- c(rep(0,N+1))
  pu <- dt*(0.5*r*(N+1-i) + 0.5*(sigma*(N+1-i))^2 )
  pm <- 1-dt*(r+(sigma*(N+1-i))^2 )
  pd <- dt*(-0.5*r*(N+1-i) + 0.5*(sigma*(N+1-i))^2 )
  dummy[(i-1):(i+1)] <- c(pu,pm,pd)
  A[i,] <- dummy
}

B <- vector("list", length = M+1)
F <- vector("list", length = M+1)

if (type == "call"){
  B[[M+1]] <- c((S[1]-S[2]),rep(0,N))
  F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,S[k]-K)))
} else {
  B[[M+1]] <- c(rep(0,N),(S[length(S)-1]-S[length(S)]))
  F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,K-S[k]))) 
}

# matrix solution
for (j in M:1){
  
  if (type=="call") B[[j]] <- c((S[1]-S[2]),rep(0,N)) else B[[j]] <- c(rep(0,N),(S[length(S)-1]-S[length(S)]))
  temp <- A %*% F[[(j+1)]] + B[[j]]
  #american option condition
  F[[j]] <- pmax(temp, F[[j+1]], na.rm = TRUE)
}

exputprice <- cbind(S,F[[1]])
return (exputprice[which(S==S0),2])
}

#Implicit Finite Method function

IFD <- function(S0,sigma,r,dt,ds,K,T,type){

# taking minimum value as 0
S <- seq((S0+10*ds),0,-ds)
N <- length(S)-1
M <- T/dt

# A matrix  
A <- matrix(NA,nrow=N+1, ncol = N+1)
A[1,] <- c(1,-1,rep(0,N-1))
A[N+1,] <- c( rep(0,N-1),1,-1)

for (i in 2:N){
  dummy <- c(rep(0,N+1))
  
  pu <- -dt*(0.5*r*(N+1-i) + 0.5*(sigma*(N+1-i))^2 )
  pm <- 1+dt*(r+(sigma*(N+1-i))^2 )
  pd <- -dt*(-0.5*r*(N+1-i) + 0.5*(sigma*(N+1-i))^2 )
  dummy[(i-1):(i+1)] <- c(pu,pm,pd)
  A[i,] <- dummy
}

B <- vector("list", length = M+1)
F <- vector("list", length = M+1)

if (type == "call"){
  F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,S[k]-K)))
  B[[M+1]] <- c(S[1]-S[2],F[[M+1]][2:(length(S)-1)],0)
} else {
  F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,K-S[k])))
  B[[M+1]] <- c(0,F[[M+1]][2:(length(S)-1)] ,-(S[length(S)-1]-S[length(S)]))
}

for (j in M:1){
  
  if(type == "call") B[[j]] <- c(S[1]-S[2],F[[j+1]][2:N],0) else B[[j]] <- c(0,F[[j+1]][2:N],-(S[length(S)-1]-S[length(S)]))
  temp <- matrix.inverse(A) %*% B[[j]]
  F[[j]] <- pmax(temp, F[[j+1]], na.rm = TRUE)
}

imputprice <- cbind(S,F[[1]])
return (imputprice[which(S==S0),2])

}

# Crank-Nicolson function
CNFD <- function(S0,sigma,r,dt,ds,K,T,type){

S <- seq((S0+10*ds),0,-ds)
N <- length(S)-1
M <- T/dt
  
A <- matrix(NA,nrow=N+1, ncol = N+1)
A[1,] <- c(1,-1,rep(0,N-1))
A[N+1,] <- c( rep(0,N-1),1,-1)

for (i in 2:N){
  dummy <- c(rep(0,N+1))
  
  pu <- -0.25*dt*( r*(N+1-i) + (sigma*(N+1-i))^2 )
  pm <- 1+0.5*dt*( r+(sigma*(N+1-i))^2 )
  pd <- -0.25*dt*( -r*(N+1-i) + (sigma*(N+1-i))^2 )
  dummy[(i-1):(i+1)] <- c(pu,pm,pd)
  A[i,] <- dummy
}

F <- vector("list", length = M+1)

if (type=="call") { 
  F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,S[k]-K))) 
} else { 
  F[[M+1]] <- c(sapply(1:length(S), function(k) max(0,K-S[k]))) 
}

for (j in M:1){
  
  if (type == "call"){
    z <- c(rep(NA,length(S)))
    z[1] <- S[1]-S[2]
    z[length(S)] <- 0
    z[2:N] <- c(sapply(2:N, function(k) -(-0.25*dt*( r*(N+1-k) + (sigma*(N+1-k))^2 ))*(F[[j+1]][k-1]) - ((1+ 0.5*dt*( r+(sigma*(N+1-k))^2 ))- 2)*(F[[j+1]][[k]]) - (-0.25*dt*(-r*(N+1-k) + (sigma*(N+1-k))^2 ))*(F[[j+1]][[k+1]])))
  } else{
    z <- c(rep(NA,length(S)))
    z[1] <- 0
    z[length(S)] <- -(S[length(S)-1] - S[length(S)])
    z[2:N] <- c(sapply(2:N, function(k) -(-0.25*dt*( r*(N+1-k) + (sigma*(N+1-k))^2 ))*(F[[j+1]][k-1]) - ((1+ 0.5*dt*( r+(sigma*(N+1-k))^2 ))- 2)*(F[[j+1]][[k]]) - (-0.25*dt*(-r*(N+1-k) + (sigma*(N+1-k))^2 ))*(F[[j+1]][[k+1]])))
  }
  temp <- matrix.inverse(A) %*% z
  F[[j]] <- pmax(temp, F[[j+1]], na.rm = TRUE)
  
}

cnputprice <- cbind(S,F[[1]])
return (cnputprice[which(S==S0),2])
}

# Main code

S0 <- 10
sigma <- 0.2
r <- 0.04
dt <- 0.002
K <- 10
T <- 0.5
ds <- c(0.5,1,1.5)

current <- seq(16,4,-1)
optprices <- cbind(current,NA,NA,NA)
colnames(optprices) <- c("Stock","EFD","IFD","CNFD")

final <- list(ds1=list(call = optprices, put = optprices),ds2=list(call = optprices, put = optprices),ds3=list(call = optprices, put = optprices))

for (j in 1:length(ds)){
  for (i in 1:length(current)){
    
  final[[j]]$call[i,2] <- EFD(current[i],0.2,0.04,0.002,ds[j],K,T,"call")
  final[[j]]$call[i,3] <- IFD(current[i],0.2,0.04,0.002,ds[j],K,T,"call")
  final[[j]]$call[i,4] <- CNFD(current[i],0.2,0.04,0.002,ds[j],K,T,"call")
  
  final[[j]]$put[i,2] <- EFD(current[i],0.2,0.04,0.002,ds[j],K,T,"put")
  final[[j]]$put[i,3] <- IFD(current[i],0.2,0.04,0.002,ds[j],K,T,"put")
  final[[j]]$put[i,4] <- CNFD(current[i],0.2,0.04,0.002,ds[j],K,T,"put")
  }
}

Pa <- sapply(1:3, function(k) EFD(S0,0.2,0.04,0.002,ds[k],K,T,"put") )
Pb <- sapply(1:3, function(k) IFD(S0,0.2,0.04,0.002,ds[k],K,T,"put") )
Pc <- sapply(1:3, function(k) CNFD(S0,0.2,0.04,0.002,ds[k],K,T,"put") )

Ca <- sapply(1:3, function(k) EFD(S0,0.2,0.04,0.002,ds[k],K,T,"call") )
Cb <- sapply(1:3, function(k) IFD(S0,0.2,0.04,0.002,ds[k],K,T,"call") )
Cc <- sapply(1:3, function(k) CNFD(S0,0.2,0.04,0.002,ds[k],K,T,"call") )

P <- rbind(Pa,Pb,Pc)
colnames(P) <- c("ds1","ds2","ds3")
rownames(P) <- c("Pa","Pb","Pc")

C <- rbind(Ca,Cb,Cc)
colnames(C) <- c("ds1","ds2","ds3")
rownames(C) <- c("Ca","Cb","Cc")

#output
#call prices
C
#put prices
P

#plots
require(ggplot2)
callprice <- data.frame(cbind(final[[1]]$call,final[[2]]$call,final[[3]]$call))

ggplot(callprice, aes(x=Stock)) + geom_line(aes(y=EFD, colour="0.5 EFD method")) + 
 geom_line(aes(y=IFD, colour="0.5 IFD method")) +
 geom_line(aes(y=CNFD, colour="0.5 CNFD method")) +
 geom_line(aes(y=EFD.1, colour="1.0 EFD method")) +
 geom_line(aes(y=IFD.1, colour="1.0 IFD method")) +
 geom_line(aes(y=CNFD.1, colour="1.0 CNFD method")) +
 geom_line(aes(y=EFD.2, colour="1.5 EFD method")) +
 geom_line(aes(y=IFD.2, colour="1.5 IFD method")) +
 geom_line(aes(y=CNFD.2, colour="1.5 CNFD method")) +
 xlab("Stock Prices") + ylab("option price") + ggtitle("American Call Option Prices")

putprice <- data.frame(cbind(final[[1]]$put,final[[2]]$put,final[[3]]$put))
ggplot(putprice, aes(x=Stock)) + geom_line(aes(y=EFD, colour="0.5 EFD method")) + 
  geom_line(aes(y=IFD, colour="0.5 IFD method")) +
  geom_line(aes(y=CNFD, colour="0.5 CNFD method")) +
  geom_line(aes(y=EFD.1, colour="1.0 EFD method")) +
  geom_line(aes(y=IFD.1, colour="1.0 IFD method")) +
  geom_line(aes(y=CNFD.1, colour="1.0 CNFD method")) +
  geom_line(aes(y=EFD.2, colour="1.5 EFD method")) +
  geom_line(aes(y=IFD.2, colour="1.5 IFD method")) +
  geom_line(aes(y=CNFD.2, colour="1.5 CNFD method")) +
  xlab("Stock Prices") + ylab("option price") + ggtitle("American Put Option Prices")

