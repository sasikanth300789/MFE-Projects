# Set directory
setwd("/Users/sruthi/Desktop")
rm(list = ls())

# Question 1

r <- 0.03
S0 <- 98
X <- 100
T <- 1
sigma <- seq(0.12,0.48,0.04)

Proj6_1func <- function(S0,r,sigma,T,X){
  n <- 1000
  m <- 1000
  
  set.seed(0)
  Z = matrix(rnorm(n*m),n,m) 
  
  # Simulate stock paths
  normal <- matrix(exp((r-0.5*sigma^2)*(T/m) + sigma*sqrt(T/m)*Z),nrow=n,ncol=m)
  pr_norm <- cbind(rep(1,n),normal)
  pr <- S0*t(apply(pr_norm,1,cumprod))
  
  Smax <- apply(pr,1,max)
  Smin <- apply(pr,1,min)
  return (list(Smax = Smax, Smin = Smin))
}

finalmax <- cbind(sapply(1:length(sigma), function(j) Proj6_1func(S0,r,sigma[j],T,X)$Smax ))
finalmin <- cbind(sapply(1:length(sigma), function(j) Proj6_1func(S0,r,sigma[j],T,X)$Smin ))

finalmax <- finalmax - X
finalmax[c(finalmax<=0)] <- 0
Callprice <- exp(-r*T)*apply(finalmax,2,mean)
png("Proj6_1a.png", width=6, height=4, units="in", res=200)
plot(sigma,Callprice, xlab="Volatility", ylab="Price", main='Call Option Price', type='o', col='Red')
dev.off()

finalmin <- X - finalmin
finalmin[c(finalmin<=0)] <- 0
Putprice <- exp(-r*T)*apply(finalmin,2,mean)
png("Proj6_1b.png", width=6, height=4, units="in", res=200)
plot(sigma,Putprice, xlab="Volatility", ylab="Price", main='Put Option Price', type='o', col='Blue')
dev.off()

# Question 2
# Jump Diffusion Process

Proj6_2func <- function(lambda1,lambda2,T){

# Define parameters  
n <- 1000
m <- 1000
r0 <- 0.02
sigma <- 0.2
mu <- -0.1
gamma <- -0.4
V0 <- 20000
L0 <- 22000
delta <- 0.25
R <- r0 + delta*lambda2
r <- R/12
N <- T*12
PMT <- L0*r/(1-(1/(1+r)^N))
a <- PMT/r
b <- (PMT/r)*(1/(1+r)^N)
c <- 1+r
epsilon <- 0.95
alpha <- 0.7
beta <- (epsilon-alpha)/T
t <- seq(0,T,T/m)
q <- alpha + beta*t
L <- a - b*(c^(12*t))

# Change the value for lambda2 when it is zero. Because rexp() fucntion in R gives NAs when lambda = 0
if (lambda2 == 0) {lambda2 <- 0.001}

# Generate random normal numbers
set.seed(0)
Z <- cbind(rep(1,n),matrix(rnorm(n*m),n,m))
V <- cbind(rep(V0,n),matrix(NA,n,m))

option <- rep(0,n)
exercise <- rep(0,n)

# loop over n paths
for (i in 1:n){
  jumps <- rexp(100,lambda1)
  jumps <- as.numeric(jumps[cumsum(jumps) <= T])
  if(length(jumps) == 0) {jumps <- 5.1}
  jumpstep <- ceiling(cumsum(jumps)*m/T)+1
  stoppage <- ceiling(rexp(1,lambda2)*m/T)+1

# loop over m time steps    
for (j in 2:(m+1)){
    if (j %in% jumpstep) {
      #V[i,j] = V[i,j-1]*exp((mu-0.5*sigma^2)*(T/m) + sigma*sqrt(T/m)*Z[i,j])*exp(gamma)
      V[i,j] = V[i,j-1]*(1 + mu*T/m + sigma*sqrt(T/m)*Z[i,j] + gamma)
    } else { 
      V[i,j] = V[i,j-1]*(1 + mu*T/m + sigma*sqrt(T/m)*Z[i,j])
      #V[i,j] = V[i,j-1]*exp((mu-0.5*sigma^2)*(T/m) + sigma*sqrt(T/m)*Z[i,j]) 
    }

    #check for Q value time step
    if (V[i,j] <= q[j]*L[j]) {
      Q <- j
    
    #check if option is exercised
    tau <- min(Q,stoppage)
      
      if(tau < m+1 ){ 
        if(tau == Q) { 
          option[i] <- (exp(-r0*t[tau]))*max(0,(L[tau]-epsilon*V[i,tau]))  
        } else { 
          option[i] <- (exp(-r0*t[tau]))*abs(L[tau]-epsilon*V[i,tau]) 
        }   
      exercise[i] <- tau
      break
      }
      
    }
    if(j==stoppage){
        option[i] <- (exp(-r0*t[j]))*abs(L[j]-epsilon*V[i,j]) 
        exercise[i] <- j
    }
  
  }
}

optionvalue <- mean(option)
defprob <- 1- sum(option==0)/n
exptime <- (T/m)*mean(exercise[exercise < m+1]) 

my_list <- list(optionvalue = optionvalue, defprob = defprob, exptime = exptime)
return(my_list) 
}

# declare the input vectors
lamb1 <- seq(0.05,0.4,0.05) 
lamb2 <- seq(0,0.8,0.1) 
T <- seq(3,8,1) 

#2.a
D1 <- sapply(1:length(lamb1), function(k) Proj6_2func(lamb1[k],0.4,5)$optionvalue)
D2 <- sapply(1:length(lamb2), function(k) Proj6_2func(0.2,lamb2[k],5)$optionvalue)
D3 <- sapply(1:length(T), function(k) Proj6_2func(0.2,0.4,T[k])$optionvalue)

#2.b
Prob1 <- sapply(1:length(lamb1), function(k) Proj6_2func(lamb1[k],0.4,5)$defprob)
Prob2 <- sapply(1:length(lamb2), function(k) Proj6_2func(0.2,lamb2[k],5)$defprob)
Prob3 <- sapply(1:length(T), function(k) Proj6_2func(0.2,0.4,T[k])$defprob)

#2.csda
Et1 <- sapply(1:length(lamb1), function(k) Proj6_2func(lamb1[k],0.4,5)$exptime)
Et2 <- sapply(1:length(lamb2), function(k) Proj6_2func(0.2,lamb2[k],5)$exptime)
Et3 <- sapply(1:length(T), function(k) Proj6_2func(0.2,0.4,T[k])$exptime)

# Define parameters for graphs
par(mar=c(4,6,2,6))
par(mfrow=c(2,1))

# plots for (a) part
label <- c("Lambda2=0","Lambda2=0.1","Lambda2=0.2","Lambda2=0.3","Lambda2=0.4","Lambda2=0.5","Lambda2=0.6","Lambda2=0.7","Lambda2=0.8")
col <- c("Red","Blue","Yellow","Green","Brown","Pink","Violet","Brown","Black")
plot(T,rep(6000,length(T)), xlab="Time",type="n", ylab="Price", main ="Plot of Default Option Prices",ylim=c(1500,6000), xlim=c(3,9))
sapply(1:length(lamb2), function(k) lines(T,sapply(1:length(T), function(l) Proj6_2func(0.2,lamb2[k],T[l])$optionvalue), type="o", col=col[k]))
legend("bottomright", label,lty=rep(1,length(lamb2)),lwd=rep(2.5,length(lamb2)),col=col,cex=0.4)

label <- c("Lambda1=0.05","Lambda1=0.1","Lambda1=0.15","Lambda1=0.2","Lambda1=0.25","Lambda1=0.3","Lambda1=0.35","Lambda1=0.4")
col <- c("Red","Blue","Yellow","Green","Brown","Pink","Violet","Brown")
plot(T,rep(6000,length(T)), type="n", xlab="Time", ylab="Price", main ="Plot of Default Option Prices",ylim=c(1500,6000), xlim=c(3,9))
sapply(1:length(lamb1), function(k) lines(T,sapply(1:length(T), function(l) Proj6_2func(lamb1[k],0.4,T[l])$optionvalue), type="o", col=col[k]))
legend("bottomright", label,lty=rep(1,length(lamb1)),lwd=rep(2.5,length(lamb1)),col=col,cex=0.4)
#dev.off()

# plots for (b) part
par(mar=c(4,6,2,6))
par(mfrow=c(2,1))

label <- c("Lambda2=0","Lambda2=0.1","Lambda2=0.2","Lambda2=0.3","Lambda2=0.4","Lambda2=0.5","Lambda2=0.6","Lambda2=0.7","Lambda2=0.8")
col <- c("Red","Blue","Yellow","Green","Brown","Pink","Violet","Brown","Black")
plot(T,rep(1,length(T)), xlab="Time",type="n", ylab="Probability", main ="Plot of Default Probability",ylim=c(0,1),xlim=c(3,9))
sapply(1:length(lamb2), function(k) lines(T,sapply(1:length(T), function(l) Proj6_2func(0.2,lamb2[k],T[l])$defprob), type="o", col=col[k]))
legend("bottomright", label,lty=rep(1,length(lamb2)),lwd=rep(2.5,length(lamb2)),col=col,cex=0.4)

label <- c("Lambda1=0.05","Lambda1=0.1","Lambda1=0.15","Lambda1=0.2","Lambda1=0.25","Lambda1=0.3","Lambda1=0.35","Lambda1=0.4")
col <- c("Red","Blue","Yellow","Green","Brown","Pink","Violet","Brown")
plot(T,rep(1,length(T)), type="n", xlab="Time", ylab="Probability", main ="Plot of Default Probability",ylim=c(0,1),xlim=c(3,9))
sapply(1:length(lamb1), function(k) lines(T,sapply(1:length(T), function(l) Proj6_2func(lamb1[k],0.4,T[l])$defprob), type="o", col=col[k]))
legend("bottomright", label,lty=rep(1,length(lamb1)),lwd=rep(2.5,length(lamb1)),col=col,cex=0.4)
#dev.off()

# plots for (c) part
par(mar=c(4,6,2,6))
par(mfrow=c(2,1))

label <- c("Lambda2=0","Lambda2=0.1","Lambda2=0.2","Lambda2=0.3","Lambda2=0.4","Lambda2=0.5","Lambda2=0.6","Lambda2=0.7","Lambda2=0.8")
col <- c("Red","Blue","Yellow","Green","Brown","Pink","Violet","Brown","Black")
plot(T,rep(5,length(T)), xlab="Time",type="n", ylab="Exercise Time", main ="Plot of Expected Exercise Times",ylim=c(0,5),xlim=c(3,9))
sapply(1:length(lamb2), function(k) lines(T,sapply(1:length(T), function(l) Proj6_2func(0.2,lamb2[k],T[l])$exptime), type="o", col=col[k]))
legend("bottomright", label,lty=rep(1,length(lamb2)),lwd=rep(2.5,length(lamb2)),col=col,cex=0.4)

label <- c("Lambda1=0.05","Lambda1=0.1","Lambda1=0.15","Lambda1=0.2","Lambda1=0.25","Lambda1=0.3","Lambda1=0.35","Lambda1=0.4")
col <- c("Red","Blue","Yellow","Green","Brown","Pink","Violet","Brown")
plot(T,rep(5,length(T)), type="n", xlab="Time", ylab="Exercise Time", main ="Plot of Expected Exercise Times",ylim=c(0,5), xlim=c(3,9))
sapply(1:length(lamb1), function(k) lines(T,sapply(1:length(T), function(l) Proj6_2func(lamb1[k],0.4,T[l])$exptime), type="o", col=col[k]))
legend("bottomright", label,lty=rep(1,length(lamb1)),lwd=rep(2.5,length(lamb1)),col=col,cex=0.4)


# Final Output 

# Default option values
D1
D2
D3

# Default Probability values
Prob1
Prob2
Prob3

# Expected Exercise Time
Et1
Et2
Et3
