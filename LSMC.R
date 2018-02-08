# remove all variables
rm(list=ls(all=TRUE))
  
# Question 1
# Input Values
S0 = c(36,40,44) # Stock price
K = 40 # Strike price
r = 0.06 # Risk free rate
sigma = 0.2 # Volatility
n = 100000 # number of paths
m = 50 # time step
T = c(0.5,1,2) # time to expiration

# function used for calculation Expected Continuation Value
getECV <- function(X,Y,npoly,type){
  
  if (type=="L"){
    L1 <- exp(-X/2)
    L2 <- (1-X)*exp(-X/2)
    L3 <- exp(-X/2) * (1-2*X + 0.5*X*X)
    L4 <- exp(-X/2) * (1 - 3*X + 1.5*X*X - X*X*X/6)
  } else if(type == "H") {
    L1 <- 1
    L2 <- 2*X
    L3 <- 4*X*X - 2
    L4 <- 8*X*X*X - 12*X
  } else {
    L1 <- 1
    L2 <- X
    L3 <- X*X
    L4 <- X*X*X
  }
  
  L <- cbind(L1,L2,L3,L4)[,1:npoly]
  #A <- crossprod(L,L)
  #b <- crossprod(L,Y)
  coef <- as.numeric(lm(Y~0+L)$coef) #Linear Regression used to avoid matrix singularity issues
  return (L%*%coef)
}

# Function for Least Squares Monte Carlo
LSMC <- function(S0,K,r,sigma,n,m,T,Z,npoly,type){

# time steps
dt = seq(0,T,T/m)

# Stock price simulation
# simulate stock paths using antithetic variates
normal <- matrix(exp((r-0.5*sigma^2)*(T/m) + sigma*sqrt(T/m)*Z),nrow=n/2,ncol=m)
antithetic <- 0.5*(normal + matrix(exp((r-0.5*sigma^2)*(T/m) + sigma*sqrt(T/m)*Z),nrow=n/2,ncol=m))                
pr_norm <- cbind(rep(1,n),rbind(normal,antithetic))
pr <- S0*t(apply(pr_norm,1,cumprod))

# create discount factor and exercise value matrices
cv <- sapply(1:(m+1), function(k) sapply(1:n, function(j) max(0,K-pr[j,k])))
disc <- t(matrix(rep(sapply(1:(m+1), function(k) exp(-r*dt[k])),n),m+1,n))
prices <- list(pr=pr,cv=cv,ecv = matrix(NA,n,m+1), Ind = matrix(NA,n,m+1), disc = disc)

# create ECV and Index matrix
prices$ecv[,m+1] <- 0
prices$Ind[,m+1] <- sapply(1:n, function(i) if(prices$cv[i,m+1] > 0) {1} else {0}) 

# Loop for calculating the Index matrix
for (j in m:2){
  
X <- prices$pr[,j]
Y <- rowSums(sapply((j+1):(m+1), function(i) exp(-r*dt[i-j])*prices$Ind[,i]*prices$cv[,i] ))
prices$ecv[,j] <- getECV(X,Y,npoly,type)
prices$Ind[,j] <- sapply(1:n, function(i) if(prices$cv[i,j] >= prices$ecv[i,j] & prices$cv[i,j] > 0) {1} else {0}) 
prices$Ind[which(prices$Ind[,j]==1),j:m+1] <- 0
}

# mean value of american put option prices from stock path simulations
V0 <- mean(rowSums(prices$disc[,-1]*prices$Ind[,-1]*prices$cv[,-1]))
return (V0)
}

#main code
# Generate random numbers of length n/2 because the antithetic will have n/2 numbers as well.
set.seed(0)
Z <- matrix(rnorm(n*m/2),n/2,m) 

# List of Results
final <- list(Laguerre=c(),Hermite=c(),Simple=c())

final$Laguerre[[1]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,2,"L"))))
final$Laguerre[[2]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,3,"L"))))
final$Laguerre[[3]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,4,"L"))))

final$Hermite[[1]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,2,"H"))))
final$Hermite[[2]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,3,"H"))))
final$Hermite[[3]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,4,"H"))))

final$Simple[[1]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,2,"S"))))
final$Simple[[2]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,3,"S"))))
final$Simple[[3]] <- cbind(sapply(1:length(T),function(t) sapply(1:length(S0), function(s) LSMC(S0[s],K,r,sigma,n,m,T[t],Z,4,"S"))))

# Output the results
#Laguerre polynomials
#k = 2
final$Laguerre[[1]]

#Laguerre polynomials
#k = 3
final$Laguerre[[2]]

#Laguerre polynomials
#k = 4
final$Laguerre[[3]]

#Hermite polynomials
#k = 2
final$Hermite[[1]]

#Hermite polynomials
#k = 3
final$Hermite[[2]]

#Hermite polynomials
#k = 4
final$Hermite[[3]]

#Simple monomials
#k = 2
final$Small[[1]]

#Simple monomials
#k = 3
final$Simple[[2]]

#Simple monomials
#k = 4
final$Simple[[3]]



# Question 2
# Forward Start options
S0 = 65 # Stock price
K = 60 # Stike price
r = 0.06 # Risk Free
sigma = 0.2 # Volatility
n = 10000 # number of paths
m = 100 # time step
t = 0.2 # Start time
T = 1 # time to expiry
dt = seq(0,T,T/m) # time steps

# Generate random normal numbers
set.seed(0)
Z = matrix(rnorm(n*m),n,m) 

# Simulate stock paths
normal <- matrix(exp((r-0.5*sigma^2)*(T/m) + sigma*sqrt(T/m)*Z),nrow=n,ncol=m)
pr_norm <- cbind(rep(1,n),normal)
pr <- S0*t(apply(pr_norm,1,cumprod))

# Create exercise value and discount factor matrix
cv <- sapply(1:(m+1), function(k) sapply(1:n, function(j) max(0,pr[j,1+(m*0.2)]-pr[j,k])))
disc <- t(matrix(rep(sapply(1:(m+1), function(k) exp(-r*dt[k])),n),m+1,n))

# European put price
Euroean_put <- exp(-r*T)*mean(cv[,m+1])
European_put

# 2.b

# We first simulate the strike prices and then simulate stock paths for each strike starting at T =0.2
# Simulate strikes
n <- 1000 
m <- 50
T <- 0.8 #option time
t <- 0.2 # Start time

# generate random number for t = 0.
set.seed(0)
strike <- S0*exp((r-0.5*sigma^2)*(t) + sigma*sqrt(t)*rnorm(n/2))
putprice <- array(NA,length(strike))

#loop for calucalting put price for each strike by simulating stock paths
for (i in 1:length(strike)){
  Z <- matrix(rnorm(n*m/2),n/2,m) 
  putprice[i] <- LSMC(strike[i],strike[i],r,sigma,n,m,T,Z,3,"H")
}

# Put price
American_put <- exp(-t*r)*mean(putprice)
American_put