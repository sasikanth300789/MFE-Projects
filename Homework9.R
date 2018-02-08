# Project 9

setwd('~/desktop/237G')
rm(list = ls())

# Question 1
# Numerix Prepayment Model

# function for explicit CIR model
explicitCIR <- function(rt,rbar,sigma,kappa,t,T){
  
  h1 <- sqrt( kappa^2 + 2*sigma^2 )
  h2 <- 0.5*(kappa+h1)
  h3 <- (2*kappa*rbar)/(sigma^2)  
  B <- (exp(h1*(T-t)) -1)/(h2*(exp(h1*(T-t)) -1) + h1)
  A <-  (h1*exp(h2*(T-t))/(h2*(exp(h1*(T-t)) -1) + h1))^h3 
  price <- A*exp(-B*rt)
  return (price)
}

# fucntion for MBS price
MBSprice <- function(AMT,WAC,r0,rbar,kappa,sigma,T,x,model){

nmonths <- T*12
time <- seq(0,nmonths/12,1/12)
r <- WAC/12
R <- 12*r
N <- 360

# CIR Model
  t <- 0
  days <- floor((T-t)*360)
  paths <- 100
  dt <- 1/360
  
  set.seed(0)
  Z <- matrix(rnorm(paths*days,0,1),nrow=paths,ncol=days)
  Rsim <- cbind(rep(r0,paths),matrix(NA,paths,days))
  
  for (i in 1:paths){
    for (j in 2:(days+1)){
      Rsim[i,j] <- Rsim[i,j-1] + kappa*(rbar - Rsim[i,j-1] )*dt + sqrt(max(0,Rsim[i,j-1]))*sigma*sqrt(dt)*Z[i,j-1] 
    }
  }
  
  #parallel shift in all rates by x bps
  Rsim <- Rsim + x
  R <- -dt*(t(apply(Rsim[,-1],1,cumsum))[,seq(30,T*360,30)])
  
# explicit CIR prices for 10year  
CIRprice <- t(as.matrix(cbind(sapply(1:paths, function(i) c(sapply(0:(nmonths-1), function(j) explicitCIR(Rsim[i,(j*30 + 1)],rbar,sigma,kappa,time[j+1],time[j+1]+10))) ))))
longrate <- -(1/10)*(log(CIRprice))

# Discount Rates
D <- exp(R)

if (model == "Numerix"){
  
# Numerix model
NPY <- list( RI = matrix(NA,paths,nmonths), BU = matrix(NA,paths,nmonths), SG = matrix(NA,paths,nmonths), SY=rep(c(0.94,0.76,0.74,0.95,0.98,0.92,0.98,1.10,1.18,1.22,1.23,0.98),T))

# Mortgage List
mortgage <- list( TPP = matrix(0,paths,nmonths), IP = matrix(0,paths,nmonths), C = matrix(0,paths,nmonths), D=D, PV = cbind(rep(AMT,paths),matrix(0,paths,nmonths)), CPR = matrix(0,paths,nmonths) )

# run loop to each month
for (j in 1:paths){
for (i in 1:nmonths){
  
  NPY$BU[[j,i]] <- 0.3 + 0.7*(mortgage$PV[[j,i]])/AMT 
  NPY$SG[[j,i]] <- min(1,i/30)
  NPY$RI[[j,i]] <- 0.28 + 0.14*atan(-8.57 + 430*(12*r - longrate[j,i]) )
  
  mortgage$CPR[[j,i]] <- NPY$RI[[j,i]] * NPY$BU[[j,i]] * NPY$SG[[j,i]] * NPY$SY[[i]]
  mortgage$C[[j,i]] <- mortgage$PV[[j,i]]*r/(1- (1+r)^(i-1-N))  + (mortgage$PV[[j,i]] - mortgage$PV[[j,i]]*r*(1/(1- (1+r)^(i-1-N)) - 1))*(1-(1-mortgage$CPR[[j,i]])^(1/12)) 
  mortgage$IP[[j,i]] <- r*mortgage$PV[[j,i]]
  mortgage$TPP[[j,i]] <- mortgage$C[[j,i]] - mortgage$IP[[j,i]]
  
  # check if PV is becoming negative or not
  if (mortgage$TPP[[j,i]] > mortgage$PV[[j,i]]){
    mortgage$TPP[[j,i]] <- mortgage$PV[[j,i]]
  }
  
  mortgage$PV[[j,(i+1)]] <- mortgage$PV[[j,i]] - mortgage$TPP[[j,i]]
  
 }
}
} else {
  
  # PSA Model
  # Mortgage List
  mortgage <- list( TPP = matrix(0,paths,nmonths), IP = matrix(0,paths,nmonths), C = matrix(0,paths,nmonths), D=D, PV = cbind(rep(AMT,paths),matrix(0,paths,nmonths)), CPR = matrix(0,paths,nmonths) )
  
  for (j in 1:paths){
    for (i in 1:nmonths){
      
      mortgage$CPR[[j,i]] <- min(0.002*i,0.06)
      mortgage$C[[j,i]] <- mortgage$PV[[j,i]]*r/(1- (1+r)^(i-1-N))  + (mortgage$PV[[j,i]] - mortgage$PV[[j,i]]*r*( 1/(1-(1+r)^(i-1-N)) - 1))*(1-(1-mortgage$CPR[[j,i]])^(1/12)) 
      mortgage$IP[[j,i]] <- r*mortgage$PV[[j,i]]
      mortgage$TPP[[j,i]] <- mortgage$C[[j,i]] - mortgage$IP[[j,i]]
      
      # check if PV is becoming negative or not
      if (mortgage$TPP[[j,i]] > mortgage$PV[[j,i]]){
        mortgage$TPP[[j,i]] <- mortgage$PV[[j,i]]
      }
      
      mortgage$PV[[j,(i+1)]] <- mortgage$PV[[j,i]] - mortgage$TPP[[j,i]]
    
    }
  }  
}
price <- mean(rowSums(mortgage$C*mortgage$D))
ioprice <- mean(rowSums(mortgage$IP*mortgage$D))
poprice <- mean(rowSums(mortgage$TPP*mortgage$D))
return( list(price = price, ioprice = ioprice, poprice = poprice) )
}

###### End of function code #######

# Question 1.a

#Input from user
T <- 30
AMT <- 100000
WAC <- 0.08
r0 <- 0.078
kappa <- 0.6
rbar <- 0.08
sigma <- 0.12
model <- "Numerix"

Numerix_model_price <- MBSprice(AMT,WAC,r0,rbar,kappa,sigma,T,0,model)$price

#1.b
K <- seq(0.3,0.9,0.1)
mbsprice <- c(sapply(1:length(K), function(i) MBSprice(AMT,WAC,r0,rbar,K[i],sigma,T,0,model)$price))
plot(K,mbsprice, xlab = "Kappa", ylab = "MBS prices", main = "Plot of MBS vs Kappa", type="o", col="red") 

#1.c
Rbar <- seq(0.03,0.09,0.01)
mbsprice2 <- c(sapply(1:length(Rbar), function(i) MBSprice(AMT,WAC,r0,Rbar[i],kappa,sigma,T,0,model)$price))
plot(Rbar,mbsprice2, xlab = "Rbar", ylab = "MBS prices", main = "Plot of MBS vs Rbar", type="o", col="red") 


# Question 2
#Input from user
T <- 30
AMT <- 100000
WAC <- 0.08
r0 <- 0.078
kappa <- 0.6
rbar <- 0.08
sigma <- 0.12
model <- "PSA"

PSA_model_price <- MBSprice(AMT,WAC,r0,rbar,kappa,sigma,T,0,model)$price

# Question 2.a
K <- seq(0.3,0.9,0.1)
mbsprice3 <- c(sapply(1:length(K), function(i) MBSprice(AMT,WAC,r0,rbar,K[i],sigma,T,0,"PSA")$price))
plot(K,mbsprice3, xlab = "Kappa", ylab = "MBS prices (PSA)", main = "Plot of MBS vs Kappa", type="o", col="red") 

# Question 3
f <- function(x)  abs(110000 - MBSprice(AMT,WAC,r0,rbar,kappa,sigma,T,x,"Numerix")$price)
OAS <- optimize(f,interval = c(-0.05,0.05))$minimum

#Question 4
x <- OAS
y <- 0.0005

p0 <- MBSprice(AMT,WAC,r0,rbar,kappa,sigma,T,x,"Numerix")$price
pplus <- MBSprice(AMT,WAC,r0,rbar,kappa,sigma,T,x+y,"Numerix")$price
pminus  <- MBSprice(AMT,WAC,r0,rbar,kappa,sigma,T,x-y,"Numerix")$price

OASDuration <- (pminus-pplus)/(2*y*p0)
OASConvex <- (pplus + pminus - 2*p0)/(2*p0*y^2)

# Question 5

IOprice <- c(sapply(1:length(Rbar), function(i) MBSprice(AMT,WAC,r0,Rbar[i],kappa,sigma,T,0,"Numerix")$ioprice ))
POprice <- c(sapply(1:length(Rbar), function(i) MBSprice(AMT,WAC,r0,Rbar[i],kappa,sigma,T,0,"Numerix")$poprice ))

############################### Thank You ###########################################