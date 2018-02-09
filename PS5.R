
# Problem Set 5

# Normal Copula Model Function
Copula <- function(Num,Size,PD,LGD,rho,alpha,N){
  
iterloss <- array(NA,dim = N)

for (i in 1:N){
  
  F <- matrix(rnorm(1,0,1),1,Num) 
  Z <- matrix(rnorm(Num,0,1),1,Num)
  U <- sqrt(rho)*F + sqrt(1-rho)*Z
  Default <- (U<qnorm(PD,0,1))
  loss <- LGD*Default*Size
  iterloss[i] <- sum(loss)
}

VaR <- quantile(iterloss,alpha)
return (VaR)
}

#Question 1
#1.a

Num <- 1000
Size <- 10000
PD <- 0.0019
rho <- 0.15
alpha <- 0.999
N <- 10000
LGD <- 0.45

WCDR <- pnorm((qnorm(PD) + sqrt(rho)*qnorm(alpha))/sqrt(1-rho),0,1)
closedformVaR <- Num*Size*LGD*(WCDR - PD)

#1.b
SimVaR <- Copula(Num,Size,PD,LGD,rho,alpha,N)

#1.c
newSize <- 100000
newNum <- 100
newSimVaR <- Copula(newNum,newSize,PD,LGD,rho,alpha,N)

#1.d
newSize <- 1000000
newNum <- 10
newSimVaR2 <- Copula(newNum,newSize,PD,LGD,rho,alpha,N)


# Question 2
# loans are paying coupon every year end
# default at the end of the year

libor <- 0.005
coupon <- Size*(libor + 0.021)

iterloss <- array(NA,dim = N)
  
for (i in 1:N){
    
    F <- matrix(rnorm(1,0,1),1,Num) 
    Z <- matrix(rnorm(Num,0,1),1,Num)
    U <- sqrt(rho)*F + sqrt(1-rho)*Z
    Default <- (U<qnorm(PD,0,1))
    loss <- Default*(LGD*Size+coupon)
    iterloss[i] <- sum(loss)
}
  
creditVaR <- quantile(iterloss,alpha) - Size*PD*LGD*Num
  
# Question 2.b
if (!require("FSA")) install.packages("FSA")
require("FSA")

ratings <- c("AAA","AA","A","BBB","BB","B","CCC","Default")
spread <- c(0.7,0.88,1.19,2.10,3.39,4.56,8.17,NA)
prob <- c(0.05,0.19,4.79,89.41,4.35,0.82,0.2,0.19) 
cumprob <- rcumsum(prob)/100
qtile <- qnorm(cumprob,0,1)

# Intial Market Value of Loan
portfolio <- Num*(coupon/(1+libor+0.021) + (Size+coupon)/(1+libor+0.021)^2)
newport <- array(NA,dim = N)

for (i in 1:N){
  
  F <- matrix(rnorm(1,0,1),1,Num) 
  Z <- matrix(rnorm(Num,0,1),1,Num)
  U <- sqrt(rho)*F + sqrt(1-rho)*Z
  newspread <-  c(sapply(1:Num,function(i) spread[which(qtile < U[i])[1]-1]/100))

# End of 1st year MTM
newport[i] <- sum(sapply(1:Num, function(k) ifelse(is.na(newspread[k]),Size*(1-LGD), coupon + ((coupon+Size)/(1+libor+newspread[k])) ) ))
}

MHV <- mean(newport)
XWC <- quantile(newport,1-alpha)
UL <- MHV-XWC


