#Question 1: Random Number Generator 

#LGM Method
options(scipen=10)
n <- 10000
a <- 7^5
b <- 0
m <- 2^31 - 1

lgm <- function(n,a,m){
 
rand <- numeric(n)
rand[1] <- 0.001
for (i in 2:n) {
    rand[i] = (rand[i-1]*a)%%m
}
return (rand/m) 
}

unirand <- lgm(n,a,m)

mean(unirand)
sd(unirand)

srand <- runif(n)
mean(srand)
sd(srand)


## Question 2: Distribution function
library(ggplot2)
X <- numeric(n)
p <- c(0.3,0.35,0.2,0.15)
for (j in 1:n){
if (unirand[j] <= p[1]){
    X[j] = -1  
} 
else if (unirand[j] > p[1] && unirand[j] <= sum(p[1:2])){
    X[j] = 0 
}
else if (unirand[j] > sum(p[1:2]) && unirand[j] <= sum(p[1:3])){
    X[j] = 1 
}
else{ 
    X[j] = 2 
}
}
X[1:10]

qplot(X,geom="histogram", fill = I("blue"), col = I("red"), binwidth = 0.1, main = "Frequency distribution of X")
mean(X)
sd(X)

# Question 4: exponentially distributed

X <- -1.5*log(1-unirand)
F[1] <- length(which(X>=1))/length(X)
F[4] <- length(which(X>=4))/length(X)
mean(X)
sd(X)
qplot(X,geom="histogram", fill = I("blue"), binwidth = 0.1, main = "Frequency distribution of X")


# Question 3: Binomial Distribution
M <- 2^31 -1
m <- 44
p <- 0.64
uniform <- as.numeric(lgm(n*m,a,M),n*m)
vec <- which(uniform<p) 
uniform[vec] <- 1
uniform[-vec] <- 0
uniformmatrix <- matrix(uniform,n,m)
binomial <- c(rowSums(uniformmatrix))   

qplot(binomial,geom="histogram", fill = I("blue"), col = I("red"), binwidth = 0.1, main = "Frequency distribution of X")
prob <- length(which(binomial>=40))/length(binomial)
bprob <- sum(dbinom(40:44,44,0.64))

# Question 5: Normal Distribution
unirand <- runif(5000)

#box-muller method
start <- Sys.time()
unimatrix <- matrix(runif(5000),2500,2)
Zone <- sqrt(-2*log(unimatrix[,1]))*cos(2*pi*unimatrix[,2])
Ztwo <- sqrt(-2*log(unimatrix[,1]))*sin(2*pi*unimatrix[,2])
Z <- rbind(Zone, Ztwo)
end <- Sys.time()
mean(Z)
sd(Z)
BMtime <- end - start

#Polar-Marsaglia
start <- Sys.time()
Urand <- matrix(runif(6500),ncol=2)
Vone <- 2*Urand[,1] - 1
Vtwo <- 2*Urand[,2] - 1
w <- Vone^2 + Vtwo^2
Urand <- cbind(Urand,Vone,Vtwo,w)
Urand <- subset(Urand, w<=1)[c(1:2500),]
polarZone <- Urand[,3] * sqrt((-2*log(Urand[,5]))/Urand[,5])
polarZtwo <- Urand[,4] * sqrt((-2*log(Urand[,5]))/Urand[,5])
polarZ <- rbind(polarZone,polarZtwo)
end <- Sys.time()
mean(polarZ)
sd(polarZ)
PMtime <- end-start
#Efficiency