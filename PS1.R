
# Set Directory
setwd("/Users/sruthi/Desktop/DA")

# Load the library
library(DataAnalytics)
data(mvehicles)
cars = mvehicles[mvehicles$bodytype != "Truck",]

# Q1 Part A

#1.a Show the logical vector
vect <- c(cars$make == "Ford")
vect[1:50]

#1.b Create a two step selection
vect <- c(cars$make == "Ford")
ford <- cars[vect,]

#1.c number of kia
kia <- cars[c(cars$make == "Kia"),]
nrow(kia)

#1.d
emv <- cars[c(cars$emv > 100000),]
nrow(emv)


# Q1 Part B
europe <- cars[grepl("Europe",cars$origin,ignore.case=TRUE) & cars$emv >75000,]
avgsale <- mean(europe$sales)
avgsale

# Q1 Part C
fourdoor <- cars[grepl("4dr",cars$style,ignore.case = TRUE),]
nrow(fourdoor)

fourdrsedan <- cars[grepl("4dr Sedan",cars$style,ignore.case = TRUE),]
nrow(fourdrsedan)


# Question 2
# example for Inf
1/0
# example for -Inf
log(1/0)

# example for NaN
0*Inf
Inf - Inf

# missing values
vec <- c(rnorm(10),NA)
# mean of vec is missing because it has a missing value in it
mean(vec)
# To avoid it, we can use na.rm=TRUE to remove missing values
mean(vec,na.rm=TRUE)

# Question 3
final <- read.table("slots.csv",header=TRUE,sep=",")
slots <- read.delim("slots.txt",header=FALSE, sep=" ")
colnames(slots) <- colnames(final)

temp <- slots[,1:3]
temp[temp==1] <- "B"
temp[temp==2] <- "BB"
temp[temp==3] <- "BBB"
temp[temp==5] <- "DD"
temp[temp==6] <- "C"

slots[,1:3] <- temp

# Question 4

numeric = sapply(cars, is.numeric)
dfn = cars[, numeric, drop = FALSE]

#The line 20 is creating a logical vector for determining which of the columns in df are having numeric values and subsetting them out using line 22.
#Drop = false is used to preserve the original dimensions since we are subsetting inside a function.

# 4.b
descstat(df,0.95,3)

# Question 5
library(ggplot2)
data(diamonds)
cutf=as.character(diamonds$cut)
cutf=as.factor(cutf)

#data visualization
qplot(carat,price,data=diamonds,col=cut,size=I(2),alpha=I(3/4), main="Prices vs Carat vs Cut", cex=0.8) + theme_bw()

#separate plots
qplot(carat,price,data=diamonds,facets=cut~.,col=I("Green"), main="Price vs Carat plot for each cut") + theme_bw()

#boxplots
qplot(cut,price,data=diamonds,geom="boxplot",fill=I("Red")) + theme_bw()

#Using log(prices)
#data visualization
qplot(carat,log(price),data=diamonds,col=cut,size=I(2),alpha=I(3/4), main="Prices vs Carat vs Cut", cex=0.8) + theme_bw()

#separate plots
qplot(carat,log(price),data=diamonds,facets=cut~.,col=I("Green"), main="Price vs Carat plot for each cut") + theme_bw()

#boxplots
qplot(cut,log(price),data=diamonds,geom="boxplot",fill=I("Red")) + theme_bw()

#Using sqrt(prices)
#data visualization
qplot(carat,sqrt(price),data=diamonds,col=cut,size=I(2),alpha=I(3/4), main="Prices vs Carat vs Cut", cex=0.8) + theme_bw()

#separate plots
qplot(carat,sqrt(price),data=diamonds,facets=cut~.,col=I("Green"), main="Price vs Carat plot for each cut") + theme_bw()

#boxplots
qplot(cut,sqrt(price),data=diamonds,geom="boxplot",fill=I("Red")) + theme_bw()

# Question 5.b
# Trying out the Multiple Linear Regression model using square root of prices
library(car)
fit <- lm(log(price) ~ carat + cutf, data=diamonds)
summary(fit)

#fitted Values
fittedval <- fitted(fit)

#residual Value
resvalue <- residuals(fit)

#Residual Plots versus the predictor variables
par(mfrow=c(1,1))
plot(diamonds$carat,fit$residuals)

#Residual plot versus fitted values
plot(fittedval,resvalue)

# Assessing Outliers
layout(matrix(c(1),1,1))
outlierTest(fit) # Bonferonni p-value for most extreme obs
qqPlot(fit, main="QQ Plot") #qq plot for studentized resid 
leveragePlots(fit) # leverage plots

#plotting residuals versus fitted and carat
qplot(diamonds$carat,resvalue)
qplot(resvalue,fittedval)

# Normal Q-Q plots
qqPlot(fit, main="QQ Plot")
qqnorm(fit$residuals)
qqline(fit$residuals)

# Histogram
hist(fit$residuals)


# distribution of studentized residuals
library(MASS)
sresid <- studres(fit) 
hist(sresid, freq=FALSE, main="Distribution of Studentized Residuals")
xfit<-seq(min(sresid),max(sresid),length=40) 
yfit<-dnorm(xfit) 
lines(xfit, yfit)
