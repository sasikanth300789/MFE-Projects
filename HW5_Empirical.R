#Loading the packages
library(DataAnalytics)
library(zoo)
library(quantmod)
library(xts)
library(sandwich)
library(lmtest)

#Read Data
excelsheet<- read.csv("C:/Users/sasikanth/Desktop/MFE Courses/Winter 2016/237E/Fama_bond_prices.csv",header=TRUE)
excelsheet <- excelsheet/100
#Getting all dates in YYYY-MM_DD format
#sampleDate= as.character(excelsheet$qdate)
#sampleDate= c(as.Date(as.character(sampleDate), "%Y%m%d"))

#Get logPrice1,logPrice2,logPrice3,logPrice4,logPrice5 from the data
price1=log(excelsheet$price1)
price2=log(excelsheet$price2)
price3=log(excelsheet$price3)
price4=log(excelsheet$price4)
price5=log(excelsheet$price5)

#Question 1(a)
#Test null hypothesis B1=1 and B0=0
#Question 1(b)
ft1=as.numeric()
ft2=as.numeric()
ft3=as.numeric()
ft4=as.numeric()
x1=as.numeric()
x2=as.numeric()
x3=as.numeric()
x4=as.numeric()
y1=as.numeric()
y2=as.numeric()
y3=as.numeric()
y4=as.numeric()

#Define left side as y, right side as x
for (n in 1:4) {
  len <- nrow(excelsheet) - n*12
  for (i in 1:len)
  {
    if (n==1)
    {
      y1[i]=-price1[i+(n*12)]+price1[i]
      ft1[i]=price1[i] - price2[i]
      x1[i]= ft1[i]+price1[i]
    }
    if (n==2)
    {
      x2[i]=(-0.5*price2[i]+price1[i])/(n-1)
      ft2[i]=price2[i] - price3[i]
      y2[i]= -price1[i+12]+0.5*price2[i]
    } 
    if (n==3)
    {
      x3[i]=(-(1/n)*price3[i]+price1[i])/(n-1)
      ft3[i]=price2[i] - price3[i]
      y3[i]= -0.5*price2[i+12]+(1/n)*price3[i]
    }     
    else
    {
      x4[i]=(-(1/4)*price4[i]+price1[i])/(n-1)
      ft4[i]=price2[i] - price3[i]
      y4[i]= -(1/3)*price3[i+12]+(1/4)*price4[i]
    }         
  }
}

#Regression for n=1
plot(y1,x1,main = "Regression for n=1",xlab="Forward Rate-Yield",ylab ="Difference between Yields")
out1=lm(y1~x1,data=excelsheet)
summary(out1)
lmSumm(out1, HAC=TRUE)
#Regression for n=2
plot(y2,x2,main = "Regression for n=2",xlab="Forward Rate-Yield",ylab ="Difference between Yields")
plot(y1,x1,main = "Regression for n=1")
out2=lm(y2~x2,data=excelsheet)
summary(out2)
lmSumm(out2, HAC=TRUE)

#Regression for n=3
plot(y3,x3,main = "Regression for n=3",xlab="Forward Rate-Yield",ylab ="Difference between Yields")
out3=lm(y3~x3,data=excelsheet)
summary(out3)
lmSumm(out3, HAC=TRUE)

#Regression for n=4
plot(y4,x4,main = "Regression for n=4",xlab="Forward Rate-Yield",ylab ="Difference between Yields")
out4=lm(y4~x4,data=excelsheet)
summary(out4)
lmSumm(out, HAC=TRUE)


#Question 2(a)
#Test hypothesis gamma1=1, gamma0=0
#Question 2(b)
a1=as.numeric()
a2=as.numeric()
a3=as.numeric()
a4=as.numeric()

#Define left side as a, right side as x
for (n in 1:4)
{
  for (i in 1:nrow(excelsheet)-12)
  {
    if (n==1)
    {
      a1[i]=price1[i+12]-price2[i]+price1[i]
      ft1[i]=price1[i]+(n*-1/n)*price2[i]
      x1[i]= ft1[i]+price1[i]
    }
    if (n==2)
    {
      a2[i]=price1[i+12]-price3[i]+price1[i]
      ft2[i]=-price2[i]+(n*-1/n)*price3[i]
      x2[i]= ft2[i]-price1[i]
    } 
    if (n==3)
    {
      a3[i]=-price2[i+12]-price4[i]-price1[i]
      ft3[i]=-price3[i]+(n*-1/n)*price4[i]
      x3[i]= ft3[i]-price1[i]
    }     
    else
    {
      a4[i]=-price4[i+12]-price5[i]-price1[i]
      ft4[i]=-price4[i]+(n*-1/n)*price5[i]
      x4[i]= ft4[i]-price1[i]
    }
  }
}
#Regression for n=1
plot(a1,x1,main = "Regression for n=1",xlab="Forward Rate-Yield",ylab ="Holding Period Return-Yield")
out1=lm(a1~x1,data=excelsheet)
summary(out1)
coeftest(out1, vcov. = vcovHAC)


#Regression for n=2
plot(a2,x2,main = "Regression for n=2",xlab="Forward Rate-Yield",ylab ="Holding Period Return-Yield")
out2=lm(a2~x2,data=excelsheet)
summary(out2)
coeftest(out2, vcov. = vcovHAC)

#Regression for n=3
plot(a3,x3,main = "Regression for n=3",xlab="Forward Rate-Yield",ylab ="Holding Period Return-Yield")
out3=lm(a3~x3,data=excelsheet)
summary(out3)
coeftest(out3, vcov. = vcovHAC)

#Regression for n=4
plot(a4,x4,main = "Regression for n=4",xlab="Forward Rate-Yield",ylab ="Holding Period Return-Yield")
out4=lm(a4~x4,data=excelsheet)
summary(out4)
coeftest(out4, vcov. = vcovHAC)

