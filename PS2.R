
# Set Directory
setwd("/Users/sruthi/Desktop/DA")

# we need the foreign package to import data in different format
require(foreign)
require(data.table)
require(ggplot2)
require(dplyr)

# Download data and set as data.table
StockRetAcct_DT <- as.data.table(read.dta("StockRetAcct_insample.dta"))

# set keys for StockRetAcct_DT In particular, it will be useful to sort on FirmID and year
setkey(StockRetAcct_DT, FirmID, year)

# create excess returns (what we really care about)
StockRetAcct_DT[,ExRet:=exp(lnAnnRet) - exp(lnRf)]

# remove NA from LnIssue
StockRetAcct_DT <- StockRetAcct_DT[!is.na(StockRetAcct_DT$lnIssue)]

# add some noise before breaking into quantiles
#StockRetAcct_DT[,lnIssue:=jitter(lnIssue, amount = 0)]

# loop through the years in the data base
for (i in 1983:2014)
{
  StockRetAcct_DT[year == i,quantile_yr:=ntile(StockRetAcct_DT[year == i,]$lnIssue,10)]
  #StockRetAcct_DT[year == i,quantile_yr:=cut(StockRetAcct_DT[year == i,]$lnIssue,breaks=quantile(StockRetAcct_DT[year == i,]$lnIssue,probs=c(0:10)/10,na.rm=TRUE, include.lowest=TRUE), labels=FALSE)]
}

# value weighted returns across portfolios
VW_MutualFunds_yr <- StockRetAcct_DT[,list(MeanExRetYr = weighted.mean(ExRet, MEwt)), by = list(quantile_yr, year) ]

# then average across years
VW_MutualFunds_yr <- VW_MutualFunds_yr[!is.na(VW_MutualFunds_yr$quantile_yr)]
VW_MutualFunds_yr <- VW_MutualFunds_yr[,list(MeanExRet = mean(MeanExRetYr)), by = quantile_yr]
VW_MutualFunds_yr[order(quantile_yr)]

qplot(quantile_yr, MeanExRet,data = VW_MutualFunds_yr, col=I("blue"), na.rm = TRUE, main = "VW Large Firm Issuance bins vs. Excess Returns") + geom_smooth(col=I("red")) + theme_bw()

#create Transformed Issuance 
StockRetAcct_DT <- StockRetAcct_DT[,transIssue:= 0]
StockRetAcct_DT$transIssue[which(StockRetAcct_DT$quantile_yr == 1)] <- -1
StockRetAcct_DT$transIssue[which(StockRetAcct_DT$quantile_yr == 10)] <- 1

# fama-mcbeth regression
port_ret = NULL

for (i in 1983:2014)
{
  temp <- StockRetAcct_DT[year==i,]
  fit_yr <- lm(temp$ExRet ~ temp$transIssue, data=temp)
  temp <- coefficients(fit_yr)
  port_ret = rbind(port_ret,temp[2])
}

fm_output = list(MeanReturn = mean(port_ret), StdReturn = sqrt(var(port_ret)), SR_Return = mean(port_ret)/sqrt(var(port_ret)), tstat_MeanRet = sqrt(1+2014-1983)*mean(port_ret)/sqrt(var(port_ret)))
fm_output



# Question 2
StockRetAcct_DT <- as.data.table(read.dta("StockRetAcct_insample.dta"))

# set keys for StockRetAcct_DT In particular, it will be useful to sort on FirmID and year
setkey(StockRetAcct_DT, FirmID, year)

# create excess returns (what we really care about)
StockRetAcct_DT[,ExRet:=exp(lnAnnRet) - exp(lnRf)]

# Remove NAs
StockRetAcct_DT <- StockRetAcct_DT[!is.na(StockRetAcct_DT$lnBM)]
StockRetAcct_DT <- StockRetAcct_DT[!is.na(StockRetAcct_DT$lnME)]

# add jitter
#StockRetAcct_DT[,lnBM:=jitter(lnBM, amount = 0.001)]

# We create the bins on the current year to calculate returns for next year. Loop through the years in the data base
for (i in 1981:2014)
{
  #StockRetAcct_DT[year == i,bm_quintile_yr:=cut(StockRetAcct_DT[year == i,]$lnBM,breaks=quantile(StockRetAcct_DT[year == i,]$lnBM,probs=c(0:5)/5,na.rm=TRUE, include.lowest=TRUE), labels=FALSE)]
  StockRetAcct_DT[year == i,bm_quintile_yr:=ntile(StockRetAcct_DT[year == i,]$lnBM,5)]
  StockRetAcct_DT[year == i,size_quintile_yr:=ntile(StockRetAcct_DT[year == i,]$lnME,5)]
  #StockRetAcct_DT[year == i,size_quintile_yr:=cut(StockRetAcct_DT[year == i,]$lnME,breaks=quantile(StockRetAcct_DT[year == i,]$lnME,probs=c(0:5)/5,na.rm=TRUE, include.lowest=TRUE), labels=FALSE)]
}

#StockRetAcct_DT[,lnME:=jitter(lnME, amount = 0.001)]

for (i in 1981:2014)
{
  #StockRetAcct_DT[year == i,bm_quintile_yr:=cut(StockRetAcct_DT[year == i,]$lnBM,breaks=quantile(StockRetAcct_DT[year == i,]$lnBM,probs=c(0:5)/5,na.rm=TRUE, include.lowest=TRUE), labels=FALSE)]
  #StockRetAcct_DT[year == i,size_quintile_yr:=cut(StockRetAcct_DT[year == i,]$lnME,breaks=quantile(StockRetAcct_DT[year == i,]$lnME,probs=c(0:5)/5,na.rm=TRUE, include.lowest=TRUE), labels=FALSE)]
}

#Using Value weighted returns across the firms in a year
StockRetAcct_DT <- StockRetAcct_DT[!is.na(StockRetAcct_DT$bm_quintile_yr)]
VW_MutualFunds_yr <- StockRetAcct_DT[,list(MeanExRetYr = weighted.mean(ExRet, MEwt)), by = list(bm_quintile_yr, size_quintile_yr, year)]

# then average returns across years
VW_MutualFunds_yr <- VW_MutualFunds_yr[,list(MeanExRet = mean(MeanExRetYr)), by = list(bm_quintile_yr, size_quintile_yr)]
VW_MutualFunds_yr[order(size_quintile_yr,bm_quintile_yr),]

#Using equal weighted returns across the firms in a year
EW_MutualFunds_yr <- StockRetAcct_DT[,list(MeanExRetYr = mean(ExRet)), by = list(bm_quintile_yr, size_quintile_yr, year)]

# then average returns across years
EW_MutualFunds_yr <- EW_MutualFunds_yr[,list(MeanExRet = mean(MeanExRetYr)), by = list(bm_quintile_yr, size_quintile_yr)]
EW_MutualFunds_yr[order(size_quintile_yr,bm_quintile_yr),]

#qplots for value weighted portfolios
qplot(bm_quintile_yr, MeanExRet,data = VW_MutualFunds_yr[size_quintile_yr == 1], col=I("blue"), na.rm = TRUE, main = "VW BM bins vs. Excess Returns for size quintile 1") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = VW_MutualFunds_yr[size_quintile_yr == 2], col=I("blue"), na.rm = TRUE, main = "VW BM bins vs. Excess Returns for size quintile 2") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = VW_MutualFunds_yr[size_quintile_yr == 3], col=I("blue"), na.rm = TRUE, main = "VW BM bins vs. Excess Returns for size quintile 3") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = VW_MutualFunds_yr[size_quintile_yr == 4], col=I("blue"), na.rm = TRUE, main = "VW BM bins vs. Excess Returns for size quintile 4") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = VW_MutualFunds_yr[size_quintile_yr == 5], col=I("blue"), na.rm = TRUE, main = "VW BM bins vs. Excess Returns for size quintile 5") + geom_smooth(col=I("red")) + theme_bw()

#qplots for equal weighted portfolios
qplot(bm_quintile_yr, MeanExRet,data = EW_MutualFunds_yr[size_quintile_yr == 1], col=I("blue"), na.rm = TRUE, main = "EW BM bins vs. Excess Returns for size quintile 1") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = EW_MutualFunds_yr[size_quintile_yr == 2], col=I("blue"), na.rm = TRUE, main = "EW BM bins vs. Excess Returns for size quintile 2") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = EW_MutualFunds_yr[size_quintile_yr == 3], col=I("blue"), na.rm = TRUE, main = "EW BM bins vs. Excess Returns for size quintile 3") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = EW_MutualFunds_yr[size_quintile_yr == 4], col=I("blue"), na.rm = TRUE, main = "EW BM bins vs. Excess Returns for size quintile 4") + geom_smooth(col=I("red")) + theme_bw()
qplot(bm_quintile_yr, MeanExRet,data = EW_MutualFunds_yr[size_quintile_yr == 5], col=I("blue"), na.rm = TRUE, main = "EW BM bins vs. Excess Returns for size quintile 5") + geom_smooth(col=I("red")) + theme_bw()
