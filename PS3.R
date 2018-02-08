# Set Directory
setwd("/Users/sruthi/Desktop/DA")

rm(list=ls())
# we need the foreign package to import data in different format
require(foreign)
require(data.table)
require(ggplot2)
require(dplyr)

# Download data and set as data.table
loan <- as.data.table(read.dta("LendingClub_LoanStats3a_v12.dta"))

# Question a
# Selecting rows with "Fully Paid" or "Charged Off"
loandata <- loan[loan$loan_status %in% c("Fully Paid","Charged Off")]

# Define new variable Default
loandata$default <- 0
loandata$default[c(loandata$loan_status %in% "Charged Off")] <- 1

# Average Default
avg_default <- sum(loandata$default)/length(loandata$default)

# Question b
# fit logit model
out=glm(default ~ grade,family=binomial(link = "logit"),data=loandata)
summary(out) 

# Null Deviance test
test_stat <- out$null.deviance - out$deviance
k <- out$df.null - out$df.residual
pvalue_chisq <- 1-pchisq(test_stat,df=k)
pvalue_chisq

#Lift table
phat <- predict(out,type="response")
deciles <- ntile(phat,10)
#deciles <- cut(phat,breaks=quantile(phat,probs=c(seq(0,1,0.1))),include.lowest=TRUE)
deciles <- as.numeric(deciles)
df <- data.frame(deciles=deciles,phat=phat,default=loandata$default)
lift <- aggregate(df,by=list(deciles),FUN="mean",data=df) # find mean default for each decile
lift <- lift[,c(2,4)]
lift[,3] <- lift[,2]/mean(loandata$default)
names(lift)=c("decile","Mean Response","Lift Factor")
lift

# Create ROC curve
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

# Using a full model
glm_simple_roc <- simple_roc(loandata$default=="1", phat)
TPR <- glm_simple_roc$TPR
FPR <- glm_simple_roc$FPR
glm_simple_roc <- cbind(glm_simple_roc, Model = "grade")

# plot the corresponding ROC curve
q <- qplot(FPR,TPR,xlab="FPR",ylab="TPR",col=I("blue"),main="ROC Curve for Logistic Regression Default Model",size=I(0.75))
# add straight 45 degree line from 0 to 1
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0)) + theme_bw()


# Question c
# fit new logit model
out_nongrade=glm(default ~ loan_amnt + annual_inc,family=binomial(link = "logit"),data=loandata)
summary(out_nongrade)

# lift tables
phat_nongrade <- predict(out_nongrade,type="response")
deciles <- ntile(phat_nongrade,10)
deciles <- as.numeric(deciles)
df <- data.frame(deciles=deciles,phat=phat_nongrade,default=loandata$default)
lift_nongrade <- aggregate(df,by=list(deciles),FUN="mean",data=df) # find mean default for each decile
lift_nongrade <- lift_nongrade[,c(2,4)]
lift_nongrade[,3] <- lift_nongrade[,2]/mean(loandata$default)
names(lift_nongrade)=c("decile","Mean Response NG","Lift Factor NG")
lift_nongrade

# Compare lift tables
lift_compare <- cbind(lift,lift_nongrade[,2:3])
lift_compare

#plot the ROC
phat_full_nongrade <- predict(out_nongrade,type="response")
glm_simple_roc_full_nongrade <- simple_roc(loandata$default=="1", phat_full_nongrade)
glm_simple_roc_full_nongrade <- cbind(glm_simple_roc_full_nongrade, Model = "non-grade")
new_roc <- rbind(glm_simple_roc, glm_simple_roc_full_nongrade)
q <- qplot(FPR,TPR,data = new_roc, colour = Model, xlab="FPR",ylab="TPR",main="ROC Curve for grade and non-grade Logistic Regression Models",size=I(0.75))
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0), col = I("black")) + theme_bw()

# Question c.2
out_term_int=glm(default ~ loan_amnt + annual_inc + term + int_rate,family=binomial(link = "logit"),data=loandata)
summary(out_term_int)

# lift tables
phat_term_int <- predict(out_term_int,type="response")
deciles <- ntile(phat_term_int,10)
deciles <- as.numeric(deciles)
df <- data.frame(deciles=deciles,phat=phat_term_int,default=loandata$default)
lift_term_int <- aggregate(df,by=list(deciles),FUN="mean",data=df) # find mean default for each decile
lift_term_int <- lift_term_int[,c(2,4)]
lift_term_int[,3] <- lift_term_int[,2]/mean(loandata$default)
names(lift_term_int)=c("decile","Mean Response NG","Lift Factor NG")
lift_term_int

# Compare lift tables
lift_compare_new <- cbind(lift,lift_term_int[,2:3])
lift_compare_new

#plot the ROC
phat_full_term_int <- predict(out_term_int,type="response")
glm_simple_roc_full_term_int <- simple_roc(loandata$default=="1", phat_full_term_int)
glm_simple_roc_full_term_int <- cbind(glm_simple_roc_full_term_int, Model = "term & int")
new_roc_term_int <- rbind(glm_simple_roc, glm_simple_roc_full_term_int)
q <- qplot(FPR,TPR,data = new_roc_term_int, colour = Model, xlab="FPR",ylab="TPR",main="ROC Curve for grade and term_int Logistic Regression Models",size=I(0.75))
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0), col = I("black")) + theme_bw()


#Question c.3
loandata$ratesquare <- loandata$int_rate^2
out_term_int_new <- glm(default ~ loan_amnt + annual_inc+ term + int_rate + ratesquare,family=binomial(link = "logit"),data=loandata)
summary(out_term_int_new)

phat_term_int_new <- predict(out_term_int_new,type="response")
