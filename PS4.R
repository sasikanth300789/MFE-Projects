# Set Directory
setwd("/Users/sruthi/Desktop/DA")

rm(list = ls())
# clean-up R script for Problem Set 4

require(tm)
require(glmnet)
require(ggplot2)
require(pROC)
require(SnowballC)

# Read in the data, use the right file path for you in the below
data <- read.csv('DJIA_Headline_News.csv', stringsAsFactors = FALSE)

# First, we will clean up the data and do some quick preprocessing. Let's also add a '<\\s>' token between
# the headlines. We don't want the last word of a headline and the first word of the next to be counted as a bigram.

# Make 'Date' column a Date object to make train/test splitting easier
data$Date <- as.Date(data$Date)

# Combine headlines into one text blob for each day and add sentence separation token
data$all <- paste(data$Top1, data$Top2, data$Top3, data$Top4, data$Top5, data$Top6,
                  data$Top7, data$Top8, data$Top9, data$Top10, data$Top11, data$Top12, 
                  data$Top13, data$Top14, data$Top15, data$Top16, data$Top17, data$Top18,
                  data$Top19, data$Top20, data$Top21, data$Top22, data$Top23, data$Top24,
                  data$Top25, sep=' <s> ')

# Get rid of those pesky b's and backslashes you see if you inspect the raw data
data$all <- gsub('b"|b\'|\\\\|\\"', "", data$all)

# Get rid of all punctuation except headline separators, alternative to cleaning done in tm-package
data$all <- gsub("([<>])|[[:punct:]]", "\\1", data$all)

# Reduce to only the three columns we need. 
data <- data[, c('Date', 'Label', 'all')]

# Part 1
main_corpus <- Corpus(VectorSource(data$all))
#summary(main_corpus)

# Part 2

# remove punctuation
main_corpus <- tm_map(main_corpus, removePunctuation)
#inspect(main_corpus)

# remove special characters
for (jj in seq(main_corpus))
{
  main_corpus[[jj]] <- gsub("/", " ", main_corpus[[jj]])
  main_corpus[[jj]] <- gsub("@", " ", main_corpus[[jj]])
  main_corpus[[jj]] <- gsub("\\|", " ", main_corpus[[jj]])
  main_corpus[[jj]] <- gsub("$", " ", main_corpus[[jj]])
  main_corpus[[jj]] <- gsub("\t", " ", main_corpus[[jj]])
}

# remove numbers (we have those better represented in CompuStat)
main_corpus <- tm_map(main_corpus, removeNumbers)

# convert all to lowercase, so word is recognized with arbitrary capitalization
main_corpus <- tm_map(main_corpus, tolower)

# remove "stopwords" (e.g., and, to, a, as, the, ...)
main_corpus <- tm_map(main_corpus, removeWords, stopwords("english"))
main_corpus <- tm_map(main_corpus, PlainTextDocument)
# stemming words, i.e., keep only the stem so as not to differentially count investing, invest, invests
# taking out common word endings such as 'ing', 'es', and 's'
require(SnowballC)
main_corpus <- tm_map(main_corpus, stemDocument)

#part 3
# next, organize words into matrix, which then can be used for analysis
dtm <- DocumentTermMatrix(main_corpus)
inspect(dtm[5:10,801:810])

#part 4
# organize words by frequency
freq <- colSums(as.matrix(dtm))
ord_corpus <- order(freq)

#  identify words that appear frequently
# create a convenient data.frame
word_freq <- data.frame(word=names(freq), freq=freq)

# plot most frequent words along with frequency
require(ggplot2)
p <- ggplot(subset(word_freq, freq>1000), aes(word, freq))
p <- p + geom_bar(stat = "identity")
p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
p

#part 5
# plot wordcloud, a net way to express the data, size of font relates to frequency
require(wordcloud)
# plot 100 most frequent words, and add some color
wordcloud(names(freq), freq, max.words=100, rot.per = 0.2, colors=brewer.pal(6, "Dark2"))

#part 6
y_data <- as.factor(data$Label)
x_data <- as.matrix(dtm)

#part 7
sample_y <- y_data[which(data$Date <= "2014-12-31")]
out_of_sample_y <- y_data[-which(data$Date <= "2014-12-31")]
sample_x <- x_data[which(data$Date <= "2014-12-31"),]
out_of_sample_x <- x_data[-which(data$Date <= "2014-12-31"),]

#part 8
words <- c("invest","growth","grow","high","strong","lead","bankrupt","good","bull","bear","interest","market","hous","rate","oil","loss","weak","low","fear","poor","risk","stock","debt","financi","fiscal","reserv","crash","war","recess")
dtm_sentiment <- dtm[,words]

# overall frequency
freq[words]
sum(freq[words])

# comparison
senti_freq <- freq[c("invest","growth","grow","high","strong","lead","bankrupt","good","bull","bear","interest","market","hous","rate","oil","loss","weak","low","fear","poor","risk","stock","debt","financi","fiscal","reserv","crash","war","recess")]
compare <- word_freq[which(word_freq$word %in% c("invest","growth","grow","high","strong","lead","bankrupt","good","bull","bear","interest","market","hous","rate","oil","loss","weak","low","fear","poor","risk","stock","debt","financi","fiscal","reserv","crash","war","recess")),]

#part 9
# fit a standard logistic regression

# create dataframe
reg_df <- data.frame(sample_y, as.matrix(dtm_sentiment)[which(data$Date <= "2014-12-31"),])
out_words2 <- glm(reg_df$sample_y ~ ., data = reg_df, family="binomial")
summary(out_words2)

##### 10.
phat <- predict(out_words2,type="response")

# compute a ROC curve
# define a function that creates the true and false positive rates
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

glm_simple_roc <- simple_roc(reg_df$sample_y=="1", phat)
TPR <- glm_simple_roc$TPR
FPR <- glm_simple_roc$FPR

# plot the corresponding ROC curve
q <- qplot(FPR,TPR,xlab="FPR",ylab="TPR",col=I("blue"),main="ROC Curve for Logistic Regression",size=I(0.75))
# add straight 45 degree line from 0 to 1
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0)) + theme_bw()

##### 11.
# run the logistic ridge regression. Ridge is done when alpha = 0.5
out_words_ridge <- cv.glmnet(as.matrix(reg_df[-1]),as.matrix(reg_df$sample_y), family=c("binomial"), alpha = 0.5, standardize = TRUE)

# plot the coefficients, as a function of the log of the constraint value given as "lambda".
# Note, this is Lagrangian formulation, high value of constraint means it binds more
plot.glmnet(out_words_ridge$glmnet.fit, "lambda", label = TRUE)
plot.cv.glmnet(out_words_ridge)
coef(out_words_ridge)

##### 12.
# let's see how the unrestricted model fares
phat_full <- predict(out_words_ridge, newx = as.matrix(reg_df[-1]), s="lambda.1se",type="response")
glm_simple_roc_full <- simple_roc(reg_df$sample_y=="1", phat_full)

glm_simple_roc <- cbind(glm_simple_roc,Model = "GLM")
glm_simple_roc_full <- cbind(glm_simple_roc_full, Model = "GLMNET")

New_ROC <- rbind(glm_simple_roc, glm_simple_roc_full)

q <- qplot(FPR,TPR,data = New_ROC, colour = Model, xlab="FPR",ylab="TPR",main="ROC Curve for GLM and GLMNET Logistic Regression Models",size=I(0.75))
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0), col = I("black")) + theme_bw()


##### 13.
reg_df2 <- data.frame(out_of_sample_y, as.matrix(dtm_sentiment)[-which(data$Date <= "2014-12-31"),])
phat_test <- as.matrix(predict(out_words2, reg_df2[-1], type="response"))
phat_full_test <- predict(out_words_ridge, newx = as.matrix(reg_df2[-1]), type="response")

glm_test_roc <- simple_roc(reg_df2$out_of_sample_y=="1", phat_test)
glm_test_roc_full <- simple_roc(reg_df2$out_of_sample_y=="1", phat_full_test)

glm_test_roc <- cbind(glm_test_roc,Model = "GLM")
glm_test_roc_full <- cbind(glm_test_roc_full, Model = "GLMNET")

New_ROC_test <- rbind(glm_test_roc, glm_test_roc_full)

q <- qplot(FPR,TPR,data = New_ROC_test, colour = Model, xlab="FPR",ylab="TPR",main="ROC Curve for GLM and GLMNET on Test Sample",size=I(0.75))
q + geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size=I(1.0), col = I("black")) + theme_bw()

# proportion of days model predicts the right decision is when it predicts actual up and actual down. 
# This is calculated using TPR + (1-FPR) = Total right predictions
roc.curve=function(S,Y,s,print=FALSE){
  Ps=(S>s)*1
  FP=sum((Ps==1)*(Y==0))/sum(Y==0)
  TP=sum((Ps==1)*(Y==1))/sum(Y==1)
  if(print==TRUE){
    print(table(Observed=Y,Predicted=Ps))
  }
  vect=c(FP,TP)
  names(vect)=c("FPR","TPR")
  return(vect)
} 
threshold = 0.5

# Total prediction for GLM model
S = phat_test
Y = reg_df2$out_of_sample_y
roc.curve(S,Y,threshold,print=TRUE)

# Total prediction for GLMNET model
S = phat_full_test
Y = reg_df2$out_of_sample_y
roc.curve(S,Y,threshold,print=TRUE)

# Part 14
#63 day moving average
pred_ma_glm <- c(sapply(63:length(phat_test), function(i) mean(phat_test[(i-62):i])))
pred_ma_glmnet <- c(sapply(63:length(phat_full_test), function(i) mean(phat_full_test[(i-62):i])))
y_ma <- c(sapply(63:length(reg_df2$out_of_sample_y), function(i) mean(as.numeric(reg_df2$out_of_sample_y[(i-62):i]))))
dates <- as.Date(data$Date[1674:length(data$Date)])  

plot(dates,pred_ma_glm, xlab= "Date", ylab = "Predicted 63 day MA", main = "Plot for GLM model", type='l')
plot(dates,pred_ma_glmnet, xlab= "Date", ylab = "Predicted 63 day MA", main = "Plot of GLMNET model", type='l')

# Standard Regression
lm(y_ma~pred_ma_glm)$coef
summary(lm(y_ma~pred_ma_glm))$r.squared

lm(y_ma~pred_ma_glmnet)$coef
summary(lm(y_ma~pred_ma_glmnet))$r.squared




x <- 1:100
filter(x, rep(1, 3))
filter(x, rep(1, 3), sides = 1)

