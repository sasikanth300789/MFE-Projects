# clean-up R script for Problem Set 4

require(tm)
require(glmnet)
require(ggplot2)
require(pROC)
require(SnowballC)

# Read in the data, use the right file path for you in the below
data <- read.csv('C:/Users/cudoh/UCLA/2016 Fall Quarter/Data Analytics/HW/HW4/DJIA_Headline_News.csv', stringsAsFactors = FALSE)

################### Clean up script #########################

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


################ End clean up script #############################





setwd('C:/Users/cudoh/UCLA/2016 Fall Quarter/Data Analytics/HW/HW4')

##### 1.
main_corpus <- Corpus(VectorSource(data$all))


##### 2.
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

# stemming words, i.e., keep only the stem so as not to differentially count investing, invest, invests
# taking out common word endings such as 'ing', 'es', and 's'
main_corpus <- tm_map(main_corpus, PlainTextDocument)
main_corpus <- tm_map(main_corpus, stemDocument)

# finally, let's get rid of all the extra white space in the document so all words are only 
# separated by one space
#main_corpus <- tm_map(main_corpus, stripWhitespace)

# we are done cleaning the document, just now make sure it is treated as plain text
#main_corpus <- tm_map(main_corpus, PlainTextDocument)


##### 3.
dtm <- DocumentTermMatrix(main_corpus)
inspect(dtm[5:10,801:810])


##### 4.
# organize words by frequency
freq <- colSums(as.matrix(dtm))
ord_corpus <- order(freq)

#  identify words that appear frequently
# create a convenient data.frame
word_freq <- data.frame(word=names(freq), freq=freq)
# taking out the separation character "<s>" for the chart
word_freq2 <- word_freq[-1,]

# plot most frequent words along with frequency
p <- ggplot(subset(word_freq2, freq>1000), aes(word, freq))
p <- p + geom_bar(stat = "identity")
p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
p


##### 5.
# plot wordcloud, a net way to express the data, size of font relates to frequency
require(wordcloud)
# plot 100 most frequent words, and add some color
wordcloud(names(freq[-1]), freq[-1], max.words=100, rot.per = 0.2, colors=brewer.pal(6, "Dark2"))


##### 6.
y_data <- as.factor(data$Label)
y_data <- data$Label
x_data <- as.matrix(dtm)


##### 7.
sample_date <- which(data$Date == "2014-12-31")
sample_x <- x_data[1:sample_date,]
out_of_sample_x <- x_data[(sample_date+1):length(data$Date),]
sample_y <- y_data[1:sample_date]
out_of_sample_y <- y_data[(sample_date+1):length(data$Date)]
#sample_length <- length(sample_x$Date)



##### 8.
words <- c("invest","growth","grow","high","strong","lead","bankrupt","good","bull","bear","interest","market","hous","rate","oil","loss","weak","low","fear","poor","risk","stock","debt","financi","fiscal","reserv","crash","war","recess")
dtm_sentiment <- dtm[,words]

freq[words]
sum(freq[words])

##### 9.
reg_df <- data.frame(sample_y, as.matrix(dtm_sentiment)[1:sample_date,])
#out_words <- glm(sample$Label ~ x_data[1:sample_length,], family="binomial")
out_words2 <- glm(reg_df$sample_y ~ ., data = reg_df, family="binomial")

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
out_words_ridge <- cv.glmnet(as.matrix(reg_df[-1]),as.vector(reg_df$sample_y), family=c("binomial"),alpha = 0.5, standardize = TRUE)

# plot the coefficients, as a function of the log of the constraint value given as "lambda".
# Note, this is Lagrangian formulation, high value of constraint means it binds more
plot.glmnet(out_words_ridge$glmnet.fit, "lambda", label = TRUE)
plot.cv.glmnet(out_words_ridge)


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
reg_df2 <- data.frame(out_of_sample_y, as.matrix(dtm_sentiment)[(sample_date+1):(sample_date+length(out_of_sample_y)),])
a <- as.matrix(predict(out_words2, (reg_df2[-1]), type="response"))
a2 <- predict(out_words_ridge, newx = as.matrix(reg_df2[-1]), s="lambda.1se",type="response")
