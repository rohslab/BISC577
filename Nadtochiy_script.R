######################################
# Install packages
######################################
# Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite()
# DNAshapeR
biocLite("DNAshapeR")
# Biostrings
biocLite("Biostrings")
# Caret
install.packages("caret")
install.packages("e1071")
install.packages("ROCR")
install.packages("ggplot2")
install.packages("grid")

library(ggplot2)
library(grid)
library(ROCR)
library(Biostrings)
library(DNAshapeR)
library(caret)

######################################
#Q4. Multiple Linear Regression (MLR)
#    for "1-mer+shape"
######################################

## set directory
setwd("C:/Users/pc/Dropbox/Education/USC/Remo/BISC577-master")

## Predict DNA shapes
fn_fasta <- "gcPBM/Max.fa"
pred <- getShape(fn_fasta)

## Encode feature vectors
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn_fasta, pred, featureType)
head(featureVector)
############

## Build MLR model by using Caret
# Data preparation
fn_exp <- "gcPBM/Max.txt"
exp_data <- read.table(fn_exp)
df <- data.frame(affinity=exp_data$V2, featureVector)
################

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction without L2-regularized
#model <- train (affinity~ ., data = df, trControl=trainControl, 
#                method = "lm", preProcess=NULL)
#summary(model)

# Prediction with L2-regularized
model2 <- train(affinity~., data = df, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model2

#Get values of R squared
#Rs_mad<-na.omit(model2$results$Rsquared)
#Rs_mad

## SEQUENCE Encode feature vectors
# (this function is unavailable)
#featureVector <- encodeKMerSeq(1, fn_fasta)
#head(featureVector)

## cut out SEQUENCE-ONLY from feature vector
featureSeq <- featureVector[,c(1:144)]
head(featureSeq)

## Build MLR model by using Caret
# Data preparation
df <- data.frame(affinity=exp_data$V2, featureSeq)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

# Prediction without L2-regularized
#model <- train (affinity~ ., data = df, trControl=trainControl, 
#                method = "lm", preProcess=NULL)
#summary(model)

# Prediction with L2-regularized
model2 <- train(affinity~., data = df, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model2

#values of R squared
#Rs_mad<-na.omit(model2$results$Rsquared)
#Rs_mad

#####################################
#Q5. Comparison of two models
#####################################

## Theme
my.theme <- theme(
  plot.margin = unit(c(0.2, 0.5, 0.2, 0.2), 'cm'),
  axis.text.x = element_text(colour='black', size=12),
  axis.text.y = element_text(colour='black', size=12),
  axis.title.x = element_text(colour='black', size=12),
  axis.title.y = element_text(colour='black', size=12),
  legend.text = element_text(colour='black', size=12),
  legend.position = c(.9, .6),
  legend.key = element_rect(fill='white'),
  legend.title=element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.text = element_text(colour ="black"),
  axis.ticks = element_line(colour = "black")
)

## Data preparation
data_seq <- c(0.7753123, 0.7852302, 0.7789938)
data_seq_shape <- c(0.8628297, 0.8640685, 0.8543180)
my_dataset<- c('Mad','Max','Myc')

## Ploting
ggplot() +
  geom_point(aes(x = data_seq, y = data_seq_shape, color = factor(my_dataset)), size=1) +
  geom_abline(slope=1) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +
  xlab(expression("1 mer R "^2))+ylab(expression("1 mer+sequence R "^2))+
  coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  my.theme 

######################################
# Emsemble plots
######################################

# Initialization
#library(DNAshapeR)

# Extract sample sequences
fnb <- "CTCF/bound_500.fa"
fnu <- "CTCF/unbound_500.fa"

# Predict DNA shapes for bound and unbound
pred_b <- getShape(fnb)
pred_u <- getShape(fnu)


# Generate ensemble plots
plotShape(pred_b$MGW, main="MGW", xlim = c(0,500), ylim = c(4.95,5.13))
par(new=T)
plotShape(pred_u$MGW, xlim = c(0,500), ylim = c(4.95,5.13), axes = F,
          colDots = rgb( 1, 0, 0, 0.1),colLine = "red")
par(new=F)


plotShape(pred_b$Roll, main="Roll", xlim = c(0,500), ylim = c(-1.15,-0.7))
par(new=T)
plotShape(pred_u$Roll, xlim = c(0,500), ylim = c(-1.15,-0.7), axes = F,
          colDots = rgb( 1, 0, 0, 0.1), colLine = "red")
par(new=F)

plotShape(pred_b$HelT,main="HelT", xlim = c(0,500), ylim = c(34.1,34.6))
par(new=T)
plotShape(pred_u$HelT, xlim = c(0,500), ylim = c(34.1,34.6), axes = F,
          colDots = rgb( 1, 0, 0, 0.1),colLine = "red")
par(new=F)

plotShape(pred_b$ProT, main= "ProT", xlim = c(0,500), ylim = c(-7.5,-5.6))
par(new=T)
plotShape(pred_u$ProT, xlim = c(0,500), ylim = c(-7.5,-5.6), axes = F,
          colDots = rgb( 1, 0, 0, 0.1),colLine = "red")
par(new=F)

######################################
#Q8. Logistic regression on ChIP-seq data
# "1-mer + shape"
######################################

workingPath <- "CTCF/"

## Generate data for the classifcation (assign Y to bound and N to non-bound)
# bound
boundFasta <- readDNAStringSet(paste0(workingPath, "bound_30.fa"))
sequences <- paste(boundFasta)
boundTxt <- data.frame(seq=sequences, isBound="Y")

# non-bound
nonboundFasta <- readDNAStringSet(paste0(workingPath, "unbound_30.fa"))
sequences <- paste(nonboundFasta)
nonboundTxt <- data.frame(seq=sequences, isBound="N")

# merge two datasets
writeXStringSet( c(boundFasta, nonboundFasta), paste0(workingPath, "ctcf.fa"))
exp_data <- rbind(boundTxt, nonboundTxt)


## DNAshapeR prediction
pred <- getShape(paste0(workingPath, "ctcf.fa"))


## Encode feature vectors
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(paste0(workingPath, "ctcf.fa"), pred, featureType)
df <- data.frame(isBound = exp_data$isBound, featureVector)


## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
summary(model)

## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance, main = {" 1-mer + shape \nAUC = 0.842"})

## Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc

######################################
# '1-mer', sequence only
######################################

## cut-out sequence-only part (30 nucleotides => first 120 enteries)
featureSeq <- featureVector[,c(1:120)]
df <- data.frame(isBound = exp_data$isBound, featureSeq)


## Logistic regression
# Set parameters for Caret
trainControl <- trainControl(method = "cv", number = 10, 
                             savePredictions = TRUE, classProbs = TRUE)
# Perform prediction
model <- train(isBound~ ., data = df, trControl = trainControl,
               method = "glm", family = binomial, metric ="ROC")
summary(model)

## Plot AUROC
prediction <- prediction( model$pred$Y, model$pred$obs )
performance <- performance( prediction, "tpr", "fpr" )
plot(performance, main = {" 1-mer \nAUC = 0.840"})

## Caluculate AUROC
auc <- performance(prediction, "auc")
auc <- unlist(slot(auc, "y.values"))
auc










