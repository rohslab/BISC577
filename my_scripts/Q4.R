library(DNAshapeR)
library(caret)

### Part 1: Mad dataset

## Predict DNA shapes
fn_fasta <- "/Users/apple/BISC577HW2/bsun/hw2/gcPBM/Mad.txt.fa"
pred <- getShape(fn_fasta)

## Encode feature vectors
featureType.1 <- c("1-mer") # "1-mer" model
featureType.2 <- c("1-mer", "1-shape") # "1-mer+1-shape" model
featureVector.1 <- encodeSeqShape(fn_fasta, pred, featureType.1)
featureVector.2 <- encodeSeqShape(fn_fasta, pred, featureType.2)
#head(featureVector)

## Build MLR model by using Caret
# Data preparation
fn_exp <- "/Users/apple/BISC577HW2/bsun/hw2/gcPBM/Mad.txt"
exp_data <- read.table(fn_exp)
df.1 <- data.frame(affinity=exp_data$V2, featureVector.1)
df.2 <- data.frame(affinity=exp_data$V2, featureVector.2)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

## The Question request L2-regularized MLR
# Prediction without L2-regularized
#model.1 <- train (affinity~ ., data = df.1, trControl=trainControl, 
#                method = "lm", preProcess=NULL)
#model.2 <- train (affinity~ ., data = df.2, trControl=trainControl, 
                  method = "lm", preProcess=NULL)
#summary(model.2)

# Prediction with L2-regularized
model.l2.1 <- train(affinity~., data = df.1, trControl=trainControl, 
                method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.l2.2 <- train(affinity~., data = df.2, trControl=trainControl, 
                    method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))

Mad.Rsquare.1 <- mean(model.l2.1$results$Rsquared[!is.na(model.l2.1$results$Rsquared)])
cat("The R^2 of 1-mer model for Mad is: ", Mad.Rsquare.1)
Mad.Rsquare.2 <- mean(model.l2.2$results$Rsquared[!is.na(model.l2.2$results$Rsquared)])
cat("The R^2 of 1-mer model for Mad is: ", Mad.Rsquare.2)

### Part 3: Myc dataset

## Predict DNA shapes
fn_fasta <- "/Users/apple/BISC577HW2/bsun/hw2/gcPBM/Myc.txt.fa"
pred <- getShape(fn_fasta)

## Encode feature vectors
featureType.1 <- c("1-mer") # "1-mer" model
featureType.2 <- c("1-mer", "1-shape") # "1-mer+1-shape" model
featureVector.1 <- encodeSeqShape(fn_fasta, pred, featureType.1)
featureVector.2 <- encodeSeqShape(fn_fasta, pred, featureType.2)
#head(featureVector)

## Build MLR model by using Caret
# Data preparation
fn_exp <- "/Users/apple/BISC577HW2/bsun/hw2/gcPBM/Myc.txt"
exp_data <- read.table(fn_exp)
df.1 <- data.frame(affinity=exp_data$V2, featureVector.1)
df.2 <- data.frame(affinity=exp_data$V2, featureVector.2)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

## The Question request L2-regularized MLR
# Prediction without L2-regularized
#model.1 <- train (affinity~ ., data = df.1, trControl=trainControl, 
#                method = "lm", preProcess=NULL)
#model.2 <- train (affinity~ ., data = df.2, trControl=trainControl, 
method = "lm", preProcess=NULL)
#summary(model.2)

# Prediction with L2-regularized
model.l2.1 <- train(affinity~., data = df.1, trControl=trainControl, 
                    method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.l2.2 <- train(affinity~., data = df.2, trControl=trainControl, 
                    method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))

Mad.Rsquare.1 <- mean(model.l2.1$results$Rsquared[!is.na(model.l2.1$results$Rsquared)])
cat("The R^2 of 1-mer model for Myc is: ", Mad.Rsquare.1)
Mad.Rsquare.2 <- mean(model.l2.2$results$Rsquared[!is.na(model.l2.2$results$Rsquared)])
cat("The R^2 of 1-mer model for Myc is: ", Mad.Rsquare.2)


### Part 2: Max dataset

## Predict DNA shapes
fn_fasta <- "/Users/apple/BISC577HW2/bsun/hw2/gcPBM/Max.txt.fa"
pred <- getShape(fn_fasta)

## Encode feature vectors
featureType.1 <- c("1-mer") # "1-mer" model
featureType.2 <- c("1-mer", "1-shape") # "1-mer+1-shape" model
featureVector.1 <- encodeSeqShape(fn_fasta, pred, featureType.1)
featureVector.2 <- encodeSeqShape(fn_fasta, pred, featureType.2)
#head(featureVector)

## Build MLR model by using Caret
# Data preparation
fn_exp <- "/Users/apple/BISC577HW2/bsun/hw2/gcPBM/Max.txt"
exp_data <- read.table(fn_exp)
df.1 <- data.frame(affinity=exp_data$V2, featureVector.1)
df.2 <- data.frame(affinity=exp_data$V2, featureVector.2)

# Arguments setting for Caret
trainControl <- trainControl(method = "cv", number = 10, savePredictions = TRUE)

## The Question request L2-regularized MLR
# Prediction without L2-regularized
#model.1 <- train (affinity~ ., data = df.1, trControl=trainControl, 
#                method = "lm", preProcess=NULL)
#model.2 <- train (affinity~ ., data = df.2, trControl=trainControl, 
method = "lm", preProcess=NULL)
#summary(model.2)

# Prediction with L2-regularized
model.l2.1 <- train(affinity~., data = df.1, trControl=trainControl, 
                    method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))
model.l2.2 <- train(affinity~., data = df.2, trControl=trainControl, 
                    method = "glmnet", tuneGrid = data.frame(alpha = 0, lambda = c(2^c(-15:15))))

Mad.Rsquare.1 <- mean(model.l2.1$results$Rsquared[!is.na(model.l2.1$results$Rsquared)])
cat("The R^2 of 1-mer model for Max is: ", Mad.Rsquare.1)
Mad.Rsquare.2 <- mean(model.l2.2$results$Rsquared[!is.na(model.l2.2$results$Rsquared)])
cat("The R^2 of 1-mer model for Max is: ", Mad.Rsquare.2)
