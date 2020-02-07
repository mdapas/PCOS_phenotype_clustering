################################################
# Notes 02.06.2020                             #
################################################

# This script is designed to read in a tab-delimited text file with 
#  normalized trait values for the following traits by default:
#  BMI, Testosterone, SHBG, Insulin, Glucose, DHEAS, LH, and FSH
#
# The input trait files should feature the sample_IDs as the first column, the classifications in the 2nd column, 
#  and the traits on which to train in the remaining columns.
# These variables can be modified as desired in the 'Variables' section
#
# This script evaluates four different models:
#  support vector machine, random forest, gaussian mixture model, and quadratic discriminant analysis
#
# The manuscript detailing these methods and applications is currently under review.
# Direct questions to Matthew Dapas at mdapas@uchicago.edu


################################################
# Load Libraries                               #
################################################

library(MASS)
library(e1071)
library(randomForest)
library(mclust)
library(gmm)
library(factoextra)
library(FactoMineR)


################################################
# Variables                                    #
################################################

train_input_file <- '~/train_input_file.txt' # File path to input table with training data
test_input_file <- '~/test_input_file.txt' # File path to input table
output_file <- '~/output_file.txt'

#train_input_file <- '/Users/mdapas/Documents/Research/Hayes/subtyping/genotyped_cohort.tsv'
#test_input_file <- '/Users/mdapas/Documents/Research/Hayes/subtyping/fam_cohort_76_Z.tsv'

model_cols <- c(2:10)
class_col <- 'cluster_id'
input_classes <- data.frame(clusters=c(1,2,3), subtype=c('Metabolic','Reproductive','Intermediate'), 
                            stringsAsFactors = FALSE)
methods <- c('svm','rf','gda','gmm')


################################################
# Functions                                    #
################################################

k_cross <- function(data, model, y, method='svm', K=10, seed=123) {
  suppressWarnings(RNGkind(sample.kind = 'Rounding'))  # To use random generator from R < 3.6.0
  set.seed(seed)
  
  n <- nrow(data)  # Number of observations
  data_y=data[, y]  # Response variable

  f <- ceiling(n/K)  # Number of samples in each of K partitions from nearest total divisible by K
  s <- sample(rep(1:K, f), n)  # Generate indices 1:K f times and sample n of them
    
  # K fold cross-validated error
  CV=NULL
  cat(c('\n',toupper(method)),sep='')
  cat('\nK: ')
  for (i in 1:K) {
    cat(paste(i,' '))
    test_index <- seq_len(n)[(s==i)]  # Test data indices
    train_index <- seq_len(n)[(s!=i)]  # Training data indices
    test_y <- data[test_index, y]
    
    if (tolower(method)=='svm') {
      fit <- svm(model, data[train_index,])
      pred_y <- predict(fit, data[test_index,])
    } else if (tolower(method)=='qda') {
      fit <- qda(model, data[train_index,])
      pred_y <- predict(fit, data[test_index,])$class
    } else if (tolower(method)=='gmm') {
      fit <- MclustDA(data[train_index,2:ncol(data)], data[train_index, y], verbose=F)
      pred_y <- predict(fit, data[test_index,2:ncol(data)])$classification
    } else if (tolower(method)=='rf') {
      fit <- randomForest(model, data[train_index,])
      pred_y <- predict(fit, data[test_index,])
    } else {
      print('invalid method')
    }

    #error rate
    error <- mean(test_y!=pred_y)
    CV <- c(CV,error)

  }
  cat(c('\nError: ',mean(CV),'\n'))
  #Output
  list(call=model, K=K, seed=seed, error_rate=mean(CV))
}


classify <- function(train_data, model, y, test_data, method, graph=TRUE, graph_cols=c(1,2), seed=123) {
  suppressWarnings(RNGkind(sample.kind = 'Rounding'))  # To use random generator from R < 3.6.0
  set.seed(seed)  
  if (tolower(method)=='svm') {
    fit <- svm(model, train_data, probability=TRUE)
    pred_y <- predict(fit, test_data, decision.values=TRUE, probability=TRUE)
    post_y <- attr(pred_y,'probabilities')
  } else if (tolower(method)=='qda') {
    fit <- qda(model, train_data)
    pred_y <- predict(fit, test_data)$class
    post_y <- pred_y$post
  } else if (tolower(method)=='gmm') {
    fit <- MclustDA(train_data[,2:ncol(train_data)], train_data[, 1], verbose=F)
    pred_y <- predict(fit, test_data[,2:ncol(test_data)])$classification
    post_y <- data.frame(pred_y$z)
  } else if (tolower(method)=='rf') {
    fit <- randomForest(model, train_data)
    pred_y <- predict(fit, test_data)
    post_y <- predict(fit, test_data, type='prob')
  } else {
    print('invalid method')
  }
  
  if (graph) {
    plot(test_data[,graph_cols], pch=as.numeric(pred_y), col=as.numeric(pred_y))
    legend('topleft',legend=c('Pred 1','Pred 2', 'Pred 3'),pch=1:3,col=1:3, cex=.7)
  }
  
  list(call=model, method=method, seed=seed, totals=table(pred_y), classes=pred_y)
}


################################################
# Main                                         #
################################################

training_df <- read.table(train_input_file, header=T)
predict_df <- read.table(test_input_file, header=T)

for (i in seq(1:length(methods)))  {
  min_err <- ifelse(i==1,1,min_err)
  model_cross <- k_cross(training_df[,model_cols], as.factor(cluster_id) ~ ., y=class_col, method=methods[i])
  if (model_cross$error_rate < min_err) {
    min_err <- model_cross$error_rate
    best_method <- methods[i]
  }
}

predicts <- classify(training_df[,model_cols], as.factor(cluster_id) ~ ., 
                     y=class_col, predict_df[,model_cols], method=best_method, graph_cols=c('SHBG','BMI'))

predict_df$clusters <- predicts$classes

# write.table(predict_df, output_file, sep='\t', quote=F, row.names=F)

predict_df$subtype <- sapply(predict_df$clusters, function(x) input_classes[input_classes$clusters==x,'subtype'])

pca <- PCA(predict_df[model_cols[-1]], graph = T)
#pca$ind$coord[,1] <- pca$ind$coord[,1]*-1
#pca$var$coord[,1] <- pca$var$coord[,1]*-1
fviz_pca_biplot(pca, col.ind=predict_df$subtype, palette = c('grey',rgb(0.81,0.3,0.25,.85), rgb(.33,.5,.77,1)),
                     addEllipses=T, label='var', col.var='black', repel=T, show.clust.cent=F,
                     legend.title='Subtype', pointsize=5, invisible='quali', #pointshape=19,
                     xlab=paste('PC1 (', round(pca$eig[1,2],2), '%)', sep=''),
                     ylab=paste('PC2 (', round(pca$eig[2,2],2), '%)', sep=''),
                     title='') + theme(legend.text=element_text(size=10))
