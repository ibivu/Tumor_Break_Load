auprcSummary <- function(data, lev = NULL, model = NULL){
  # performance estimation function for the implementation of the AUC of 
  # the precision recall curve
  
  library(PRROC)
  
  # Select lable of interest  
  index_class2 <- data$obs == "abr"
  index_class1 <- data$obs == "wt"
  
  # calculate the AUC of the precision recall curve
  the_curve <- pr.curve(data$abr[index_class2], data$abr[index_class1], curve = FALSE)
  out <- the_curve$auc.integral
  names(out) <- "AUPRC"
  
  out
  
}

# function to run and select best performing features
modelTraining <- function(train_set, train_labels, method, folds, seq_length, tune_length, features, min_features) {
  # train_set: NxM expression dataframe for parameter tuning and training
  # train_labels: vector of classification labels for train data
  # method: Model method to use
  # folds: amount of folds for cross validation (int)
  # seq_length: number of steps to take for feature selection sequence
  # tune_length: number of steps to take for tuning sequence
  # features: max features to consider for feature selection
  # min_features: minimum features to consider for feature selection
  
  print("start feature selection")
  
  
  rfFuncs$summary <- auprcSummary
  
  # setup training and feature selection strategystrategy
  filter.control = rfeControl(functions = rfFuncs, method = "repeatedCV", 
                              number = folds, verbose = F)
  train.control = trainControl(method = "cv", number = folds, 
                               summaryFunction = auprcSummary, classProbs = T)
  
  # feature selection using recursive feature elimination
  feature.results = rfe(train_set, train_labels, rfeControl = filter.control, 
                        sizes = round(seq(min_features, features, length = seq_length)), 
                        metric = "AUPRC", maximize = T, method = method)
  
  top_features = feature.results$optVariables
  
  print("end feature selection")
  
  train_set <- train_set[top_features]
  
  print("start tuning")
  
  # parameter tuning with the selected feature set
  train.results = train(train_set, train_labels, method = method, trControl = train.control,
                        tuneLength = tune_length, metric = "AUPRC")
  
  print("end tuning")
  
  return(list(features = top_features, train_result = train.results, feature_results = feature.results))
}

# function to subset data from partition sets
subset_data <- function(data, mapping, labels) {
  sample_names <- mapping[match(names(labels), mapping$File.ID), "File.ID.RNA"]
  subset_data <- data[match(sample_names, data$X),]
  rownames(subset_data) <- names(labels)
  subset_data <- subset(subset_data, select = -c(X))
  return(subset_data)
}

run_model <- function(training_data, training_labels, 
                      method, folds, tune_length, features, 
                      seq_length, iteration_length, cores=4, min_features = 0.1*features, parallel = T) {
  # training_data: NxM expression dataframe for parameter tuning and training
  # training_labels: vector of classification labels for train data
  # method: Model method to use
  # folds: amount of folds for cross validation (int)
  # seq_length: number of steps to take for feature selection sequence
  # tune_length: number of steps to take for tuning sequence
  # features: max features to consider for feature selection
  # min_features: minimum features to consider for feature selection
  # iteration_length: number of iterations to rerun the model training
  # cores: number of cores to use for parallel processing
  # parallel: boolean value to initiate parallel processing
  
  # setup storage for model metrics
  best_model <- NULL
  best_features <- NULL
  performance_list <- c(0)
  max_performance <- 0
  
  # if parallel true initiate cores
  if (parallel) {
    core.cluster = makeCluster(cores)
    registerDoParallel(core.cluster)
  }
  
  print("start")
  for (i in 1:iteration_length){
    set.seed(i)
    
    # train model
    model_result <- modelTraining(training_data,  
                                  training_labels, method = method, 
                                  seq_length =  seq_length, folds = folds, 
                                  tune_length = tune_length, 
                                  features = features, min_features = min_features)
    
    performance_list <- c(performance_list, 
                          max(model_result$train_result$results$AUPRC))
    
    # select best performing trained model
    if (max(model_result$train_result$results$AUPRC) > max_performance) {
      best_model <- model_result$train_result
      best_features <- model_result$feature_results
      max_performance <- max(model_result$train_result$results$AUPRC)
    }
    print(i)
  }
  print("end")
  
  # close cores
  if (parallel) {
    stopCluster(core.cluster)
    registerDoSEQ()
  }
  
  return(list(best_model = best_model, performance_list = performance_list[-1], 
              max_performance = max_performance, best_features = best_features))
}

random_permutation_predictor <- function(labels, data, permutations = 1000, 
                                         folds = 5, tune_length = 3, features = 300,
                                         seq_length = 5, iteration_length = 1,
                                         cores = 45, min_features = 10) {
  # Assess the random model strategy performance
  # labels: vector of classification labels
  # data: NxM expression dataframe
  # permutations: number of permutations to try, influences the atainable pvalue
  # folds: amount of folds for cross validation (int)
  # seq_length: number of steps to take for feature selection sequence
  # tune_length: number of steps to take for tuning sequence
  # features: max features to consider for feature selection
  # min_features: minimum features to consider for feature selection
  # iteration_length: number of iterations to rerun the model training
  # cores: number of cores to use for parallel processing
  
  permutation_results <- c()
  
  for (i in 1:permutations) {
    print(i)
  
    names(labels) <- sample(names(labels))
    
    # create partitions
    splits <- createDataPartition(labels, times = 1, p = 0.65, list = TRUE)
    train_labels <- labels[splits$Resample1]
    test_labels <- labels[-splits$Resample1]
    
    # subsample rna-seq data
    RNA_train <- data[names(train_labels),]
    RNA_test <- data[names(test_labels),]
    
    # train random model
    results <- run_model(RNA_train, train_labels, method = "rf", folds = folds,
                         tune_length = tune_length, features = features, 
                         seq_length = seq_length, iteration_length = iteration_length,
                         cores = cores, min_features = min_features, parallel = F)
    
    prediction <- predict(results$best_model, RNA_test, type = "prob")
    
    roc_bin <- rocit(score = prediction$abr, class = test_labels, negref = "wt")
    
    permutation_results <- c(permutation_results, roc_bin$AUC)
  }
  return(permutation_results)
}