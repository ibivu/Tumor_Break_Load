auprcSummary <- function(data, lev = NULL, model = NULL){
  
  library(PRROC)
  
  index_class2 <- data$obs == "abr"
  index_class1 <- data$obs == "wt"
  
  the_curve <- pr.curve(data$abr[index_class2], data$abr[index_class1], curve = FALSE)
  out <- the_curve$auc.integral
  names(out) <- "AUPRC"
  
  out
  
}

# function to run and select best performing features
modelTraining <- function(train_set, train_labels, method, folds, seq_length, tune_length, features, min_features) {
  print("start feature selection")
  # feature_set: nxm expression dataframe for feature selection
  # train_set: nxm expression dataframe for parameter tuning and training
  # feature_labels: vector of classification labels for feature data
  # train_labels: vector of classificiation lables for train data
  # method: Model method to use
  # folds: ammount of folds for cross validation (int)
  # seq_length: number of steps to take for feature selection sequence
  # tune_length: number of steps to take for tuning sequence
  # features: max features to consider for feature selection
  # min_features: minimum features to consider for feature selection
  
  rfFuncs$summary <- auprcSummary
  
  filter.control = rfeControl(functions = rfFuncs, method = "repeatedCV", number = folds, verbose = F)
  train.control = trainControl(method = "cv", number = folds, summaryFunction = auprcSummary, classProbs = T)
  
  feature.results = rfe(train_set, train_labels, rfeControl = filter.control, 
                        sizes = round(seq(min_features, features, length = seq_length)), 
                        metric = "AUPRC", maximize = T, method = method)
  
  top_features = feature.results$optVariables
  
  print("end feature selection")
  
  train_set <- train_set[top_features]
  
  print("start tuning")
  
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
                      seq_length, iteration_length, cores, min_features = 0.1*features, parallel = T) {
  
  best_model <- NULL
  best_features <- NULL
  performance_list <- c(0)
  max_performance <- 0
  
  if (parallel) {
    core.cluster = makeCluster(cores)
    registerDoParallel(core.cluster)
  }
  
  print("start")
  for (i in 1:iteration_length){
    set.seed(i)
    
    # # run model training only if there are enough bps
    # splits <- createDataPartition(training_labels, times = 1, p = 0.5, list = TRUE)
    # feature_labels <- training_labels[splits$Resample1]
    # train_labels <- training_labels[-splits$Resample1]
    # 
    # # subsample rna-seq data
    # feature_data <- training_data[names(feature_labels),]
    # train_data <- training_data[names(train_labels),]
    
    model_result <- modelTraining(training_data,  
                                  training_labels, method = method, 
                                  seq_length =  seq_length, folds = folds, 
                                  tune_length = tune_length, 
                                  features = features, min_features = min_features)
    
    performance_list <- c(performance_list, 
                          max(model_result$train_result$results$AUPRC))
    
    if (max(model_result$train_result$results$AUPRC) > max_performance) {
      best_model <- model_result$train_result
      best_features <- model_result$feature_results
      max_performance <- max(model_result$train_result$results$AUPRC)
    }
    print(i)
  }
  print("end")
  
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
  permutation_results <- c()
  
  for (i in 1:permutations) {
    print(i)
  
    names(labels) <- sample(names(labels))
    
    # run model training only if there are enough bps
    splits <- createDataPartition(labels, times = 1, p = 0.65, list = TRUE)
    train_labels <- labels[splits$Resample1]
    test_labels <- labels[-splits$Resample1]
    
    # subsample rna-seq data
    RNA_train <- data[names(train_labels),]
    RNA_test <- data[names(test_labels),]
    
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