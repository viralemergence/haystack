#
# Utility function for feature selection
# 

library(stringr)

# USAGE:
# select_features(FinalData, predict_column = PREDICT_COLUMN, positive_name = POSITIVE_NAME,
# 								removecols = removeCols,
# 								n_features = INPUT$nfeatures, train_proportion = INPUT$trainProportion)


# Select the top n features predicting a response variable
# - data: the final dataset, including features
# - predict_column: the column containing the response variable
# - positive_name: value in predict_column to be considered as the 'TRUE' or positive class
# - removecols: columns in data which are not to be considered features
# - n_features: number of features to retain
# - n_repeats: number of random test-train splits to perform
# - train_proportion: proportion of data to select for training
select_best_features <- function(data, predict_column, positive_name, removecols, 
																 n_features = 100, n_repeats = 100, train_proportion = 0.7) {
	
	# These are the defaults for xgboost
	tuning_params <- data.frame(eta = 0.3,
															max_depth = 6,
															subsample = 1,
															colsample_bytree = 1,
															nrounds = 150, # Has no default in xgboost, taken from Babayan et al. (2018) instead
															min_child_weight = 1,
															gamma = 0)
	
	varimps <- data.frame()
	
	for (i in 1:n_repeats) {  
		# Select training data
		trainData <- data %>%
			group_by_at(predict_column) %>%
			sample_frac(size = train_proportion, replace = FALSE) %>%
			ungroup()
		
		# Training
		caretData <- trainData %>%
			as_caret_data(labelcol = predict_column, removecols = removecols)
		
		xgboost_seed <- sample(1e8, size = 1)
		
		trainedModel <- train(x = caretData[, !colnames(caretData) == predict_column],
													y = caretData[[predict_column]],
													method = 'xgbTree',
													tuneGrid = tuning_params,
													trControl = trainControl(method = 'none', number = 1),
													seed = xgboost_seed, nthread = 1)
		
		current_importance <- varImp(trainedModel)
		varimps <- bind_rows(varimps, data.frame(Iteration = i,
																						 Feature = rownames(current_importance$importance),
																						 Importance = current_importance$importance$Overall,
																						 stringsAsFactors = FALSE))
		
	}
	
	## Summarise importance and select feature:
	#   - Taxonomy features need special handling, since most are one-hot encoded
	#   - Importance in this case is the sum of importance values across all one-hot encoded columns
	keep_features <- varimps %>%
		mutate(Feature = str_replace(.data$Feature, "(Taxonomy_[[:alpha:]]+)[.]*[[:alpha:]_]*$", "\\1")) %>% 
		group_by(.data$Iteration, .data$Feature) %>% 
		summarise(Importance = sum(.data$Importance)) %>% 
		
		# Summarise across iterations
		group_by(.data$Feature) %>%
		summarise(Importance = mean(.data$Importance)) %>%
		ungroup() %>%
		
		# Select top features
		top_n(n = n_features, wt = .data$Importance) %>%
		pull(.data$Feature)
		
	keep_features
}
