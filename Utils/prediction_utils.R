## Utility functions: prediction of binary labels / zoonotic potential categories

# Calculate Youden's J
j_stat <- function(cutoff, obs, prob) {
  stopifnot(length(obs) == length(prob))
  
  pred <- prob > cutoff
  sensitivity = sum(obs & (obs == pred)) / sum(obs)
  specificity = sum(!obs & (obs == pred)) / sum(!obs)
  
  sensitivity + specificity - 1
}

# Calculate balance: the absolute distance between sensitivity and specificity
# (negative, so this can be maximized)
balance_stat <- function(cutoff, obs, prob) {
  stopifnot(length(obs) == length(prob))
  
  pred <- prob > cutoff
  sensitivity = sum(obs & (obs == pred)) / sum(obs)
  specificity = sum(!obs & (obs == pred)) / sum(!obs)
  
  -abs(sensitivity - specificity)
}


# Find the optimal cutoff balancing sensitivity and specificity
# - Achieved by minimizing a give test statistic: either Youdens' J (informedness) or 
#   the best balance between sensitivity and specificity
find_best_cutoff <- function(observed_labels, predicted_score, 
                             positive_value = "True", 
                             increment_size = 0.0001,
                             stat = c("informedness", "balance")) {
  
  stat <- match.arg(stat)
  stat_fun <- switch(stat,
                     informedness = j_stat,
                     balance = balance_stat)
  
  obs <- observed_labels == positive_value
  try_cutoffs <- seq(0, 1, by = increment_size)
  
  distances <- vapply(try_cutoffs, stat_fun, 
                      FUN.VALUE = numeric(1),
                      obs = obs, prob = predicted_score)
  
  try_cutoffs[which.max(distances)]
}



# Convert scores into zoonotic potential / priority categories
prioritize <- function(lower_bound, upper_bound, median, cutoff) {
  stopifnot(length(cutoff) == 1)
  stopifnot(cutoff > 0 & cutoff < 1)
  
  p <- if_else(lower_bound > cutoff, 'Very high',
               if_else(median > cutoff, 'High',
                       if_else(upper_bound > cutoff, 'Medium',
                               'Low')))
  factor(p, levels = c('Low', 'Medium', 'High', 'Very high'))
}
