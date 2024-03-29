#! Rscript
# Part 3 of Zoonosis prediction pipeline
#	 - Combine results from models fitted to different train/test splits to
#	   make model-averaged predictions

# Input arguments:
library(argparse)

parser <- ArgumentParser(description = 'Calculate bagged predictions from test splits')
parser$add_argument('RandomSeed', type = 'integer', 
										help = 'a random seed (used by both R and xgboost)')
parser$add_argument('RunID', type = 'character', 
										help = 'an identifier for this run, used as a base to name output files')

parser$add_argument('--Ntop', type = 'integer', default = 15,
										help = 'Number of top models to use (default: 15)')

parser$add_argument('--nthreads', type = 'integer', default = 8,
										help = 'Number of threads/cores to use (default: 8)')

parser$add_argument('--skip_checks', default = FALSE, action='store_true',
										help = 'Do not check if specified number of top models can be achieved for all viruses (default: false). A warning will still be emited.')

INPUT <- parser$parse_args()


# Constants	based on input args:
set.seed(INPUT$RandomSeed)

RUN_ID <- INPUT$RunID  # Identify a run from TrainAndValidate.R to analyse
N_TOP <- INPUT$Ntop  # Predictions for each virus will be based on the top N_TOP models containing predictions for that species


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Libraries / utils --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(dplyr)
library(tidyr)
library(ModelMetrics)
library(parallel)

# betacal can't be installed via conda:
if (!require("betacal")) {
  install.packages("betacal", repos = "https://cloud.r-project.org", verbose = FALSE)
  library(betacal)
}


RunDataDir <- file.path('RunData', RUN_ID)

source(file.path('Utils', 'rundata_utils.R'))
source(file.path('Utils', 'calibration_utils.R'))
source(file.path('Utils', 'prediction_utils.R'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Data ---------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
predictions_all <- file.path(RunDataDir, paste0(RUN_ID, '_Predictions.rds')) %>% 
	readRDS()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Calculate AUC for each iteration -----------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# - For each virus, calculate AUCs for all models in which it occurs in the test set,
#   but *without* using that virus in the calculation
get_auc <- function(virus_name, predictions) {
	keep_models <- predictions$Iteration[predictions$LatestSppName == virus_name]
	
	predictions %>% 
		filter(.data$Iteration %in% keep_models) %>% 
		filter(.data$LatestSppName != virus_name) %>% 
		group_by(.data$Iteration) %>% 
		summarise(AUC = auc(actual = .data$InfectsHumans, predicted = .data$RawScore)) %>% 
		mutate(LatestSppName = virus_name)
}

preds <- predictions_all %>% 
	filter(.data$Dataset == 'test') %>% 
	mutate(InfectsHumans = .data$InfectsHumans == 'True')

AUCs <- mclapply(unique(preds$LatestSppName), get_auc, predictions = preds, mc.cores = INPUT$nthreads) %>% 
	bind_rows()


predictions_all <- predictions_all %>% 
	left_join(AUCs, by = c('Iteration', 'LatestSppName'))



# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Calibrate probabilities --------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
calibration_preds <- predictions_all %>% 
	filter(.data$Dataset == 'calibration')

predictions_all <- predictions_all %>% 
	group_by(.data$Iteration) %>% 
	group_modify(~ calibrate_preds(.x, calibration_preds = calibration_preds))


# Use (calibrated) calibration set to find a cutoff too:
calibration_preds <- predictions_all %>% 
  filter(.data$Dataset == 'calibration')

cutoff <- find_best_cutoff(calibration_preds$InfectsHumans, calibration_preds$CalibratedScore)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Calculate averaged predictions -------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# The bagged prediction for each virus is based on the N_TOP best models which had that virus in the test set:

# Check that all viruses occur in at least N_TOP models:
ModelCounts <- predictions_all %>% 
	filter(.data$Dataset == 'test') %>% 
	group_by(.data$LatestSppName) %>% 
	summarise(N = n()) %>% 
	.$N

if (min(ModelCounts) < N_TOP) {
	err_message <- paste0("Not all viruses occur in the test sets of enough models (min = ", min(ModelCounts), 
												"; requested = ", N_TOP, " top models)")
	
	if (! INPUT$skip_checks)
		stop(err_message)
	
	# If asked to continue despite this issue, raise a warning instead:
	warning(err_message)
}

message("Viruses occur in the test sets of ", min(ModelCounts), " to ", max(ModelCounts),
				" models. Keeping predictions from the best ", N_TOP, " models")


# Make predictions
bagged_predictions <- predictions_all %>%
	filter(.data$Dataset == 'test') %>% 
	group_by(.data$LatestSppName) %>% 
	mutate(Rank = rank(-.data$AUC, ties.method = 'random')) %>%     # Reverse rank, i.e. the best model has rank 1
	filter(.data$Rank <= N_TOP) 

if (any(bagged_predictions$CalibrationMethod == 'ab'))
	warning('Three-paramater calibration failed in some cases. Using a mixture of two- and three-parameter calibrated probabilities for bagging')

bagged_predictions <- bagged_predictions %>%
  group_by(.data$LatestSppName) %>% 
	summarise(N = n(),
						InfectsHumans = unique(.data$InfectsHumans),
						BagScore_Lower = quantile(.data$CalibratedScore, probs = 0.05/2),
						BagScore_Upper = quantile(.data$CalibratedScore, probs = 1 - 0.05/2),
						BagScore = mean(.data$CalibratedScore))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Add binary prediction / zoonotic potential level -------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
bagged_predictions <- bagged_predictions %>% 
  group_by(.data$LatestSppName) %>% 
  mutate(bagged_prediction = .data$BagScore >= cutoff,
         zoonotic_potential = prioritize(lower_bound = .data$BagScore_Lower,
                                         upper_bound = .data$BagScore_Upper,
                                         median = .data$BagScore,
                                         cutoff = cutoff)) %>%
  ungroup()


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output results -----------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
out_dir <- file.path('RunData', 'BaggedModels')
dir.create(out_dir)

saveRDS(AUCs, file.path(out_dir, paste0(RUN_ID, '_Bagging_AUCs.rds')))

saveRDS(bagged_predictions, 
        file.path(out_dir, paste0(RUN_ID, '_Bagged_predictions.rds')))

saveRDS(cutoff, file.path(out_dir, paste0(RUN_ID, '_Bagging_cutoff.rds')))