##
## Part of Zoonosis prediction pipeline:
## 	- Process input options passed to training script, loading optional feature data as needed
## 	
## 

library(argparse)

parser <- ArgumentParser(description = 'Train a model with specific feature sets. Validation is against subsets of the test set only - this script does not see the holdout data.')
parser$add_argument('RandomSeed', type = 'integer', 
										help = 'a random seed (used by both R and xgboost)')

parser$add_argument('RunID', type = 'character', 
										help = 'an identifier for this run, used as a base to name output files')


## Optional arguments:
# Feature types to include in training:
featureGroup <- parser$add_argument_group('Feature sets')
featureGroup$add_argument('--includeVirusFeatures', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include direct genomic features (calculated directly from virus genomes)')
featureGroup$add_argument('--includeISG', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include genomic features relative to human ISGs')
featureGroup$add_argument('--includeHousekeeping', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include genomic features relative to human housekeeping genes')
featureGroup$add_argument('--includeRemaining', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include genomic features relative to remaining human genes')
featureGroup$add_argument('--includeSVD_clover', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include SVD embeddings of observed host range (i.e. the clover database)')
featureGroup$add_argument('--includeSVD_trefle', action = 'store_const', const = TRUE, default = FALSE,
													help = 'include SVD embeddings of imputed host range (i.e. clover plus imputations)')
featureGroup$add_argument('--SVD_max_rank', type = 'integer', default = 12, 
													help = 'maximum rank to include for SVD features')


# Options to control training
trainGroup <- parser$add_argument_group('Training options')
trainGroup$add_argument('--trainProportion', metavar = 'p', type = 'double', default = 0.7, 
												help = 'proportion of virus species to select for each training set - all remaining virus species go to test set (default: 0.85)')
trainGroup$add_argument('--topFeatures', metavar = 'f', type = 'double', default = 100, 
												help = 'Number of features to retain for final training')
trainGroup$add_argument('--nboot', metavar = 'b', type = 'integer', default = 100, 
												help = 'number of iterations to perform (default: 100)')
trainGroup$add_argument('--nseeds', metavar = 's', type = 'integer', default = 10, 
												help = 'number of xgboost random seeds to use (default: 10). If e.g. nboot = 100 and nseeds = 10, there will be 10 iterations using each seed.')

trainGroup$add_argument('--nthread', metavar = 't', type = 'integer', default = 1, 
												help = 'number of parallel threads allowed (default: 1)')


# Parse input:
# - These may be specified by the calling script (as a vector of strings named OVERRIDE_EXTERNAL_COMMANDS),
#   but will be taken from arguments used to invoke the calling script if this is left unspecified
if (exists("OVERRIDE_EXTERNAL_COMMANDS")) {
	INPUT <- parser$parse_args(OVERRIDE_EXTERNAL_COMMANDS)
} else {
	INPUT <- parser$parse_args()
}

if (INPUT$trainProportion <= 0 | INPUT$trainProportion >= 1) 
	stop('trainProportion should be between 0 and 1') # But exactly 0 or 1 makes no sense either


if (! any(INPUT$includeVirusFeatures, INPUT$includeISG, 
					INPUT$includeHousekeeping, INPUT$includeRemaining, 
					INPUT$includeSVD_clover, INPUT$includeSVD_trefle))
	stop('At least one of the feature set flags must be set. See TrainAndValidate.R --help')


if (INPUT$nboot <= 0)
	stop('nboot must be positive')

if (INPUT$nseeds <= 0)
	stop('nseeds must be positive')

if (INPUT$nseeds > INPUT$nboot)
	stop('nseeds must be less than or equal to nboot')

if (INPUT$nboot %% INPUT$nseeds != 0)
	warning('nseeds is not a multiple of nboot. nboot will be increased')


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Core data ----------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
set.seed(INPUT$RandomSeed)

InputData <- readRDS(file.path('CalculatedData', 'FinalData_Cleaned.rds'))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Load/Add optional data based on input options ----------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
FinalData <- InputData  # Will be modified below


if (INPUT$includeVirusFeatures) {
	# Virus features (created by CalculateGenomicFeatures.R)
	VirusFeaturesDirect <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Virus.rds'))
	
	VirusFeaturesDirect <- VirusFeaturesDirect %>% 
		rename_at(vars(-LatestSppName), ~ paste('VirusDirect', ., sep = '_'))
	
	FinalData <- FinalData %>% 
		left_join(VirusFeaturesDirect, by = 'LatestSppName')
}


if (INPUT$includeISG | INPUT$includeHousekeeping | INPUT$includeRemaining) {
	# Genomic features (created by CalculateGenomicFeatures.R)
	GenomicDistances <- readRDS(file.path('CalculatedData', 'GenomicFeatures-Distances.rds'))
	
	if (!INPUT$includeISG) GenomicDistances <- select(GenomicDistances, -contains('_ISG_'))
	if (!INPUT$includeHousekeeping) GenomicDistances <- select(GenomicDistances, -contains('_Housekeeping_'))
	if (!INPUT$includeRemaining) GenomicDistances <- select(GenomicDistances, -contains('_Remaining_'))
	
	FinalData <- FinalData %>% 
		left_join(GenomicDistances, by = 'LatestSppName')
}


if (INPUT$includeSVD_clover | INPUT$includeSVD_trefle) {
	raw_embeddings <- read.csv("ExternalData/svd_embeddings.csv", stringsAsFactors = FALSE)
	
	raw_embeddings <- raw_embeddings %>% 
		filter(rank <= INPUT$SVD_max_rank)
	
	if (INPUT$includeSVD_clover) {
		clover_feats <- raw_embeddings %>% 
			select(virus, rank, clover_nohuman) %>% 
			mutate(rank = paste0("svd_clover_", rank)) %>% 
			spread(key = rank, value = clover_nohuman)
		
		FinalData <- FinalData %>% 
			left_join(clover_feats, by = c("LatestSppName" = "virus"))
	}
	
	if (INPUT$includeSVD_trefle) {
		trefle_feats <- raw_embeddings %>% 
			select(virus, rank, trefle_nohuman) %>% 
			mutate(rank = paste0("svd_trefle_", rank)) %>% 
			spread(key = rank, value = trefle_nohuman)
		
		FinalData <- FinalData %>% 
			left_join(trefle_feats, by = c("LatestSppName" = "virus"))
	}
}

