##
## Utility functions for reading and processing run data files
## 
library(readr)
library(dplyr)

read_run_rds <- function(runid, suffix, run_dir = 'RunData') {
	## Read an rds file from a specific run
	print(runid)
	paste0(runid, suffix) %>% 
		file.path(run_dir, runid, .) %>% 
		read_rds()
}

read_run_csv <- function(runid, suffix, ..., run_dir = 'RunData') {
	## Read a csv file from a specific run
	print(runid)
	paste0(runid, suffix) %>% 
		file.path(run_dir, runid, .) %>% 
		read.csv(stringsAsFactors = FALSE, ...) %>% 
		as_tibble()
}


read_multiple_runs <- function(runids, suffix, file_type = c('rds', 'csv'), simplify = TRUE, ...) {
	# Read data for multiple runs
	#   if simplify = TRUE, returns a single dataframe, with a 'RunID' column added;
	#   otherwise, a named list is returned
	file_type <- match.arg(file_type)
	
	if (file_type == 'rds') {
		res <- sapply(runids, FUN = read_run_rds, suffix = suffix, ..., simplify = FALSE)
	} else {
		res <- sapply(runids, FUN = read_run_csv, suffix = suffix, ..., simplify = FALSE)
	}
	
	if (simplify) {
		bind_rows(res, .id = 'RunID')
	} else {
		res
	}
}

