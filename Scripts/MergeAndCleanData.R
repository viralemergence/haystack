#! Rscript
# Combine zoonotic status and accession data and perform final data cleaning steps


ALLOW_INDIRECT_DETECTION <- FALSE  # Whether serological detections should be
																	 # included when determining zoonotic status and known hosts


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Dependencies and data ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(readr)
library(dplyr)
library(tidyr)

options(stringsAsFactors = FALSE)

FinalAccessions <- read_csv("InternalData/svd_curated_accessions.csv", 
                            col_types = cols(.default = "c"))

ZoonoticStatusData <- readRDS("CalculatedData/ZoonoticStatus_Merged.rds")


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Merge data ---------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
FinalData <- FinalAccessions %>% 
  select(-.data$notes) %>% 
  left_join(ZoonoticStatusData, by = c("LatestSppName" = "Virus"))

stopifnot(!any(is.na(FinalData$IsZoonotic)))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Main cleanup -------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Modify zoonotic status of viruses thought to be zoonotic based only on indirect 
# detection methods (if requested):
if (!ALLOW_INDIRECT_DETECTION) {
	FinalData <- FinalData %>% 
		mutate(IsZoonotic = if_else(.data$IndirectDetectionOnly, FALSE, .data$IsZoonotic))
}


## Other cleanup steps:
#		- Remove viruses which have no available sequence information
FinalData <- FinalData %>% 
	filter(! is.na(accession))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output final data --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
saveRDS(FinalData, "CalculatedData/FinalData_Cleaned.rds")
