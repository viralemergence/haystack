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

HostData <- read_csv("ExternalData/clover.csv")


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
#   - Remove viruses whose only host is human - since embeddings are calculated while excluding
#     humans, these viruses become the only included entries with no connections in the network,
#     resulting in a data leak if included
human_only <- HostData %>% 
  group_by(.data$Virus) %>% 
  summarise(n_nonhuman = sum(.data$Host != "Homo sapiens")) %>% 
  filter(n_nonhuman == 0) %>% 
  pull(.data$Virus)

FinalData <- FinalData %>% 
	filter(!is.na(accession)) %>% 
  filter(!.data$LatestSppName %in% human_only)


## Rename columns to match downstream scripts
FinalData <- FinalData %>% 
  rename(InfectsHumans = IsZoonotic,
         Accessions = accession)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output final data --------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
saveRDS(FinalData, "CalculatedData/FinalData_Cleaned.rds")
