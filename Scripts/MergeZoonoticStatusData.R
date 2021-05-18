##
## - Extract zoonotic status from the clover dataset (https://github.com/viralemergence/clover)
##

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Dependencies and data ----------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
library(readr)
library(dplyr)
library(tidyr)

options(stringsAsFactors = FALSE)

clover <- read_csv("ExternalData/clover.csv",
                   col_types = cols(
                     .default = col_character(),
                     Year = col_double(),
                     Detection_NotSpecified = col_logical(),
                     Detection_Serology = col_logical(),
                     Detection_Genetic = col_logical(),
                     Detection_Isolation = col_logical(),
                     Host_NCBIResolved = col_logical(),
                     Virus_NCBIResolved = col_logical()
                   ))


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Extract zoonotic status --------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
zoonoses_direct <- clover %>% 
  filter(.data$Detection_Genetic | .data$Detection_Isolation) %>% 
  group_by(.data$Virus) %>% 
  filter("Homo sapiens" %in% .data$Host) %>% 
  pull(.data$Virus) %>% 
  unique()

zoonoses_indirect_only <- clover %>% 
  filter(.data$Detection_Serology | .data$Detection_NotSpecified) %>% 
  filter(!.data$Virus %in% zoonoses_direct) %>% 
  group_by(.data$Virus) %>% 
  filter("Homo sapiens" %in% .data$Host) %>% 
  pull(.data$Virus) %>% 
  unique()

zoonotic_status <- clover %>% 
  group_by(.data$Virus) %>% 
  summarise(IsZoonotic = unique(.data$Virus) %in% c(zoonoses_direct, zoonoses_indirect_only),
            IndirectDetectionOnly = unique(.data$Virus) %in% zoonoses_indirect_only)


# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ---- Output -------------------------------------------------------------------------------------
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.create("CalculatedData")
saveRDS(zoonotic_status, "CalculatedData/ZoonoticStatus_Merged.rds")
