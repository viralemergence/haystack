# Notes:
# $(@D) means the directory of the target file


#? Usage: 'make all' updates all files as needed. To re-run the entire analysis, use 'make clean all'.
#? To change configuration options, edit the 'options.config' file in this directory
include options.config
$(info Running on $(N_CORES) cores with random seed $(RANDOM_SEED). Edit options.config to change.)
$(info )


.PHONY: all
all: download \
     get_sequences \
	 get_transcripts \
	 merge_zoonotic_status \
	 merge_and_clean_data \
	 calculate_genomic \
	 train \
	 bag_predictions \
	 make_plots



#?
#? To run individual steps in the pipeline, combine 'make' with the command given in brackets below:

# ----------------------------------------------------------------------------------------
#?	 1. Download external data (download)
# ----------------------------------------------------------------------------------------
EXTERNALDATAFILES = ExternalData/HousekeepingGenes.txt \
					ExternalData/svd_embeddings.csv \
					ExternalData/clover.csv

.PHONY: download
download: $(EXTERNALDATAFILES)


# Gene sets:
ExternalData/HousekeepingGenes.txt:
	mkdir -p ExternalData
	curl -L -o $@ 'https://www.tau.ac.il/~elieis/HKG/HK_genes.txt'

ExternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv:
	curl -L -o $@ 'https://github.com/Nardus/zoonotic_rank/raw/main/InternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv'

ExternalData/Shaw2017_raw/ISG_PublishedData_Web.csv:
	curl -L -o $@ 'https://github.com/Nardus/zoonotic_rank/raw/main/InternalData/Shaw2017_raw/ISG_PublishedData_Web.csv'


# Data from main LF-SVD project:
ExternalData/svd_embeddings.csv:
	mkdir -p ExternalData
	curl -L -o $@ 'https://raw.githubusercontent.com/viralemergence/trefle/main/artifacts/viral_subspace.csv'

ExternalData/clover.csv:
	mkdir -p ExternalData
	curl -L -o $@ 'https://github.com/viralemergence/trefle/data/clover.csv'




# ----------------------------------------------------------------------------------------
#?	 2. Download virus sequences from GenBank (get_sequences)
# ----------------------------------------------------------------------------------------
# Some rules below actually need the individual gb files in ExternalData/Sequences/, 
# but CombinedSequences.fasta is always the last file to be created during download,
# so if it is up to date, all sequences must be present:
ExternalData/Sequences/CombinedSequences.fasta: InternalData/svd_curated_accessions.csv
	python3 Misc/DownloadSequences.py

.PHONY: get_sequences
get_sequences: ExternalData/Sequences/CombinedSequences.fasta



# ----------------------------------------------------------------------------------------
#?	 3. Download human transcript sequences from Ensembl (get_transcripts)
# ----------------------------------------------------------------------------------------
CalculatedData/HumanGeneSets/TranscriptSequences.fasta: ExternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv \
                                                        ExternalData/HousekeepingGenes.txt
	python3 Misc/DownloadGeneSets.py

.PHONY: get_transcripts
get_transcripts: CalculatedData/HumanGeneSets/TranscriptSequences.fasta



# ----------------------------------------------------------------------------------------
#?	 4. Create zoonotic status data (merge_zoonotic_status)
# ----------------------------------------------------------------------------------------
CalculatedData/ZoonoticStatus_Merged.rds: ExternalData/clover.csv
	Rscript Scripts/MergeZoonoticStatusData.R

.PHONY: merge_zoonotic_status
merge_zoonotic_status: CalculatedData/ZoonoticStatus_Merged.rds



# ----------------------------------------------------------------------------------------
#?	 5. Merge and clean final dataset (merge_and_clean_data)
# ----------------------------------------------------------------------------------------
# This has multiple outputs: using a pattern rule ensures the command is
# run just once (see https://www.cmcrossroads.com/article/rules-multiple-outputs-gnu-make)
CalculatedData/FinalData_%.rds CalculatedData/FinalData_%.csv: InternalData/svd_curated_accessions.csv \
															   CalculatedData/ZoonoticStatus_Merged.rds
	Rscript Scripts/MergeAndCleanData.R

.PHONY: merge_and_clean_data
merge_and_clean_data:  CalculatedData/FinalData_Cleaned.rds



# ----------------------------------------------------------------------------------------
#?	 6. Calculate genomic features (calculate_genomic)
# ----------------------------------------------------------------------------------------
# This actually creates multiple output files, but as they are always required
# together, simply ensuring the first of them gets updated
CalculatedData/GenomicFeatures-Virus.rds: CalculatedData/FinalData_Cleaned.rds \
                                          ExternalData/Sequences/CombinedSequences.fasta \
										  ExternalData/Shaw2017_raw/ISG_PublishedData_Web.csv \
										  ExternalData/Shaw2017_raw/ISG_CountsPerMillion_Human.csv \
										  CalculatedData/HumanGeneSets/TranscriptSequences.fasta
	Rscript Scripts/CalculateGenomicFeatures.R


.PHONY: calculate_genomic
calculate_genomic: CalculatedData/GenomicFeatures-Virus.rds



# ----------------------------------------------------------------------------------------
#?	 7. Train models (train)
# ----------------------------------------------------------------------------------------
# NOTE: $(notdir $(@)) means the last part of the target, i.e. 'RunID' in 'RunData/RunID'

# Feature sets /  data required to calculate them:
DIRECT_GENOMIC = CalculatedData/GenomicFeatures-Virus.rds
RELATIVE_GENOMIC = CalculatedData/GenomicFeatures-Distances.rds

# Data common to all possible train calls:
TRAIN_REQUIREMENTS = CalculatedData/FinalData_Cleaned.rds

N_FEATS = 125  # Optimized in preious manuscript (zoonotic_rank)


# Genome features only
RunData/AllGenomeFeatures_SVD: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC) ExternalData/svd_embeddings.csv
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeVirusFeatures --includeISG --includeHousekeeping --includeRemaining \
		--topFeatures $(N_FEATS)

# - SVD runs
RunData/SVD_clover_rank12: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC) ExternalData/svd_embeddings.csv
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeSVD_clover \
		--SVD_max_rank 12 \
		--topFeatures $(N_FEATS)
		
RunData/SVD_trefle_rank12: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC) ExternalData/svd_embeddings.csv
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeSVD_trefle \
		--SVD_max_rank 12 \
		--topFeatures $(N_FEATS)
		
# Combinations with SVD
RunData/AllGenomeFeatures_and_SVD_clover_rank12: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC) ExternalData/svd_embeddings.csv
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeVirusFeatures --includeISG --includeHousekeeping --includeRemaining \
		--includeSVD_clover \
		--SVD_max_rank 12 \
		--topFeatures $(N_FEATS)

RunData/AllGenomeFeatures_and_SVD_trefle_rank12: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC) ExternalData/svd_embeddings.csv
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeVirusFeatures --includeISG --includeHousekeeping --includeRemaining \
		--includeSVD_trefle \
		--SVD_max_rank 12 \
		--topFeatures $(N_FEATS)


## Long run for best model:
RunData/AllGenomeFeatures_and_SVD_trefle_rank12_LongRun: $(TRAIN_REQUIREMENTS) $(RELATIVE_GENOMIC) ExternalData/svd_embeddings.csv
	Rscript Scripts/TrainAndValidate.R $(RANDOM_SEED) $(notdir $(@)) --nthread $(N_CORES) \
		--includeVirusFeatures --includeISG --includeHousekeeping --includeRemaining \
		--includeSVD_trefle \
		--SVD_max_rank 12 \
		--topFeatures $(N_FEATS) \
		--nboot 1000 --nseeds 100


ALL_RUN_IDS = AllGenomeFeatures_SVD SVD_clover_rank12 SVD_trefle_rank12 \
				AllGenomeFeatures_and_SVD_clover_rank12 \
				AllGenomeFeatures_and_SVD_trefle_rank12 \
				AllGenomeFeatures_and_SVD_trefle_rank12_LongRun

TRAIN_OUTPUT_FOLDERS = $(patsubst %, RunData/%, $(ALL_RUN_IDS))

.PHONY: train
train: $(TRAIN_OUTPUT_FOLDERS)



# ----------------------------------------------------------------------------------------
#?	8. Bagged predictions (bag_predictions)
# ----------------------------------------------------------------------------------------
# Currently only using bagging for long runs - need each virus to occur enough test sets:
RunData/AllGenomeFeatures_and_SVD_trefle_rank12_LongRun/AllGenomeFeatures_and_SVD_trefle_rank12_LongRun_Bagged_predictions.rds: | RunData/AllGenomeFeatures_and_SVD_trefle_rank12_LongRun
	Rscript Scripts/CalculateBaggedPredictions.R $(RANDOM_SEED) AllGenomeFeatures_and_SVD_trefle_rank12_LongRun --Ntop 100


.PHONY: bag_predictions
bag_predictions: RunData/AllGenomeFeatures_and_SVD_trefle_rank12_LongRun/AllGenomeFeatures_and_SVD_trefle_rank12_LongRun_Bagged_predictions.rds



# ----------------------------------------------------------------------------------------
#?	9. Plot (make_plots)
# ----------------------------------------------------------------------------------------

# TODO

.PHONY: make_plots
make_plots: 



# ----------------------------------------------------------------------------------------
# Cleanup
# ----------------------------------------------------------------------------------------
.PHONY: confirm clean

confirm:
	@echo -n "Removing generated files - are you sure? [y/N] " && read ans && [ $${ans:-N} = y ]

#?
#? Other commands:
#?	clean: Remove all intermediate files, including those required for predictions (which are distributed)
clean: confirm
	-rm -rfv ExternalData
	-rm -rfv CalculatedData
	-rm -rfv RunData
	-rm -rfv Plots


# ----------------------------------------------------------------------------------------
# Auto document this file 
#  - Comments above that start with a ? become help strings
# ----------------------------------------------------------------------------------------
.PHONY: help
help: Makefile
	@grep "^#?" $< | cut -c4-


# ----------------------------------------------------------------------------------------
# Make options
# ----------------------------------------------------------------------------------------
.DELETE_ON_ERROR:
.SECONDARY:
.NOTPARALLEL:  # Individual scripts are already parallel

#?
#?
