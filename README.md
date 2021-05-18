# Zoonotic risk assessment from viral genomes

Code used for the zoonotic potential part of Poisot _et al._ (2021) "Imputing the mammalian virome with linear filtering and singular value decomposition". This repository is derived from [zoonotic_rank](https://github.com/Nardus/zoonotic_rank), described in [Mollentze et al. (2020)](https://doi.org/10.1101/2020.11.12.379917).



## Table of contents
- [Requirements](#requirements)
- [Repeating all analyses](#repeating-all-analyses)
- [File structure](#file-structure)



## Requirements
- [R](https://www.r-project.org/) (tested using version 3.5.1)
   - Most required R libraries can be installed using `Rscript -e "renv::restore()"`
   - Install `ggtree` from bioconductor using: `Rscript -e "install.packages('BiocManager'); BiocManager::install('ggtree')"`
- [Python](https://www.python.org/) (version >=3.6)
   - [Biopython](https://biopython.org/)
   - [Pandas](https://pandas.pydata.org/)
   - [xlrd](https://xlrd.readthedocs.io/en/latest/)
- [Java JDK](https://www.oracle.com/uk/java/technologies/javase-downloads.html) (version >=8)



## Usage

Follow instructions below to repeat the analyses described in the manuscript. Note that making predictions for novel viruses is not currently possible with these models; but can be done with the original genome feature based models (see [zoonotic_rank](https://github.com/Nardus/zoonotic_rank)).

#### Basic
These steps will download any missing source data and automatically create/update files as needed.

_Using Rstudio:_
1. Open `haystack_zoonotic.Rproj` in RStudio
2. On the `Build` tab, select `More` > `Clean and Rebuild`

_Using the command-line:_
```
make clean all
```

#### Advanced options (command-line only)

- Use `make help` to see individual steps in the pipeline. Upstream steps are run automatically if needed. For example, using `make prepare` will run the data cleanup step, but also downloads the raw data if needed.
- `make <path to file>` runs all steps neccesary to produce/update the specified file (e.g. `make Plots/Figure1.pdf`).
- `make as_distributed` resets the project to the state in which it was distributed.
- `make clean` removes all run-related files, allowing a complete re-run (in contrast to `as_distributed`, this includes removing the pre-trained models required for prediction).


## File structure
(files/folders which will be created during a full run are indicated by `[*]`)

```
└─zoonotic_rank/
   ├─Makefile ................................. Record of workflow and dependencies
   │                                            between files
   ├─options.config ........................... Runtime options (speciefies number
   │                                            of parrallel threads allowed and 
   │                                            the random seed)
   ├─InternalData/ ............................ All data unique to this project
   │   ├─example_files/ ....................... Example input files for 
   │   │                                        predicting novel viruses 
   │   ├─Shaw2017_raw/ ........................ Raw ISG data from Shaw et al. 
   │   │                                        2017 (see https://isg.data.cvr.ac.uk/)
   │   ├─FinalData_Cleaned.csv .................Final dataset, as used for training. 
   │   │                                        Created by merging files below 
   │   │                                        (see Scripts/MergeAndCleanData.R)
   │   ├─AllInternalData_Checked.csv .......... Metadata for the viruses used 
   │   │                                        as training data
   │   ├─Final_Accessions_Unique_Spp.csv ...... Accession numbers of sequences 
   │   │                                        used for training (replaces 
   │   │                                        those in the metadata file)
   │   ├─NameMatches_All.csv .................. Manually curated list used to 
   │   │                                        match virus names to unique 
   │   │                                        species across external datasets
   │   ├─SourcesOfZoonoses_BabayanZoonotic.csv  Additional zoonotic status data 
   │   │                                        for species not available in 
   │   │                                        external data sources
   │   └─Taxonomy_UnclassifiedViruses.csv ..... Taxonomic information for 
   │                                            unclassified viruses in the 
   │                                            metadata (unused)
   │
   ├─CalculatedData/ .......................... Intermediate calculations ([*], except 
   │                                            for files required by PredictNovel.R)
   ├─ExternalData/ ............................ [*] Data from external sources, 
   │                                            dowloaded as needed (see Makefile)
   ├─Misc/ .................................... Miscelaneous scripts to download 
   │                                            external data
   ├─Plots/ ................................... [*] Final plots generated
   ├─Predictions/ ............................. [*] Predictions for case studies
   ├─renv/ .................................... Record of R libraries required
   ├─RunData/ ................................. Trained models ([*], except for 
   │                                            files required by PredictNovel.R)
   ├─Scripts/ ................................. Main analysis, prediction, and 
   │   │                                        plotting scripts
   │   └─Plotting/ ............................ Scripts to generate published plots
   ├─Tests/ ................................... Unit tests for basic functionality 
   │                                            of utility scripts
   └─Utils/ ................................... Utility functions and tools called 
                                                by other scripts
```
