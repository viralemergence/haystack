# Use imputed host-virus associations to predict human infection risk

Code used for the zoonotic potential part of Poisot _et al._ (2021) "Imputing the mammalian virome with linear filtering and singular value decomposition". See [viralemergence/trefle](https://github.com/viralemergence/trefle) for the main code of this paper.

This repository is derived from [zoonotic_rank](https://github.com/Nardus/zoonotic_rank), described in [Mollentze et al. (2020)](https://doi.org/10.1101/2020.11.12.379917).


## Requirements
- Install the [conda package manager](https://conda.io/)
- All other requirements can be installed in a self-contained environment using:
```
conda env create -f environment.yml
```


## Usage

The commands below will repeat all analyses described in the manuscript. Note that making predictions for novel viruses is not currently possible with these models, but can be done with the original genome feature-based models (see [zoonotic_rank](https://github.com/Nardus/zoonotic_rank)). Running the full pipeline will take several days (approximately 1 week on a 16-core AMD 5950X processor). 
```
conda activate haystack_zoonotic
make
```

Use `make help` to see individual steps in the pipeline and further instructions. 


## Folder structure
(folders which will be created during model training are indicated by `[*]`)

```
└─haystack_zoonotic/
   ├─Makefile ................................. Workflow specification
   ├─environment.yml .......................... Record of software/package versions 
   │                                            used
   ├─options.config ........................... Runtime options (specifies number
   │                                            of parallel threads allowed and 
   │                                            the random seed)
   ├─InternalData/ ............................ Data unique to this project
   │   └─svd_curated_accessions.csv ........... Accession numbers for each virus  
   │                                            species in the clover database
   │                                            (where a full genome is available)
   ├─CalculatedData/ .......................... [*] Intermediate calculations 
   ├─ExternalData/ ............................ [*] Data from external sources, 
   │                                            dowloaded as needed (see Makefile)
   ├─Misc/ .................................... Scripts to download external data
   ├─Plots/ ................................... [*] Final plots generated
   ├─RunData/ ................................. [*] Trained models
   ├─Scripts/ ................................. Main analysis and plotting scripts
   └─Utils/ ................................... Utility functions and tools called 
                                                by other scripts, including CPB_Machine.jar
                                                from https://github.com/rjorton/VirusFeatures
```
