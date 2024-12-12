##############################################
# Main file for running species experiments

# Nahian Ahmed
# December 10, 2024
##############################################


library(here)

# set working directory to current/source directory
setwd(here::here())

# Import required files
source("helper/BayesOptHelper.R")
source("helper/clustGeoHelper.R")
source("helper/DBSCHelper.R")
source("helper/helpers.R")
source("helper/kmsq.R")
source("helper/mainHelper.R")



#######
# Main function for running site clustering experiments
#######
run_site_clustering_experiments = function(species_names, method_names, runs, preprocess=TRUE, experiments=TRUE, summarize=TRUE, plot_results=TRUE){

    assign("species_names", species_names, envir = .GlobalEnv)
    assign("method_names", method_names, envir = .GlobalEnv)


    if (preprocess){
        source("species_experiments/preprocess_data.R")
    }


    if (experiments){

        # Run experiments
        for (species_name in species_names){
            
            assign("spec_name", species_name, envir = .GlobalEnv)
            
            # Run experiments on specific species data
            set.seed(1) # Ensures that order of species in species_names does not affect reproducibility
            source("species_experiments/run_experiments.R")

        }
    }

    if (summarize){
        # Summarize results
        source("species_experiments/summarize_results.R")
    }

    if (plot_results){
        # Plot results
        source("species_experiments/plot_results.R")
    }
}




#######
# Specify species data to run experiments on
# Please see species_descr_ext.csv for full names and other information
#######
species_names <- c(
    "AMCR",
    "AMRO",
    "BAEA",
    "BKHGRO",
    "BRCR",
    "BUTI",
    "CASC",
    "CHBCHI",
    "COHA",
    "HAFL",
    "HAWO",
    "HEWA",
    "MAWA",
    "MOQU",
    "NOFL",
    "NOOW",
    "OLFL",
    "PAFL",
    "PAWR",
    "PIWO",
    "REHA",
    "SOSP",
    "SPTO",
    "SWTH",
    "WAVI",
    "WEPE",
    "WETA",
    "WIWA",
    "WRENTI",
    "YEBCHA",
    "YEWA"
)

#######
# Specify site clustering methods
# Please see helper/mainHelper.R for method strings
#######
method_names <- c(
    "two_to_10",               #2to10 
    "two_to_10_sameObs",       #2to10-sameObs
    "kmSq-1000",               #1-kmSq
    "lat_long",                #lat-long
    "rounded-4",            
    "svs",                     #SVS
    "one_UL",                  #1-UL
    "clustGeo_25_60",
    "clustGeo_50_60",
    "clustGeo_75_60",
    "clustGeo_25_70",
    "clustGeo_50_70",
    "clustGeo_75_70",
    "clustGeo_25_80",
    "clustGeo_50_80",
    "clustGeo_75_80",
    "clustGeo_25_90",
    "clustGeo_50_90",
    "clustGeo_75_90",
    "DBSC",
    "BayesOptClustGeo"
)


#######
# Specify number of experiment runs/repeats
#######
runs <- 25


#######
# Run site clustering experiments
#######
run_site_clustering_experiments(
    species_names = species_names,
    method_names = method_names,
    runs = runs,
    preprocess = TRUE,
    experiments = TRUE,
    summarize = TRUE,
    plot_results = TRUE
)
