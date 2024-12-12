##############################################
# Preprocess training and test data

# December 10, 2024
##############################################

library(dggridR) # spatial subsampling operations


dirs_to_create = c(
    "results",
    "results/species_experiments",
    "results/species_experiments/raw",
    "results/species_experiments/raw/test_data_splits",
    "results/species_experiments/raw/clusterings",
    "results/species_experiments/raw/clustering_parameters",
    "results/species_experiments/raw/model_parameters",
    "results/species_experiments/raw/predictions",
    "results/species_experiments/raw/metrics",
    "results/species_experiments/summarized",
    "results/species_experiments/plotted",
    "results/species_experiments/plotted/plots",
    "results/species_experiments/plotted/maps"
)

for (f.path in dirs_to_create){
    if (!dir.exists(f.path)) {
        dir.create(f.path)
    }
}



data_path <- "checklist_data/species/"

spec_names <- list.files(path = data_path, include.dirs = TRUE)

for (spec_name in spec_names){

    # For code reproducibility of test splits
    set.seed(1)

    if (!dir.exists(paste0("results/species_experiments/raw/test_data_splits/", spec_name))) {
        dir.create(paste0("results/species_experiments/raw/test_data_splits/", spec_name))
    }  

    # f.name_train <- paste0(spec_name, "/", spec_name, "_zf_filtered_region_2017.csv")
    f.name_test <- paste0(spec_name, "/", spec_name, "_zf_filtered_region_2018.csv")

    # train.df.og <- read.delim(paste0(data_path, f.name_train), sep = ",", header = T)
    test.df.og <- read.delim(paste0(data_path, f.name_test), sep = ",", header = T)

    # train.df.og <- train.df.og[!is.na(train.df.og$duration_minutes),]
    test.df.og <- test.df.og[!is.na(test.df.og$duration_minutes),]

    # train.df.og <- train.df.og[train.df.og$observation_date >= "2017-05-15" & train.df.og$observation_date <= "2017-07-09",]
    test.df.og <- test.df.og[test.df.og$observation_date >= "2018-05-15" & test.df.og$observation_date <= "2018-07-09",]


    for(exp_run in seq(1:runs)){

        test.df <- test.df.og
        
        ########
        # construct validation dataset as specified by this paper:
        # https://onlinelibrary.wiley.com/doi/epdf/10.1111/ddi.13271
        ########
        hexagons <- dgconstruct(spacing = 5, topology = "HEXAGON")
        test.df$cell <- dgGEO_to_SEQNUM(hexagons, test.df$latitude, test.df$longitude)$seqnum
        
        test.det.df <- test.df[test.df$species_observed == T,]
        test.undet.df <- test.df[test.df$species_observed == F,]
        
        # subsampling detected checklists
        cell_names <- names(table(test.det.df$cell))
        det.valid.df <- spatial.subsample(test.det.df, cell_names)
        
        # subsampling undetected checklists
        undet.cell_names <- names(table(test.undet.df$cell))
        undet.valid.df <- spatial.subsample(test.undet.df, undet.cell_names)
        
        # downsampling to the df with fewer checklists
        if(nrow(undet.valid.df) > nrow(det.valid.df)){
            idx <- sample(seq(1:nrow(undet.valid.df)), nrow(det.valid.df))
            undet.valid.df <- undet.valid.df[idx,]
        } else if (nrow(undet.valid.df) < nrow(det.valid.df)){
            idx <- sample(seq(1:nrow(det.valid.df)), nrow(undet.valid.df))
            det.valid.df <- det.valid.df[idx,]
            
        }
        test.df <- rbind(det.valid.df, undet.valid.df)


        write.csv(test.df, paste0("results/species_experiments/raw/test_data_splits/", spec_name, "/run=", exp_run,".csv"), row.names=FALSE)
    }
}