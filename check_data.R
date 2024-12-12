############################################
# Check prevalence rates and dataset sizes

# December 10, 2024
############################################

library(here)
library(dplyr)

# set working directory to current/source directory
setwd(here::here())



spec_names <- c(
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

for (spec_name in spec_names) {

    cat("Species: " , spec_name, "\n")
    f.name_train <- paste0(spec_name, "/", spec_name, "_zf_filtered_region_2017.csv")
    f.name_test <- paste0(spec_name, "/", spec_name, "_zf_filtered_region_2018.csv")

    train.df.og <- read.delim(paste0("checklist_data/species/", f.name_train), sep = ",", header = T)
    test.df.og <- read.delim(paste0("checklist_data/species/", f.name_test), sep = ",", header = T)

    cat("Train checklists: ", nrow(train.df.og), "; Test checklists: ", nrow(test.df.og),"\n")
    
    train.df.og <- train.df.og[!is.na(train.df.og$duration_minutes),]
    test.df.og <- test.df.og[!is.na(test.df.og$duration_minutes),]

    cat("After removing checklists with missing duration_minutes => Train checklists: ", nrow(train.df.og), "; Test checklists: ", nrow(test.df.og), "\n")
    
    train.df.og <- train.df.og[train.df.og$observation_date >= "2017-05-15" & train.df.og$observation_date <= "2017-07-09",]
    test.df.og <- test.df.og[test.df.og$observation_date >= "2018-05-15" & test.df.og$observation_date <= "2018-07-09",]
    
    cat("After filtering by observation_date -> Train checklists: ", nrow(train.df.og), "; Test checklists: ", nrow(test.df.og), "\n")

    train.df.og.c <- train.df.og
    train.df.og.c$lat_long <- paste0(train.df.og.c$latitude, "_", train.df.og.c$longitude)
    uniq_loc_df_train <- dplyr::distinct(train.df.og.c, lat_long, .keep_all = T)

    test.df.og.c <- test.df.og
    test.df.og.c$lat_long <- paste0(test.df.og.c$latitude, "_", test.df.og.c$longitude)
    uniq_loc_df_test <- dplyr::distinct(test.df.og.c, lat_long, .keep_all = T)
    

    train_ul = nrow(uniq_loc_df_train)
    test_ul = nrow(uniq_loc_df_test)
    train_prevalence = round(mean(train.df.og$species_observed)*100, 4)
    test_prevalence = round(mean(test.df.og$species_observed)*100, 4)
    prevalence = round(mean(c(train_prevalence, test_prevalence)), 4)
    
    cat("Train unique locations: ", train_ul, "; Test unique locations: ", test_ul, "; Train prevalence: ", train_prevalence, "%; Test Prevalence: ", test_prevalence, "%; Overall Prevalence: ", prevalence,"%\n")


    # write.csv(train.df.og, paste0("checklist_data/species/", spec_name, "/2017_final.csv"), row.names=FALSE)
    # write.csv(test.df.og, paste0("checklist_data/species/", spec_name, "/2018_final.csv"), row.names=FALSE)


    cat("\n")    

}
