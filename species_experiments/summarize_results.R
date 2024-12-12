#############################################
# Summarize species results

# December 10, 2024
#############################################

# Read the species description dataframe
species_descr_df <- read.csv("species_descr_ext.csv")

# Format species traits information
species_descr_df[species_descr_df$Prevalence.Level=="l",]$Prevalence.Level <- "low"
species_descr_df[species_descr_df$Prevalence.Level=="m",]$Prevalence.Level <- "medium"
species_descr_df[species_descr_df$Prevalence.Level=="h",]$Prevalence.Level <- "high"
species_descr_df$Prevalence.Level <- factor(species_descr_df$Prevalence.Level, levels = c("low", "medium", "high"))

species_descr_df[species_descr_df$Habitat=="f",]$Habitat <- "forested"
species_descr_df[species_descr_df$Habitat=="e",]$Habitat <- "seral"
species_descr_df$Habitat <- factor(species_descr_df$Habitat, levels = c("forested", "seral"))

species_descr_df[species_descr_df$Generalist.Specialist=="g",]$Generalist.Specialist <- "generalist"
species_descr_df[species_descr_df$Generalist.Specialist=="s",]$Generalist.Specialist <- "specialist"
species_descr_df$Generalist.Specialist <- factor(species_descr_df$Generalist.Specialist, levels = c("generalist", "specialist"))


species_descr_df[species_descr_df$Home.Range=="s",]$Home.Range <- "small"
species_descr_df[species_descr_df$Home.Range=="m",]$Home.Range <- "medium"
species_descr_df[species_descr_df$Home.Range=="l",]$Home.Range <- "large"
species_descr_df$Home.Range <- factor(species_descr_df$Home.Range, levels = c("small", "medium", "large"))



# Initialize results dataframes

best_pars_summarized <- data.frame(
    species = character(0),
    method = character(0),
    par_name = character(0),
    best_par = character(0)
)

results_clust <- data.frame(
    species = character(0),
    method = character(0),
    n_points = numeric(0),
    n_clusters = numeric(0),
    min_size = numeric(0),
    max_size = numeric(0),
    mean_size = numeric(0),
    sd_size = numeric(0),
    perc_svs = numeric(0)
)

results_sum <- data.frame(
    species = character(0),
    prevalence = character(0),
    habitat = character(0),
    specialist = character(0),
    home_range = character(0),
    run = character(0),
    method = character(0),
    auc = numeric(0)
)

data_loss <- data.frame(
    method = character(0),
    original_train_size = numeric(0),
    clustered_on = numeric(0),
    thrown = numeric(0),
    percentage_thrown = numeric(0)
)

placeholder_spec_name <- "AMCR"

# Loop over species
for (spec_name in species_names) {

    
    # Get species description
    species_descr <- species_descr_df[species_descr_df$Abbreviation == spec_name,]
    

    prevalence <- species_descr$Prevalence.Level
    habitat <- species_descr$Habitat
    specialist <- species_descr$Generalist.Specialist
    home_range <- species_descr$Home.Range
    
    # Loop over clustering methods
    for (m_name in method_names) {
        
        method_name <- parseMethodName(m_name)

        if (startsWith(method_name, "BayesOpt")){
            
            bo_pars <- read.delim(paste0("results/species_experiments/raw/clustering_parameters/", spec_name,"_", method_name,".csv"), sep = ",", header = T)
            n_pars <- nrow(bo_pars)

            t_df <- data.frame(
                species = rep(spec_name, n_pars),
                method = rep(method_name, n_pars),
                par_name = bo_pars$parameter,
                best_par = bo_pars$value
            )
            best_pars_summarized <- rbind(best_pars_summarized, t_df)

        }

        clust_metrics <- read.delim(paste0("results/species_experiments/raw/clusterings/", spec_name,"_", method_name, "_clust_descr_stats.csv"), sep = ",", header = T)
        
        t_df <- data.frame(
                species = spec_name,
                method = method_name,
                n_points = clust_metrics$n_points,
                n_clusters = clust_metrics$n_clusters,
                min_size = clust_metrics$min_size,
                max_size = clust_metrics$max_size,
                mean_size = clust_metrics$mean_size,
                sd_size = clust_metrics$sd_size,
                perc_svs = clust_metrics$perc_svs
            )
        results_clust <- rbind(results_clust, t_df)

        # Loop over experimental runs/repeats
        for (exp_run in seq(1:runs)){

            perf_metrics <- read.delim(paste0("results/species_experiments/raw/metrics/", spec_name,"_", method_name, "_run=", exp_run, ".csv"), sep = ",", header = T)

            t_df <- data.frame(
                species = spec_name,
                prevalence = prevalence,
                habitat = habitat,
                specialist = specialist,
                home_range = home_range,
                run = exp_run,
                method = method_name,
                auc = perf_metrics$auc,
                auprc = perf_metrics$auprc
            )

            results_sum <- rbind(results_sum, t_df)
        }


        # Summarize data loss stats

        train.df.og <- read.delim(paste0("checklist_data/species/", spec_name, "/", spec_name, "_zf_filtered_region_2017.csv"), sep = ",", header = T)
        
        train.df.og <- train.df.og[!is.na(train.df.og$duration_minutes),]

        train.df.og <- train.df.og[train.df.og$observation_date >= "2017-05-15" & train.df.og$observation_date <= "2017-07-09",]

        og_train_size <- nrow(train.df.og)

        clustering <- read.delim(paste0("results/species_experiments/raw/clusterings/", spec_name,"_", method_name,".csv"), sep = ",", header = T)
        clust_on <- nrow(clustering)


        if ((clust_on < og_train_size) && (spec_name == placeholder_spec_name)){
            t_df <- data.frame(
                method = method_name,
                original_train_size = og_train_size,
                clustered_on = clust_on,
                thrown = og_train_size - clust_on,
                percentage_thrown = round((((og_train_size-clust_on)/og_train_size) * 100), 4)
                
            )
            data_loss <- rbind(data_loss, t_df)
        }
    }

    

}


# Write the summarized results to CSV
write.csv(best_pars_summarized, "results/species_experiments/summarized/best_pars_summarized.csv", row.names=FALSE)
write.csv(results_clust, "results/species_experiments/summarized/results_clust_summarized.csv", row.names=FALSE)
write.csv(results_sum, "results/species_experiments/summarized/results_summarized.csv", row.names=FALSE)
write.csv(data_loss, "results/species_experiments/summarized/data_loss_summarized.csv", row.names=FALSE)
