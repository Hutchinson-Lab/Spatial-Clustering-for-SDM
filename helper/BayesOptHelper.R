###########################
# BayesOptClustGeo helper

# December 10, 2024
###########################

library(rBayesianOptimization) # bayesian optimization library
library(cluster) # silhouette width calculation




train_data <- NA
occ_covs <- NA
det_covs <- NA
fit_func <- NA


normalize <- function(x) {
    return((x- min(x)) /(max(x)-min(x)))
}

get_pairwise_distances <- function(df, sp_features, env_features, normalize = TRUE) {
    
    
    if (normalize){
        for (sp_feature in sp_features){
            df[,sp_feature] <- normalize(df[,sp_feature])
        }
    }
    df <- df[,c(sp_features,env_features)]
    


    m_dist <- distances::distances(df)
    
    return(m_dist)
}

# Fitness function for clustGeo
clustGeo_fit <- function(alpha, lambda){


    og_data_cgul <- train_data

    og_data_cgul$lat_long <- paste0(og_data_cgul$latitude, "_", og_data_cgul$longitude)
    uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
    
    
    percent <- lambda/100
    
    clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
    
    # link un labeled, but lat-long duplicated checklists to the correct site
    og_data_cgul$site <- -1
    for(j in seq(1:nrow(clustGeo_df_i))){
        og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
    }

    # Calculate Silhouette width
    m_dist <- get_pairwise_distances(og_data_cgul, c("latitude","longitude"), occ_covs)
    m_sil <- silhouette(og_data_cgul$site, m_dist)

    silh <- mean(m_sil[, 3]) # third column is silhoutte widths of instances
    if (fit_func == "silhouette"){ # Silhouette width

        result <- list(Score = silh, 
                Pred = 0)

        return(result)
    }
    
    else if (fit_func == "prevalence"){ # sum of "silhouette" and "prevalence"

        n_sites <- length(unique(og_data_cgul$site))
        n_occ_sites <- og_data_cgul %>%
            group_by(site) %>%
            summarize(has_true = any(species_observed == TRUE)) %>%
            filter(has_true) %>%
            nrow()

        site_occ_rate <- n_occ_sites/n_sites

        
        prevalence <- mean(og_data_cgul$species_observed)

        sp_objective <- -abs(prevalence - site_occ_rate)  

        # Fitness function = -|(Prevalence) - (Site Occupancy Rate)| + Average Silhouette width
        result <- list(Score =  sp_objective + silh,
                Pred = 0)

        return(result)
    }
    
}



bayesianOptimizedClustGeo <- function(train_data_p, occ_covs_p, det_covs_p, fit_func_p){
    
    train_data <<- train_data_p
    occ_covs <<- occ_covs_p
    det_covs <<- det_covs_p
    fit_func <<- fit_func_p



    
    n_init <- 20 # Number of initial samples
    n_iter <- 10 # Number of iterations

 

    # Define the search boundary
    search_bound <- list(alpha = c(0.01, 0.99),
                           lambda = c(10, 90))


    # Define initial search sample
    # set.seed(1)
    search_grid <- data.frame(alpha = runif(n_init, search_bound$alpha[1], search_bound$alpha[2]),
                                lambda = runif(n_init, search_bound$lambda[1], search_bound$lambda[2]))
    
    
    # Bayesian Optimization
    bayesianOptimized <- rBayesianOptimization::BayesianOptimization(FUN = clustGeo_fit, bounds = search_bound, 
                     init_points = 0, init_grid_dt = search_grid, 
                     n_iter = n_iter, acq = "ucb")
    

    return (list(Best_Value = bayesianOptimized$Best_Value,
                 Best_Pars = list(alpha = bayesianOptimized$Best_Par[1], 
                                  lambda = bayesianOptimized$Best_Par[2])
                )
            )
}