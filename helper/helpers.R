####################
# Helper functions

# December 10, 2024
####################

library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(sf) # spatial geometry operations
library(auk) # ebird data processing
library(unmarked) # occupancy modeling
library(rje) # expit function
library(terra) # geospatial operations




#########
# Rounding Lat/Long
#########
roundLatLong <- function(df, rounding_degree){
    df$rounded_lat <- round(df$latitude, digits = rounding_degree)
    df$rounded_long <- round(df$longitude, digits = rounding_degree)
    df$rounded_locality_id <- paste(as.character(df$rounded_long), as.character(df$rounded_lat), sep = "_")
    
    return(df)
}


##########
# site closure for Occ Model:
#     1. constant site covariates
#     2. no false positives (detected only if occupied)
##########
enforceClosure <- function(sites_df, occ_cov_list, sites_list){
    j<-1
    closed_df <- NA

    # print(dim(sites_df))
    # print(length(sites_list))
    for(eBird_site in sites_list){
        
    
        checklists_at_site <- sites_df[sites_df$site == eBird_site,]
        
        for(occCov_i in occ_cov_list){
            checklists_at_site[occCov_i] <- mean(checklists_at_site[[occCov_i]])
        
        }

        
        if(j==1){
            closed_df = checklists_at_site
        } else {
            closed_df = rbind(closed_df, checklists_at_site)
        }
        j = j+1
    
    }
    return(closed_df)
}



calcDescriptiveClusteringStats <- function(clustered_df){
      
    # Calculate clustering and species specific stats
    clust_freq_df <- clustered_df %>% group_by(site) %>% dplyr::summarise(freq=n())

    descr_stats <- list(
        n_points = nrow(clustered_df),
        n_clusters = nrow(clust_freq_df),
        min_size = min(clust_freq_df$freq),
        max_size = max(clust_freq_df$freq),
        mean_size = round(mean(clust_freq_df$freq), 4),
        sd_size = round(sd(clust_freq_df$freq), 4),
        perc_svs = round((nrow(clust_freq_df[clust_freq_df$freq==1,])/nrow(clust_freq_df)) *100, 4)
     )
    return(descr_stats)

}




#######
# extract environmental features
# at checklist locations
#######
extractEnvFeat <- function(df, OR.tif, obs_covs) {
    # Convert the dataframe to a SpatVector
    df.pts <- vect(df, geom = c("longitude", "latitude"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    # Extract environmental features
    env_vars.df <- data.frame(
        checklist_id = df$checklist_id,
        terra::extract(OR.tif, df.pts)
    )
    return(env_vars.df)
}


#######
# normalize dataset
#######
norm_ds <- function(df, det_covs, occ_covs, test=FALSE, norm.list=list()){

    if(length(norm.list) == 0){
        for(name in c(det_covs, occ_covs)){
            # calc mean/var for each cov, if training
            ma <- max(df[[name]])
            mi <- min(df[[name]])
            norm.list[[name]] <- c(ma, mi)
        }
    }
    
    # xi - min(x)/(max(x) - min(x))
    for(cov in names(norm.list)){
        df[[cov]] <- (df[[cov]] - norm.list[[cov]][[2]])/(norm.list[[cov]][[1]] - norm.list[[cov]][[2]])
    }
    
    return(list(df=df, n_l=norm.list))
    
}

#######
# spatial subsampling as defined by: 
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/ddi.13271
#######
spatial.subsample <- function(df, cell.names){
    valid.df <- data.frame()
    i <- 0
    for(freq in table(df$cell)){ 
        i <- i + 1
        if(freq > 1){
            chklsts <- df[df$cell == cell.names[i],]
            samp <- chklsts[sample(seq(1:nrow(chklsts)), 1),]
        } else {
            samp <- df[df$cell == cell.names[i],]
        }
        valid.df <- rbind(valid.df, samp)
    }    
    return(valid.df)
}



#######
# calculates the occupancy model from a given dataset
# containing checklists and sites
#######
# 1. enforces closure
# 2. formats it w/r/t eBird data
# 3. runs through occupancy model
#######
calcOccModel <- function(df, occ_covs, det_covs, skip_closure=FALSE, occu_model_type="occu"){
    sites_occ <- subset(df, !duplicated(site))$site
    # this (v v) function is synthetic species specific
    
    closed_df <- df
    if(!skip_closure){
        closed_df <- enforceClosure(df, occ_covs, sites_occ)
    } 
    
    
    umf_AUK <- auk::format_unmarked_occu(
        closed_df,
        site_id = "site",
        response = "species_observed",
        # response = "species_observed_syn",
        site_covs = occ_covs,
        obs_covs = det_covs
    )
    
    det_cov_str <- paste("", paste(det_covs, collapse="+"), sep=" ~ ")
    occ_cov_str <- paste("", paste(occ_covs, collapse="+"), sep=" ~ ")
    
    species_formula <- paste(det_cov_str, occ_cov_str, sep = " ")
    species_formula <- as.formula(species_formula)
    
    occ_um <- unmarked::formatWide(umf_AUK, type = "unmarkedFrameOccu")
    

    if (occu_model_type == "occu"){
   
        og_syn_gen_form <- unmarked::occu(formula = species_formula, data = occ_um)
    }
    
    return(og_syn_gen_form)
}
#######




########
predict_sdm_map <- function(occ_eq, region){
    
    valid_boundary <- terra::vect("occupancy_feature_raster/boundary/boundary.shp")
    crs(valid_boundary) <- crs(region)
    region <- terra::crop(region, valid_boundary, mask = TRUE)

    # print(paste0("SDM of occupancy probability for ", clustering_name))
    region$occ_prob <- expit(
                    occ_eq[[1]] + 
                    region$elevation * occ_eq[[2]] +
                    region$TCB * occ_eq[[3]] +
                    region$TCG * occ_eq[[4]] +
                    region$TCW * occ_eq[[5]] + 
                    region$TCA * occ_eq[[6]]
                    )
    
    return(region)
}
########


########
get_model_parameter_list <- function(df, i){

    occ.estimates_list <- list(
        occ_intercept = df[i,]$occ_intercept,
        elevation = df[i,]$elevation,
        TCB = df[i,]$TCB,
        TCG = df[i,]$TCG,
        TCW = df[i,]$TCW,
        TCA = df[i,]$TCA
    )

    det.estimates_list <- list(
        det_intercept = df[i,]$det_intercept,
        day_of_year = df[i,]$day_of_year,
        time_observations_started = df[i,]$time_observations_started,
        duration_minutes = df[i,]$duration_minutes,
        effort_distance_km = df[i,]$effort_distance_km,
        number_observers = df[i,]$number_observers
    )

    return (list(occ.estimates_list = occ.estimates_list, det.estimates_list = det.estimates_list))
}
########




