####################
# Main helper

# December 10, 2024
####################

library(auk)


# functionn to make test flow smoother
# test_names is a list of strings specifying
# the methods you want to run. the parameters
# of each test are separated with a '-'
#
# for example: clustGeo-.8-850, kmSq-1000, rounded-4
genExp <- function(method_names){
    
    methods <- list(
                                # base=NA,
                                ground_truth=NA, 

                                one_to_10=NA,
                                two_to_10=NA,
                                two_to_10_sameObs=NA,
                                kmSq=list(),
                                lat_long=NA,
                                rounded=list(), 
                                svs=NA, 
                                one_UL=NA,

                                clustGeo_25_60=NA,
                                clustGeo_50_60=NA,
                                clustGeo_75_60=NA,
                                clustGeo_25_70=NA,
                                clustGeo_50_70=NA,
                                clustGeo_75_70=NA, 
                                clustGeo_25_80=NA,
                                clustGeo_50_80=NA,
                                clustGeo_75_80=NA,
                                clustGeo_25_90=NA,
                                clustGeo_50_90=NA,
                                clustGeo_75_90=NA,
                                clustGeo_25_100=NA,
                                clustGeo_50_100=NA,
                                clustGeo_75_100=NA,

                                clustGeo_50_20=NA,
                                clustGeo_50_40=NA,
                                
                                DBSC=NA,
                                
                                BayesOptClustGeo=NA    

    )
    
    for(i in 1:length(method_names)){
        
        m_name <- strsplit(method_names[[i]], "-")[[1]]
        
        if(m_name[1] %in% c( "kmSq", "rounded")){
            
            len <- length(methods[[m_name[1]]]) + 1
            methods[[m_name[1]]][len] <- list(as.double(m_name[2:length(m_name)]))
        } else {
            methods[[m_name[[1]][1]]] <- T    
        }
    }
    
    return(methods)
}





# experiment to test similarity and mse of sp. clustering algs
# against the ground truth
# this runs the clustering aspect of most of the algorithms.
# evaluation occurs at a later step
getClusterings <- function(methods, og_data, occ_covs, det_covs, truth_df=data.frame()){
    
    
    results <- list()

    if(!is.na(methods$ground_truth)){
        set.seed(1)
        results[["ground-truth"]] <- truth_df 
    }
    ###############################################

    ###############################################
    # 1to10
    if(!is.na(methods$one_to_10)){
        set.seed(1)
        df_1to10 <- filter_repeat_visits(
            og_data,
            min_obs = 1,
            max_obs = 10,
            annual_closure = TRUE,
            date_var = "formatted_date",
            site_vars = c("locality_id")
        )
        results[["1to10"]] <- df_1to10
    }
    ###############################################
    

    ###############################################
    # 2to10
    if(!is.na(methods$two_to_10)){
        set.seed(1)
        df_2to10 <- filter_repeat_visits(
            og_data,
            min_obs = 2,
            max_obs = 10,
            annual_closure = TRUE,
            date_var = "formatted_date",
            site_vars = c("locality_id")
        )
        results[["2to10"]] <- df_2to10
    }
    ###############################################
    

    ###############################################
    # 2to10-sameObs
    if(!is.na(methods$two_to_10_sameObs)){
        set.seed(1)
        df_2to10_sameObs <- filter_repeat_visits(
            og_data,
            min_obs = 2,
            max_obs = 10,
            annual_closure = TRUE,
            date_var = "formatted_date",
            site_vars = c("locality_id", "observer_id")
        )
        results[["2to10-sameObs"]] <- df_2to10_sameObs
    }
    ###############################################



    ###############################################
    # kmSq
    if(length(methods$kmSq) > 0){
        set.seed(1)
        for(i in 1:length(methods$kmSq)){
            kmSq.df <- kmsq.Sites(og_data, rad_m = methods$kmSq[[i]])
            km_len <- round_any ((methods$kmSq[[i]] / 1000)^2, 0.125, ceiling)
            # results[[paste0("kmSq-", as.character(methods$kmSq[[i]]))]] <- kmSq.df
            results[[paste(km_len,"-kmSq",sep="")]] <- kmSq.df
        }
    }
    ###############################################

    
    ###############################################
    # lat-long
    if(!is.na(methods$lat_long)){
        set.seed(1)
        simple_df <- filter_repeat_visits(
            og_data,
            min_obs = 1,
            max_obs = 1000000,
            annual_closure = TRUE,
            date_var = "formatted_date",
            site_vars = c("locality_id")
        )
        results[["lat-long"]] <- simple_df
    }
    ###############################################


    
    ###############################################
    # rounded-4
    if(length(methods$rounded) > 0){
        set.seed(1)
        rounded_data <- roundLatLong(og_data, 4)
        rounded_df <- filter_repeat_visits(
            rounded_data,
            min_obs = 1,
            max_obs = 1000000,
            annual_closure = TRUE,
            date_var = "formatted_date",
            site_vars = c("rounded_locality_id")
        )
        results[[paste0("rounded-", as.character(4))]] <- rounded_df
    }
    ###############################################

    
    
    ###############################################
    # SVS
    if(!is.na(methods$svs)){
        set.seed(1)
        
        svs_df <- og_data

        # n_sites = dim(svs_df)[1]
        svs_df$site <- svs_df$checklist_id 

        results[["SVS"]] <- svs_df
    }
    ###############################################
    

    ###############################################
    # 1-UL
    if(!is.na(methods$one_UL)){
        set.seed(1)
     
        df_one_UL_t <- filter_repeat_visits(
            og_data,
            min_obs = 1,
            max_obs = 1000000,
            annual_closure = TRUE,
            date_var = "formatted_date",
            site_vars = c("locality_id")
        )
    
        one_UL_df <- # dplr, 
        df_one_UL_t %>%
            group_by(site) %>% #group rows by "site"
            filter(row_number()==1) #select first row from each group

        results[["1-UL"]] <- one_UL_df
    }
    ###############################################

    ###############################################
    # clustGeo-.25-60
    if(!is.na(methods$clustGeo_25_60)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.25
        percent <- 0.60
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-25-60"]] <- og_data_cgul
    
        
    }
    ###############################################

    ###############################################
    # clustGeo-.5-60
    if(!is.na(methods$clustGeo_50_60)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.5
        percent <- 0.60
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-50-60"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.75-60
    if(!is.na(methods$clustGeo_75_60)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.75
        percent <- 0.60
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-75-60"]] <- og_data_cgul    
    }
    ###############################################


    ###############################################
    # clustGeo-.25-70
    if(!is.na(methods$clustGeo_25_70)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.25
        percent <- 0.70
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-25-70"]] <- og_data_cgul
    
        
    }
    ###############################################

    ###############################################
    # clustGeo-.5-70
    if(!is.na(methods$clustGeo_50_70)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.5
        percent <- 0.70
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-50-70"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.75-70
    if(!is.na(methods$clustGeo_75_70)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.75
        percent <- 0.70
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-75-70"]] <- og_data_cgul    
    }
    ###############################################


    ###############################################
    # clustGeo-.25-80
    if(!is.na(methods$clustGeo_25_80)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.25
        percent <- 0.80
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-25-80"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.5-80
    if(!is.na(methods$clustGeo_50_80)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.5
        percent <- 0.80
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-50-80"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.75-80
    if(!is.na(methods$clustGeo_75_80)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.75
        percent <- 0.80
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-75-80"]] <- og_data_cgul    
    }
    ###############################################


    ###############################################
    # clustGeo-.25-90
    if(!is.na(methods$clustGeo_25_90)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.25
        percent <- 0.90
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-25-90"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.5-90
    if(!is.na(methods$clustGeo_50_90)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.5
        percent <- 0.90
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-50-90"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.75-90
    if(!is.na(methods$clustGeo_75_90)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.75
        percent <- 0.90
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-75-90"]] <- og_data_cgul    
    }
    ###############################################




    ###############################################
    # clustGeo-.25-100
    if(!is.na(methods$clustGeo_25_100)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.25
        percent <- 1.0
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-25-100"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.5-100
    if(!is.na(methods$clustGeo_50_100)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.5
        percent <- 1.0
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-50-100"]] <- og_data_cgul
    
        
    }
    ###############################################

        ###############################################
    # clustGeo-.75-100
    if(!is.na(methods$clustGeo_75_100)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.75
        percent <- 1.0
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-75-100"]] <- og_data_cgul    
    }
    ###############################################

    ###############################################
    # clustGeo-.50-20
    if(!is.na(methods$clustGeo_50_20)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.50
        percent <- 0.20
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-50-20"]] <- og_data_cgul    
    }
    ###############################################

    ###############################################
    # clustGeo-.50-40
    if(!is.na(methods$clustGeo_50_40)){
        set.seed(1)
        og_data_cgul <- og_data
        og_data_cgul$lat_long <- paste0(og_data$latitude, "_", og_data$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        alpha <- 0.50
        percent <- 0.40
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        results[["clustGeo-50-40"]] <- og_data_cgul    
    }
    ###############################################

    ###############################################
    # DBSC
    if(!is.na(methods$DBSC)){
        set.seed(1)
        DBSC_df <- runDBSC(og_data, occ_covs, det_covs)
    
        results[["DBSC"]] <- DBSC_df
    }
    ###############################################
    


    ###############################################
    # BayesOptClustGeo
    if(!is.na(methods$BayesOptClustGeo)){
        
        set.seed(1)

        og_data_cgul <- og_data

        bayesianOptimized <- bayesianOptimizedClustGeo(og_data_cgul, occ_covs, det_covs, "silhouette")
        alpha <- bayesianOptimized$Best_Pars$alpha
        percent <- bayesianOptimized$Best_Pars$lambda/100
        
        og_data_cgul$lat_long <- paste0(og_data_cgul$latitude, "_", og_data_cgul$longitude)
        uniq_loc_df <- dplyr::distinct(og_data_cgul, lat_long, .keep_all = T)
        
        clustGeo_df_i <- clustGeoSites(alpha = alpha, uniq_loc_df, occ_covs, det_covs, num_sites = round(nrow(uniq_loc_df)*percent))
        
        # link un labeled, but lat-long duplicated checklists to the correct site
        og_data_cgul$site <- -1
        for(j in seq(1:nrow(clustGeo_df_i))){
            og_data_cgul[og_data_cgul$lat_long == clustGeo_df_i[j,]$lat_long,]$site <- clustGeo_df_i[j,]$site
        }
        

        results[["BayesOptClustGeo"]] <- list(result_df = og_data_cgul, Best_Pars=bayesianOptimized$Best_Pars)
    
    }
 
    
    return(results)
}


parseMethodName <- function(method){
    if (method=="ground_truth"){
        return ("ground-truth")
    }
    else if (method=="one_to_10"){
        return ("1to10")
    }
    else if (method=="two_to_10"){
        return ("2to10")
    }
    else if (method=="two_to_10_sameObs"){
        return ("2to10-sameObs")
    }
    else if (method=="kmSq-1000"){
        return ("1-kmSq")
    }
    else if (method=="lat_long"){
        return ("lat-long")
    }
    else if (method=="svs"){
        return ("SVS")
    }
    else if (method=="one_UL"){
        return ("1-UL")
    }
    else if (startsWith(method, "clustGeo")){
        return (gsub("_", "-", method)) 
    }
    else {
        return (method)
    }
}