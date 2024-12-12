###########################################################
# File to run spatial clustering algorithms on eBird data

# December 10, 2024
###########################################################


library(terra) # raster operations
library(rje) # expit function
library(PRROC) # AUC and AUPRC calculation

det_covs <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")
occ_covs <- c("elevation", "TCB", "TCG", "TCW", "TCA")


OR.init <- rast("occupancy_feature_raster/occupancy_features.tif")
names(OR.init) <- occ_covs



f.name_train <- paste0(spec_name, "/", spec_name, "_zf_filtered_region_2017.csv")


train.df.og <- read.delim(paste0("checklist_data/species/", f.name_train), sep = ",", header = T)


train.df.og <- train.df.og[!is.na(train.df.og$duration_minutes),]


train.df.og <- train.df.og[train.df.og$observation_date >= "2017-05-15" & train.df.og$observation_date <= "2017-07-09",]





method_names <- as.list(method_names)

methods <- genExp(method_names)



test_split_path = "results/species_experiments/raw/test_data_splits/"



model_list <- list()


for(exp_run in seq(1:runs)){

    set.seed(1)

    train.df <- train.df.og
    test.df <- read.delim(paste0(test_split_path, "/", spec_name, "/run=", exp_run,".csv"), sep = ",", header = T)

    cat(paste0("\n\n\n", spec_name, " running experiment #", exp_run, "\n\n"))
    
    
    train_env.df <- extractEnvFeat(train.df, OR.init, obs_covs = det_covs)
    train.df <- inner_join(train.df, train_env.df, by = "checklist_id")
    norm_res <- norm_ds(train.df, det_covs, occ_covs)
    train.df <- norm_res$df
    norm.list <- norm_res$n_l
    
    train.df$occupied_prob <- -1
    train.df$det_prob <- -1
    
    
    test_env.df <- extractEnvFeat(test.df, OR.init, obs_covs = det_covs)
    test.df <- inner_join(test.df, test_env.df, by = "checklist_id")
    norm_res <- norm_ds(test.df, det_covs, occ_covs, test=TRUE, norm.list = norm.list)
    test.df <- norm_res$df
    
    test.df$occupied_prob <- -1
    test.df$det_prob <- -1

    train.df$formatted_date <- train.df$observation_date
    
    
    if(exp_run == 1){
        groupedSite <- getClusterings(methods, train.df, occ_covs, det_covs)
    }


    for(method_name in names(groupedSite)){
        
        set.seed(1)

        if (exp_run==1){

            
            if (startsWith(method_name, "BayesOpt")){

                best_par_df <- data.frame(parameter = names(groupedSite[[method_name]]$Best_Pars), value = as.numeric(groupedSite[[method_name]]$Best_Pars))
                
                write.csv(best_par_df, paste0("results/species_experiments/raw/clustering_parameters/", spec_name,"_", method_name,".csv"), row.names=FALSE)

                groupedSite[[method_name]] <- groupedSite[[method_name]]$result_df 

            }
            


            write.csv(groupedSite[[method_name]][,c("checklist_id", "longitude", "latitude", "site")], paste0("results/species_experiments/raw/clusterings/", spec_name,"_", method_name,".csv"), row.names=FALSE)

            cl_descr_stats <- calcDescriptiveClusteringStats(groupedSite[[method_name]])
            write.csv(data.frame(cl_descr_stats), paste0("results/species_experiments/raw/clusterings/", spec_name,"_", method_name,"_clust_descr_stats.csv"), row.names=FALSE)

            model_list[[method_name]] <- list()

            test.formula <- calcOccModel(groupedSite[[method_name]], occ_covs, det_covs)
       
            model_list[[method_name]][["occ.estimates_list"]] <- test.formula@estimates@estimates$state@estimates
            model_list[[method_name]][["det.estimates_list"]] <- test.formula@estimates@estimates$det@estimates

            occ.estimates_list <- model_list[[method_name]][["occ.estimates_list"]]
            det.estimates_list <- model_list[[method_name]][["det.estimates_list"]]

        
            model_pars_df <- data.frame(
                species = spec_name,
                method = method_name,
                occ_intercept = occ.estimates_list[[1]],
                elevation = occ.estimates_list[[2]],
                TCB = occ.estimates_list[[3]],
                TCG = occ.estimates_list[[4]],
                TCW = occ.estimates_list[[5]],
                TCA = occ.estimates_list[[6]],
                det_intercept = det.estimates_list[[1]],
                day_of_year = det.estimates_list[[2]],
                time_observations_started = det.estimates_list[[3]],
                duration_minutes = det.estimates_list[[4]],
                effort_distance_km = det.estimates_list[[5]],
                number_observers = det.estimates_list[[6]]
                
            )
            write.csv(model_pars_df, paste0("results/species_experiments/raw/model_parameters/", spec_name,"_", method_name, ".csv"), row.names=FALSE)
        }

 
        occ.estimates_list <- model_list[[method_name]][["occ.estimates_list"]]
        det.estimates_list <- model_list[[method_name]][["det.estimates_list"]]

        # check occupancy probability at each occupied site
        test.df$occupied_prob <- 
            occ.estimates_list[[1]] + 
            occ.estimates_list[[2]] * test.df[[names(occ.estimates_list)[[2]]]] +
            occ.estimates_list[[3]] * test.df[[names(occ.estimates_list)[[3]]]] +
            occ.estimates_list[[4]] * test.df[[names(occ.estimates_list)[[4]]]] +
            occ.estimates_list[[5]] * test.df[[names(occ.estimates_list)[[5]]]]    +
            occ.estimates_list[[6]] * test.df[[names(occ.estimates_list)[[6]]]]
        test.df$occupied_prob <- expit(test.df$occupied_prob)
        
        test.df$det_prob <-
            det.estimates_list[[1]] +
            det.estimates_list[[2]] * test.df[[names(det.estimates_list)[[2]]]] +
            det.estimates_list[[3]] * test.df[[names(det.estimates_list)[[3]]]] +
            det.estimates_list[[4]] * test.df[[names(det.estimates_list)[[4]]]] +
            det.estimates_list[[5]] * test.df[[names(det.estimates_list)[[5]]]] +
            det.estimates_list[[6]] * test.df[[names(det.estimates_list)[[6]]]]
        test.df$det_prob <- expit(test.df$det_prob)
        
        pred_observ <- unlist(test.df$occupied_prob*test.df$det_prob)
        
        dets <- pred_observ[test.df$species_observed == TRUE]
        nondets <- pred_observ[test.df$species_observed == FALSE]


        # AUC and AUPRC calculation

        # ROC Curve    
        roc <- roc.curve(scores.class0 = dets, scores.class1 = nondets, curve = F)

        # PR Curve
        pr <- pr.curve(scores.class0 = dets, scores.class1 = nondets, curve = F)


        predictions_df <- data.frame(
            checklist_id = test.df$checklist_id,
            longitude = test.df$longitude,
            latitude = test.df$latitude,
            species_observed = test.df$species_observed,
            prediction = pred_observ
        )
        write.csv(predictions_df, paste0("results/species_experiments/raw/predictions/", spec_name,"_", method_name, "_run=", exp_run,".csv"), row.names=FALSE)

        metrics_df <- data.frame(
            species = spec_name,
            run = exp_run,
            method = method_name,
            auc = roc$auc,
            auprc = pr$auc.integral
        )
            
        write.csv(metrics_df, paste0("results/species_experiments/raw/metrics/", spec_name,"_", method_name,"_run=", exp_run,".csv"), row.names=FALSE)


        cat(paste0(method_name, " has AUC: ", round(roc$auc, 4), " and AUPRC: ", round(pr$auc.integral, 4),"\n")) 
    }

}