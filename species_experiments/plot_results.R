#########################
# Plot results and maps

# December 10, 2024
#########################

library(dplyr)
library(ggplot2)
library(ggpubr) # ggarrange function
library(ggforce) # drawing ellipses over clusters
library(gridExtra) # plot grids
library(patchwork) # arranging plots of multiple sizes on same page
library(stringr) # string formatting
library(lmerTest) # linear mixed effect models
library(terra) # raster operations
library(FSA) # dunn's test

det_covs <- c("day_of_year", "time_observations_started", "duration_minutes", "effort_distance_km", "number_observers")
occ_covs <- c("elevation", "TCB", "TCG", "TCW", "TCA")

OR.normalized <- terra::rast("occupancy_feature_raster/occupancy_features_normalized.tif")
names(OR.normalized) <- occ_covs


placeholder_spec_name <- "AMCR"

cg_benchmark_method <- "clustGeo-25-80"
benchmark_method <- "lat-long"


results_clust_sum <- read.delim("results/species_experiments/summarized/results_clust_summarized.csv", sep = ",", header = T)
results_clust_sum <- subset(results_clust_sum[results_clust_sum$species==placeholder_spec_name,], select = -c(species))
write.csv(results_clust_sum, "results/species_experiments/summarized/results_clust_formatted.csv", row.names=FALSE)


results_sum <- read.delim("results/species_experiments/summarized/results_summarized.csv", sep = ",", header = T)
sp_descr <- read.delim("species_descr_ext.csv", sep = ",", header = T)



sp_descr$Prevalence.Level.Ext <- ""
sp_descr[sp_descr$Prevalence.Level == "l",]$Prevalence.Level.Ext <- "Low Prevalence"
sp_descr[sp_descr$Prevalence.Level == "m",]$Prevalence.Level.Ext <- "Medium Prevalence"
sp_descr[sp_descr$Prevalence.Level == "h",]$Prevalence.Level.Ext <- "High Prevalence"
sp_descr$Prevalence.Level.Ext = factor(sp_descr$Prevalence.Level.Ext, levels=c('Low Prevalence','Medium Prevalence','High Prevalence'))

sp_descr$Habitat.Ext <- ""
sp_descr[sp_descr$Habitat == "f",]$Habitat.Ext <- "Forested"
sp_descr[sp_descr$Habitat == "e",]$Habitat.Ext <- "Early Seral and Grasslands"
sp_descr$Habitat.Ext = factor(sp_descr$Habitat.Ext, levels=c('Forested','Early Seral and Grasslands'))

sp_descr$Generalist.Specialist.Ext <- ""
sp_descr[sp_descr$Generalist.Specialist == "g",]$Generalist.Specialist.Ext <- "Generalist"
sp_descr[sp_descr$Generalist.Specialist == "s",]$Generalist.Specialist.Ext <- "Specialist"
sp_descr$Generalist.Specialist.Ext = factor(sp_descr$Generalist.Specialist.Ext, levels=c('Generalist','Specialist'))


sp_descr$Home.Range.Ext <- ""
sp_descr[sp_descr$Home.Range == "s",]$Home.Range.Ext <- "Small Home Range"
sp_descr[sp_descr$Home.Range == "m",]$Home.Range.Ext <- "Medium Home Range"
sp_descr[sp_descr$Home.Range == "l",]$Home.Range.Ext <- "Large Home Range"
sp_descr$Home.Range.Ext = factor(sp_descr$Home.Range.Ext, levels=c('Small Home Range','Medium Home Range', 'Large Home Range'))


results_sum <- merge(results_sum, sp_descr, by.x = "species", by.y = "Abbreviation")


results_sum$Species <- as.factor(results_sum$Species)
results_sum$Prevalence.Level <- as.factor(results_sum$Prevalence.Level)



colors <- c( "forestgreen",  "darkgrey", "red",  "blue", "yellow", "orange", "purple", "green", "pink",  "cyan" )
colors <- rev(colors)

colors_cg <- c("darkcyan", "red", "blue", "yellow", "orange", "purple", "green", "pink",  "cyan", "brown", "forestgreen",  "darkgrey")

species_all <- unique(results_sum$species)



results_sum <- results_sum[(results_sum$method != "clustGeo-25-100" & results_sum$method != "clustGeo-50-100" & results_sum$method != "clustGeo-75-100"),]



results_sum_cg <- results_sum[startsWith(results_sum$method, "clustGeo"), ]
results_sum_cg$method = factor(results_sum_cg$method, levels=c("clustGeo-25-60","clustGeo-50-60","clustGeo-75-60","clustGeo-25-70","clustGeo-50-70","clustGeo-75-70","clustGeo-25-80","clustGeo-50-80","clustGeo-75-80","clustGeo-25-90","clustGeo-50-90","clustGeo-75-90"))

alg_all <- unique(results_sum$method)

alg_all_o <- c("2to10-sameObs", "2to10" , "rounded-4", "1-UL", "SVS",  "lat-long",  "1-kmSq",  "DBSC", "clustGeo", "BayesOptClustGeo")


##############################################################
##############################################################
# AUC results
##############################################################
##############################################################


# Calculate the mean AUC for each species
species_mean_auc <- results_sum_cg %>%
  group_by(Species) %>%
  dplyr::summarize(mean_auc = mean(auc, na.rm = TRUE)) %>%
  arrange(desc(mean_auc))


results_sum_cg_t <- results_sum_cg



# Reorder the species factor levels based on mean AUC
results_sum_cg_t$Species <- factor(results_sum_cg_t$Species, levels = rev(species_mean_auc$Species))


# clustGeo variants
png("results/species_experiments/plotted/plots/clustGeo_auc.png", units="in", width=10, height=16, res=400)


plot = ggplot(results_sum_cg_t, aes(x = Species, y = auc, fill=method)) +
    theme_classic() +
    geom_boxplot(outlier.alpha=0.25, outlier.size = 0.5, lwd = 0.35, width=0.9) +
    scale_fill_manual(values=rev(colors_cg)) + 
    labs(y= "AUC", x="Species", fill = "Algorithm") +
    coord_flip() + 
    theme(legend.position='right', legend.text=element_text(size=16), legend.title = element_text(size=20)) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=22)) 
   

show(plot)
dev.off()





calculateDiff <- function(df, benchmark){


    species_auc_diff_df <- data.frame(
        species = character(0),
        run = numeric(0),
        method = character(0),
        auc_diff = numeric(0),
        auc_perc_diff = numeric(0),
        auc_diff_from_species_mean = numeric(0),
        auc_perc_diff_from_species_mean = numeric(0)
    )

    auc_diff_df <- data.frame(
        run = numeric(0),
        method = character(0),
        auc_diff = numeric(0),
        auc_perc_diff = numeric(0),
        auc_diff_from_species_mean = numeric(0),
        auc_perc_diff_from_species_mean = numeric(0)
    )

    species_all <- unique(df$species)
    alg_all <- unique(df$method)
    runs <- as.numeric(as.character(unique(df$run)))


    for (r in runs){

        t_s <- data.frame(
            method = alg_all, 
            auc_diff = rep(0,length(alg_all)),
            auc_perc_diff = rep(0,length(alg_all)),
            auc_diff_from_species_mean = rep(0,length(alg_all)),
            auc_perc_diff_from_species_mean = rep(0,length(alg_all))
            
            
        )


        for (s in species_all){
            
            t_df <- df[df$run==r & df$species==s, ]

            auc_sp_mean <- mean(t_df$auc)

            species_benchmark <- t_df[t_df$method == benchmark,]$auc
            
            for (m in alg_all){
            
                cur_diff <- t_df[t_df$method==m, "auc"] - species_benchmark 
                cur_perc_diff <- (cur_diff/species_benchmark)*100

                cur_diff_from_species_mean <- t_df[t_df$method==m, "auc"] - auc_sp_mean 
                cur_perc_diff_from_species_mean <- (cur_diff_from_species_mean/auc_sp_mean)*100

                t_s[t_s$method==m, "auc_diff"] = t_s[t_s$method==m, "auc_diff"] + cur_diff
                t_s[t_s$method==m, "auc_perc_diff"] = t_s[t_s$method==m, "auc_perc_diff"] + cur_perc_diff

                t_s[t_s$method==m, "auc_diff_from_species_mean"] = t_s[t_s$method==m, "auc_diff_from_species_mean"] + cur_diff_from_species_mean
                t_s[t_s$method==m, "auc_perc_diff_from_species_mean"] = t_s[t_s$method==m, "auc_perc_diff_from_species_mean"] + cur_perc_diff_from_species_mean

                r_to_add <- data.frame(
                    species=s, 
                    run=r, 
                    method=m, 
                    auc_diff=cur_diff, 
                    auc_perc_diff=cur_perc_diff, 
                    auc_diff_from_species_mean=cur_diff_from_species_mean, 
                    auc_perc_diff_from_species_mean=cur_perc_diff_from_species_mean
                    )
                species_auc_diff_df <- rbind(species_auc_diff_df, r_to_add)

            
            }


        }
        # Average by number of species
        t_s[, c("auc_diff","auc_perc_diff","auc_diff_from_species_mean","auc_perc_diff_from_species_mean")] <- t_s[, c("auc_diff","auc_perc_diff","auc_diff_from_species_mean","auc_perc_diff_from_species_mean")]/length(species_all)
    
        for (m in alg_all){
            r_to_add <- cbind(r, m, t_s[t_s$method==m, c("auc_diff","auc_perc_diff","auc_diff_from_species_mean","auc_perc_diff_from_species_mean")])
            auc_diff_df <- rbind(auc_diff_df, r_to_add)
        }


    }

    colnames(species_auc_diff_df) <- c("species", "run", "method", "auc_diff", "auc_perc_diff", "auc_diff_from_species_mean", "auc_perc_diff_from_species_mean")
    colnames(auc_diff_df) <- c("run", "method", "auc_diff", "auc_perc_diff", "auc_diff_from_species_mean", "auc_perc_diff_from_species_mean")

    return (list(auc_diff_df = auc_diff_df, species_auc_diff_df = species_auc_diff_df))
}


tuneDf <- function (df, width = 1.0){

	sp_df <- data.frame(species = character(0), method = character(0), auc_mean = numeric(0), auc_sd = numeric(0))

	species_all <- unique(df$species)
    alg_all <- unique(df$method)


	for (s in species_all){

        auc_l <- data.frame(method=character(0), auc_mean=numeric(0), auc_sd=numeric(0))
        for (m in alg_all){
            auc_l <- rbind(auc_l, cbind(m, round(mean(df[df$species==s & df$method==m, ]$auc), 4), round(sd(df[df$species==s & df$method==m, ]$auc), 4)))
        }
        colnames(auc_l) <- c("method", "auc_mean", "auc_sd")

        
        w_variant <- auc_l[which.max(auc_l$auc_mean),]$method
        
  

        sp_df <- rbind(sp_df, cbind(s, w_variant , auc_l[auc_l$method == w_variant,]$auc_mean, auc_l[auc_l$method == w_variant,]$auc_sd))    
    }

    colnames(sp_df) <- c("species", "method", "auc_mean", "auc_sd")


    return(sp_df)
}




setTunedClustGeo <- function(df, tuned_df){

    species_all <- unique(df$species)

    for (s in species_all){

        df[(df$species==s) & (df$method == tuned_df[tuned_df$species==s,]$method),]$method  <- "best-clustGeo"

    }
    df <- subset(df, !(startsWith(df$method, "clustGeo-")))
    df$method <- as.factor(df$method)
    df$method <- factor(df$method, levels=c("2to10", "2to10-sameObs", "1-kmSq", "lat-long", "rounded-4", "SVS", "1-UL",  "best-clustGeo", "DBSC", "BayesOptClustGeo"))

    return(df)

}


diff_dfs <- calculateDiff(results_sum_cg, cg_benchmark_method)
auc_diff_df_cg <- diff_dfs$auc_diff_df
auc_diff_df_cg$method = factor(auc_diff_df_cg$method, levels=c("clustGeo-25-60","clustGeo-50-60","clustGeo-75-60","clustGeo-25-70","clustGeo-50-70","clustGeo-75-70","clustGeo-25-80","clustGeo-50-80","clustGeo-75-80","clustGeo-25-90","clustGeo-50-90","clustGeo-75-90"))



# clustGeo AUC difference plot
png("results/species_experiments/plotted/plots/auc_diff_clustGeo.png", units="in", width=5, height=6, res=400)

plot = ggplot(auc_diff_df_cg, aes(x = reorder(method, auc_diff), y = auc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=rev(colors_cg)) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= paste0("AUC Improvement over ", cg_benchmark_method), x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


# clustGeo percentage AUC difference plot
png("results/species_experiments/plotted/plots/auc_perc_diff_clustGeo.png", units="in", width=5, height=6, res=400)

plot = ggplot(auc_diff_df_cg, aes(x = reorder(method, auc_perc_diff), y = auc_perc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=rev(colors_cg)) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= paste0("Percentage AUC Improvement over ", cg_benchmark_method), x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


# clustGeo percentage AUC difference against species mean plot
png("results/species_experiments/plotted/plots/auc_perc_diff_species_mean_clustGeo.png", units="in", width=5, height=6, res=400)

plot = ggplot(auc_diff_df_cg, aes(x = reorder(method, auc_perc_diff_from_species_mean), y = auc_perc_diff_from_species_mean, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=rev(colors_cg)) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= "Percentage Improvement over Species Mean AUC", x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()



sink("results/species_experiments/summarized/tests_auc_clustGeo.txt")
cat("\nEffects of clustGeo clustering parameters on percentage AUC improvement over ", cg_benchmark_method,"\n\n")


# Kruskal-wallis test
kruskal_test_cg_auc <- kruskal.test(auc_perc_diff ~ method, data = auc_diff_df_cg) 
print(kruskal_test_cg_auc)

# Dunn's test
dunn_test_cg_auc <- dunnTest(auc_perc_diff ~ method, data = auc_diff_df_cg, method = "bh")
print(dunn_test_cg_auc)

sink()



# Save tuned clustGeo pars
tuned_cg_df <- tuneDf(results_sum_cg)
write.csv(tuned_cg_df, "results/species_experiments/summarized/clustGeo_tuned_auc.csv", row.names=FALSE)

# calculate differences between best-clustGeo and BayesOptClustGeo clustering parameters
tuned_cg_df_formatted <- tuned_cg_df

tuned_cg_df_formatted[,c("base","best_alpha","best_lambda")] <- str_split_fixed(tuned_cg_df_formatted$method, "-", 3)

tuned_cg_df_formatted$best_alpha <- as.numeric(tuned_cg_df_formatted$best_alpha)/100
tuned_cg_df_formatted$best_alpha <- as.numeric(tuned_cg_df_formatted$best_alpha)

pars_results <- read.delim("results/species_experiments/summarized/best_pars_summarized.csv", sep = ",", header = T)

tuned_cg_df_formatted$bayesopt_alpha <- ""
tuned_cg_df_formatted$bayesopt_lambda <- ""

for (s in species_all){

    tuned_cg_df_formatted[tuned_cg_df_formatted$species==s,]$bayesopt_alpha <- pars_results[(pars_results$species==s) & (pars_results$par_name=="alpha"),]$best_par
    tuned_cg_df_formatted[tuned_cg_df_formatted$species==s,]$bayesopt_lambda <- pars_results[(pars_results$species==s) & (pars_results$par_name=="lambda"),]$best_par

}

tuned_cg_df_formatted$diff_alpha <- round(as.numeric(tuned_cg_df_formatted$best_alpha) - as.numeric(tuned_cg_df_formatted$bayesopt_alpha), 4)
tuned_cg_df_formatted$diff_lambda <- round(as.numeric(tuned_cg_df_formatted$best_lambda) - as.numeric(tuned_cg_df_formatted$bayesopt_lambda))


col_order <- c("species", "best_alpha", "best_lambda", "diff_alpha", "diff_lambda", "auc_mean", "auc_sd")
tuned_cg_df_formatted <- tuned_cg_df_formatted[, col_order]

write.csv(tuned_cg_df_formatted, "results/species_experiments/summarized/clustGeo_tuned_auc_formatted.csv", row.names=FALSE)




results_sum_t <- results_sum



results_sum_t <- setTunedClustGeo(results_sum_t, tuned_cg_df)


tuned_df <- tuneDf(results_sum_t)


write.csv(tuned_df, "results/species_experiments/summarized/tuned_auc.csv", row.names=FALSE)



# # All algorithms (tuned clustGeo)
results_sum_t_t <- results_sum_t

species_mean_auc <- results_sum_t_t %>%
  group_by(Species) %>%
  dplyr::summarize(mean_auc = mean(auc, na.rm = TRUE)) %>%
  arrange(desc(mean_auc))



# Reorder the species factor levels based on mean AUC
results_sum_t_t$Species <- factor(results_sum_t_t$Species, levels = rev(species_mean_auc$Species))


png("results/species_experiments/plotted/plots/auc.png", units="in", width=12, height=16, res=400)


plot = ggplot(results_sum_t_t, aes(x = Species, y = auc, fill=method)) +
    theme_classic() +
    geom_boxplot(outlier.alpha=0.25, outlier.size = 0.5, lwd = 0.35, width=0.9) +
    scale_fill_manual(values=colors) + 
    labs(y= "AUC", x="Species", fill = "Algorithm") +
    coord_flip() + 
    theme(legend.position='right', legend.text=element_text(size=16), legend.title = element_text(size=20)) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=22))



show(plot)
dev.off()


diff_dfs <- calculateDiff(results_sum_t, benchmark_method)
auc_diff_df <- diff_dfs$auc_diff_df

auc_diff_df$method = factor(auc_diff_df$method, levels=c("2to10", "2to10-sameObs", "1-kmSq", "lat-long", "rounded-4", "SVS", "1-UL", "best-clustGeo", "DBSC", "BayesOptClustGeo"))



# AUC difference plot (clustGeo = tuned)
png("results/species_experiments/plotted/plots/auc_diff.png", units="in", width=5, height=6, res=400)


plot = ggplot(auc_diff_df, aes(x = reorder(method, auc_diff), y = auc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=colors) +
    ylim(NA, 0.025) +
    theme(axis.text.x = element_text(angle = 55,vjust = 1, hjust=1)) +
    labs(y= paste0("AUC Improvement over ", benchmark_method), x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


# Percentage AUC difference plot (clustGeo = tuned)
png("results/species_experiments/plotted/plots/auc_perc_diff.png", units="in", width=5, height=6, res=400)

plot = ggplot(auc_diff_df, aes(x = reorder(method, auc_perc_diff), y = auc_perc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=colors) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= paste0("Percentage AUC Improvement over ", benchmark_method), x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()

# Percentage AUC difference compared to species mean plot (clustGeo = tuned) 
png("results/species_experiments/plotted/plots/auc_perc_diff_species_mean.png", units="in", width=5, height=6, res=400)

plot = ggplot(auc_diff_df, aes(x = reorder(method, auc_perc_diff_from_species_mean), y = auc_perc_diff_from_species_mean, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=colors) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= "Percentage Improvement over Species Mean AUC", x="Algorithm", fill = "Algorithm") +
    scale_y_continuous(limits = c(-8.0, 4.0), breaks = c(-8.0,-4.0,-6.0,-2.0,0,2.0,4.0))  +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


sink("results/species_experiments/summarized/tests_auc.txt")
cat("\nEffects of clustering algorithms on percentage AUC improvement over ", benchmark_method,"\n\n")

# Kruskal-wallis test
kruskal_test_auc <- kruskal.test(auc_perc_diff ~ method, data = auc_diff_df) 
print(kruskal_test_auc)

# Dunn's test
dunn_test_auc <- dunnTest(auc_perc_diff ~ method, data = auc_diff_df, method = "bh")
print(dunn_test_auc)
sink()




##############################################################
# # Mixed-effect models - AUC
##############################################################

auc_diff_df_n <- diff_dfs$species_auc_diff_df

auc_diff_df_n <- merge(auc_diff_df_n, results_sum_t, by.x = c("species","method","run"), by.y = c("species","method","run"))


auc_diff_df_n$method <- as.factor(auc_diff_df_n$method)
auc_diff_df_n$prevalence <- as.factor(auc_diff_df_n$prevalence)
auc_diff_df_n$habitat <- as.factor(auc_diff_df_n$habitat)
auc_diff_df_n$specialist <- as.factor(auc_diff_df_n$specialist)
auc_diff_df_n$home_range <- as.factor(auc_diff_df_n$home_range)


# print(unique(auc_diff_df_n$method))
auc_diff_df_n$method <- factor(auc_diff_df_n$method, levels=c( "2to10-sameObs" , "SVS", "1-UL",  "lat-long", "rounded-4", "1-kmSq",  "DBSC", "BayesOptClustGeo", "best-clustGeo", "2to10"))

auc_diff_df_n$prevalence <- factor(auc_diff_df_n$prevalence, levels = c("low", "medium", "high"))
auc_diff_df_n$habitat <- factor(auc_diff_df_n$habitat, levels = c("forested", "seral"))
auc_diff_df_n$specialist <- factor(auc_diff_df_n$specialist, levels = c("generalist", "specialist"))
auc_diff_df_n$home_range <- factor(auc_diff_df_n$home_range, levels = c("small", "medium", "large"))



# Set reference levels of categorical predictors

auc_diff_df_n$prevalence = relevel(auc_diff_df_n$prevalence, ref="low")
auc_diff_df_n$habitat = relevel(auc_diff_df_n$habitat, ref="forested")
auc_diff_df_n$specialist = relevel(auc_diff_df_n$specialist, ref="generalist")
auc_diff_df_n$home_range = relevel(auc_diff_df_n$home_range, ref="small")
auc_diff_df_n$method = relevel(auc_diff_df_n$method, ref="2to10-sameObs")


auc_diff_df_n$genspec = auc_diff_df_n$specialist



lmer_algortihm <- lmer(auc_perc_diff ~ method * + (1 | species), data = auc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm <- summary(lmer_algortihm)

lmer_algorithm_prevalence <- lmer(auc_perc_diff ~ method:prevalence + (1 | species), data = auc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_prevalence <- summary(lmer_algorithm_prevalence)

lmer_algorithm_habitat <- lmer(auc_perc_diff ~ method:habitat + (1 | species), data = auc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_habitat <- summary(lmer_algorithm_habitat)

lmer_algorithm_specialist <- lmer(auc_perc_diff ~ method:genspec + (1 | species), data = auc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_specialist <- summary(lmer_algorithm_specialist)

lmer_algorithm_hrange <- lmer(auc_perc_diff ~ method:home_range + (1 | species), data = auc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_hrange <- summary(lmer_algorithm_hrange)


# Custom function to format and print coefficient matrix with significance stars
get_coefmat <- function(coefmat, digits = 6) {
    # Round the coefficients matrix to the specified number of digits
    rounded_coefmat <- round(coefmat, digits = digits)
    
    # Add a column for significance stars
    p_values <- rounded_coefmat[, "Pr(>|t|)"]
    stars <- ifelse(p_values < 0.001, "***",
            ifelse(p_values < 0.01, "**",
            ifelse(p_values < 0.05, "*",
            ifelse(p_values < 0.1, ".", " "))))
    
    # Combine rounded coefficients and stars into a data frame
    df_coefmat <- as.data.frame(rounded_coefmat)
    df_coefmat$signif_code <- stars
    
    return(df_coefmat)
}

# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auc_algorithm_mixed_coefficients.csv", row.names=TRUE)




# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_prevalence$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auc_algorithm_prevalence_mixed_coefficients.csv", row.names=TRUE)


# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_habitat$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auc_algorithm_habitat_mixed_coefficients.csv", row.names=TRUE)


# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_specialist$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auc_algorithm_specialist_mixed_coefficients.csv", row.names=TRUE)


# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_hrange$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auc_algorithm_homerange_mixed_coefficients.csv", row.names=TRUE)

# Plot mixed effect model results


plotCoefficients <- function(df_path, output_path, dimensions, title_str, rem_str, title_fontsize) {
    
    
    # Load the data
    coef_data <- read.csv(df_path, row.names = 1)
    
    # Remove the intercept term
    coef_data <- coef_data[-1,]
    
    coef_data$Coefficient <- rownames(coef_data)
    
    # Remove the word "method" from coefficient names
    coef_data$Coefficient <- gsub("method", "", coef_data$Coefficient)
    
    # Remove rem_str from coefficient names
    coef_data$Coefficient <- gsub(rem_str, "", coef_data$Coefficient)
    
    # Rename 'Std. Error' to 'Std_Error' to avoid space issues
    if ("Std..Error" %in% colnames(coef_data)) {
        colnames(coef_data)[colnames(coef_data) == "Std..Error"] <- "Std_Error"
    }
    
    # Reorder the data for plotting
    coef_data <- coef_data[order(coef_data$Estimate), ]
    
    # Define the possible significance levels and their properties
    signif_levels <- c("***", "**", "*", ".", " ")
    signif_colors <- c("***" = "red", "**" = "orange", "*" = "green", "." = "blue", " " = "black")
    signif_labels <- c("***" = "p-value < 0.001", "**" = "p-value < 0.01", "*" = "p-value < 0.05", "." = "p-value < 0.1", " " = "p-value >= 0.1")
    
    # Filter the significance levels present in the data
    present_signif_levels <- intersect(signif_levels, unique(coef_data$signif_code))
    
    # Calculate the number of rows for the legend dynamically
    legend_rows <- ceiling(length(present_signif_levels) / 2)
    
    # Create the coefficient plot
    plot <- ggplot(coef_data, aes(x = Estimate, y = reorder(Coefficient, Estimate))) +
        geom_errorbarh(aes(xmin = Estimate - Std_Error, xmax = Estimate + Std_Error), color = "darkgrey", height = 0.15) +
        geom_point(aes(color = signif_code), size = 2) +
        labs(title = paste0("", title_str),
             x = "Coefficient",
             y = "Effect") +
        scale_color_manual(
            name = "Significance",
            values = signif_colors[present_signif_levels],
            labels = signif_labels[present_signif_levels]
        ) +
        theme_bw() +
        theme(axis.text.y = element_text(size = 10, angle = 0, hjust = 1),
              legend.position = "bottom",
              legend.justification = "center", # Center the legend
              legend.box = "horizontal",
              legend.title = element_text(size = 10, hjust = 0),  # Align legend title to the left
              legend.text = element_text(size = 9, hjust = 0),   # Align legend text to the left
              plot.title = element_text(size = title_fontsize, hjust = 0)) +  # Move the title further to the left
        guides(color = guide_legend(nrow = legend_rows, byrow = TRUE))
    
    # Save the plot to a 300 dpi PNG file with specified dimensions
    ggsave(output_path, plot = plot, dpi = 300, width = dimensions[1], height = dimensions[2])
}




plotCoefficients(
    "results/species_experiments/summarized/auc_algorithm_mixed_coefficients.csv", 
    "results/species_experiments/plotted/plots/mixed_coefficients_algorithms_auc.png", 
    c(6,5),
    "% AUC Improvement ~ Algorithm",
    "",
    9
)


plotCoefficientsTraits <- function(df_paths, output_path, dimensions, title_strs, rem_strs, title_fontsize) {
    
    plots <- list()
    combined_signif_levels <- c()
    
    for (i in 1:4) {
        # Load the data
        coef_data <- read.csv(df_paths[i], row.names = 1)
        
        # Remove the intercept term
        coef_data <- coef_data[-1,]
        
        coef_data$Coefficient <- rownames(coef_data)
        
        # Remove the word "method" from coefficient names
        coef_data$Coefficient <- gsub("method", "", coef_data$Coefficient)
        
        # Remove rem_str from coefficient names
        coef_data$Coefficient <- gsub(rem_strs[i], "", coef_data$Coefficient)
        
        # Rename 'Std. Error' to 'Std_Error' to avoid space issues
        if ("Std..Error" %in% colnames(coef_data)) {
            colnames(coef_data)[colnames(coef_data) == "Std..Error"] <- "Std_Error"
        }
        
        # Reorder the data for plotting
        coef_data <- coef_data[order(coef_data$Estimate), ]
        
        # Define the possible significance levels and their properties
        signif_levels <- c("***", "**", "*", ".", " ")
        signif_colors <- c("***" = "red", "**" = "orange", "*" = "green", "." = "blue", " " = "black")
        signif_labels <- c("***" = "p-value < 0.001", "**" = "p-value < 0.01", "*" = "p-value < 0.05", "." = "p-value < 0.1", " " = "p-value >= 0.1")
        
        # Filter the significance levels present in the data
        present_signif_levels <- intersect(signif_levels, unique(coef_data$signif_code))
        combined_signif_levels <- unique(c(combined_signif_levels, present_signif_levels))
        
        # Create the coefficient plot
        p <- ggplot(coef_data, aes(x = Estimate, y = reorder(Coefficient, Estimate))) +
            geom_errorbarh(aes(xmin = Estimate - Std_Error, xmax = Estimate + Std_Error), color = "darkgrey", height = 0.15) +
            geom_point(aes(color = signif_code), size = 2) +
            labs(title = title_strs[i],
                 x = "Coefficient",
                 y = "Effect") +
            scale_color_manual(
                name = "Significance",
                values = signif_colors[signif_levels],
                labels = signif_labels[signif_levels]
            ) +
            theme_bw() +
            theme(axis.text.y = element_text(size = 15, angle = 0, hjust = 1),
                  axis.title.y = element_text(size = 15),
                  axis.text.x = element_text(size = 14),
                  axis.title.x = element_text(size = 14), 
                  plot.title = element_text(size = title_fontsize-2, hjust = 0),
                  legend.position = "bottom",
                  legend.box = "horizontal",
                  legend.title = element_text(size = 18, hjust = 0),
                  legend.text = element_text(size = 16, hjust = 0),
                #   legend.spacing.x = unit(0.5, 'in')
                ) +
            guides(color = guide_legend(nrow = 1, byrow = TRUE))
        
        plots[[i]] <- p
    }
    
    # Extract the combined legend
    combined_legend <- ggplot(data.frame(signif_code = combined_signif_levels), aes(x = 1, y = 1, color = signif_code)) +
        geom_point() +
        scale_color_manual(
            name = "Significance",
            values = signif_colors[combined_signif_levels],
            labels = signif_labels[combined_signif_levels]
        ) +
        theme_void() +
        theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.title = element_text(size = 18, hjust = 0, margin = margin(r = 0.35, unit='in')),
            legend.text = element_text(size = 16, hjust = 0, margin = margin(t = 0, unit='in')),
            legend.key.spacing.x = unit(0.35, 'in'),
            legend.key.width = unit(0.5, 'in')) +
        guides(color = guide_legend(nrow = 1, byrow = TRUE))
    
    # Extract the legend from the combined legend plot
    g_legend <- function(a.gplot) {
        tmp <- ggplotGrob(a.gplot)
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
    }
    
    legend <- g_legend(combined_legend)
    
    # Remove legends from individual plots
    plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
    
    # Combine the plots into a 2x2 grid
    combined_plot <- arrangeGrob(
        grobs = plots,
        ncol = 2,
        nrow = 2,
        # top = textGrob("Interactions of Species Traits, Clustering Algorithm, and AUC\n", gp = gpar(fontsize = title_fontsize)),
        bottom = legend
    )
    
    # Save the combined plot to a 300 dpi PNG file with specified dimensions
    ggsave(output_path, plot = combined_plot, dpi = 300, width = dimensions[1], height = dimensions[2])
}


plotCoefficientsTraits(
    c(
        "results/species_experiments/summarized/auc_algorithm_prevalence_mixed_coefficients.csv", 
        "results/species_experiments/summarized/auc_algorithm_homerange_mixed_coefficients.csv",
        "results/species_experiments/summarized/auc_algorithm_habitat_mixed_coefficients.csv", 
        "results/species_experiments/summarized/auc_algorithm_specialist_mixed_coefficients.csv"
    ),
    "results/species_experiments/plotted/plots/mixed_coefficients_traits_auc.png",
    c(15, 14),
    c("Prevalence",
      "Home Range",
      "Habitat",
      "Generalist/Specialist"),
    c("prevalence", "home_range", "habitat", "genspec"),
    20
)

##############################################################



##############################################################
##############################################################
# AUPRC results
##############################################################
##############################################################

# Calculate the mean AUPRC for each species
species_mean_auprc <- results_sum_cg %>%
  group_by(Species) %>%
  dplyr::summarize(mean_auprc = mean(auprc, na.rm = TRUE)) %>%
  arrange(desc(mean_auprc))


results_sum_cg_t <- results_sum_cg



# Reorder the species factor levels based on mean AUPRC
results_sum_cg_t$Species <- factor(results_sum_cg_t$Species, levels = rev(species_mean_auprc$Species))


# clustGeo variants
png("results/species_experiments/plotted/plots/clustGeo_auprc.png", units="in", width=10, height=16, res=400)


plot = ggplot(results_sum_cg_t, aes(x = Species, y = auprc, fill=method)) +
    theme_classic() +
    geom_boxplot(outlier.alpha=0.25, outlier.size = 0.5, lwd = 0.35, width=0.9) +
    scale_fill_manual(values=rev(colors_cg)) + 
    labs(y= "AUPRC", x="Species", fill = "Algorithm") +
    coord_flip() + 
    theme(legend.position='right', legend.text=element_text(size=16), legend.title = element_text(size=20)) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=22)) 

show(plot)
dev.off()



calculateDiff <- function(df, benchmark){


    species_auprc_diff_df <- data.frame(
        species = character(0),
        run = numeric(0),
        method = character(0),
        auprc_diff = numeric(0),
        auprc_perc_diff = numeric(0),
        auprc_diff_from_species_mean = numeric(0),
        auprc_perc_diff_from_species_mean = numeric(0)
    )

    auprc_diff_df <- data.frame(
        run = numeric(0),
        method = character(0),
        auprc_diff = numeric(0),
        auprc_perc_diff = numeric(0),
        auprc_diff_from_species_mean = numeric(0),
        auprc_perc_diff_from_species_mean = numeric(0)
    )

    species_all <- unique(df$species)
    alg_all <- unique(df$method)
    runs <- as.numeric(as.character(unique(df$run)))


    for (r in runs){

        t_s <- data.frame(
            method = alg_all, 
            auprc_diff = rep(0,length(alg_all)),
            auprc_perc_diff = rep(0,length(alg_all)),
            auprc_diff_from_species_mean = rep(0,length(alg_all)),
            auprc_perc_diff_from_species_mean = rep(0,length(alg_all))
            
            
        )


        for (s in species_all){
            
            t_df <- df[df$run==r & df$species==s, ]

            auprc_sp_mean <- mean(t_df$auprc)

            species_benchmark <- t_df[t_df$method == benchmark,]$auprc
            
            for (m in alg_all){
            
                cur_diff <- t_df[t_df$method==m, "auprc"] - species_benchmark 
                cur_perc_diff <- (cur_diff/species_benchmark)*100

                cur_diff_from_species_mean <- t_df[t_df$method==m, "auprc"] - auprc_sp_mean 
                cur_perc_diff_from_species_mean <- (cur_diff_from_species_mean/auprc_sp_mean)*100

                t_s[t_s$method==m, "auprc_diff"] = t_s[t_s$method==m, "auprc_diff"] + cur_diff
                t_s[t_s$method==m, "auprc_perc_diff"] = t_s[t_s$method==m, "auprc_perc_diff"] + cur_perc_diff

                t_s[t_s$method==m, "auprc_diff_from_species_mean"] = t_s[t_s$method==m, "auprc_diff_from_species_mean"] + cur_diff_from_species_mean
                t_s[t_s$method==m, "auprc_perc_diff_from_species_mean"] = t_s[t_s$method==m, "auprc_perc_diff_from_species_mean"] + cur_perc_diff_from_species_mean

                r_to_add <- data.frame(
                    species=s, 
                    run=r, 
                    method=m, 
                    auprc_diff=cur_diff, 
                    auprc_perc_diff=cur_perc_diff, 
                    auprc_diff_from_species_mean=cur_diff_from_species_mean, 
                    auprc_perc_diff_from_species_mean=cur_perc_diff_from_species_mean
                    )
                species_auprc_diff_df <- rbind(species_auprc_diff_df, r_to_add)

            
            }


        }
        # Average by number of species
        t_s[, c("auprc_diff","auprc_perc_diff","auprc_diff_from_species_mean","auprc_perc_diff_from_species_mean")] <- t_s[, c("auprc_diff","auprc_perc_diff","auprc_diff_from_species_mean","auprc_perc_diff_from_species_mean")]/length(species_all)
    
        for (m in alg_all){
            r_to_add <- cbind(r, m, t_s[t_s$method==m, c("auprc_diff","auprc_perc_diff","auprc_diff_from_species_mean","auprc_perc_diff_from_species_mean")])
            auprc_diff_df <- rbind(auprc_diff_df, r_to_add)
        }


    }

    colnames(species_auprc_diff_df) <- c("species", "run", "method", "auprc_diff", "auprc_perc_diff", "auprc_diff_from_species_mean", "auprc_perc_diff_from_species_mean")
    colnames(auprc_diff_df) <- c("run", "method", "auprc_diff", "auprc_perc_diff", "auprc_diff_from_species_mean", "auprc_perc_diff_from_species_mean")

    return (list(auprc_diff_df = auprc_diff_df, species_auprc_diff_df = species_auprc_diff_df))
}



tuneDf <- function (df, width = 1.0){

	sp_df <- data.frame(species = character(0), method = character(0), auprc_mean = numeric(0), auprc_sd = numeric(0))

	species_all <- unique(df$species)
    alg_all <- unique(df$method)


	for (s in species_all){

        auprc_l <- data.frame(method=character(0), auprc_mean=numeric(0), auprc_sd=numeric(0))
        for (m in alg_all){
            auprc_l <- rbind(auprc_l, cbind(m, round(mean(df[df$species==s & df$method==m, ]$auprc), 4), round(sd(df[df$species==s & df$method==m, ]$auprc), 4)))
        }
        colnames(auprc_l) <- c("method", "auprc_mean", "auprc_sd")

        
        w_variant <- auprc_l[which.max(auprc_l$auprc_mean),]$method
        
  

        sp_df <- rbind(sp_df, cbind(s, w_variant , auprc_l[auprc_l$method == w_variant,]$auprc_mean, auprc_l[auprc_l$method == w_variant,]$auprc_sd))    
    }

    # print(dim(sp_df))
    colnames(sp_df) <- c("species", "method", "auprc_mean", "auprc_sd")
    # sp_df <- "auprc_mean_difference_from_tied"

    return(sp_df)
}




setTunedClustGeo <- function(df, tuned_df){

    species_all <- unique(df$species)

    for (s in species_all){

        df[(df$species==s) & (df$method == tuned_df[tuned_df$species==s,]$method),]$method  <- "best-clustGeo"

    }
    df <- subset(df, !(startsWith(df$method, "clustGeo-")))
    df$method <- as.factor(df$method)
    df$method <- factor(df$method, levels=c("2to10", "2to10-sameObs", "1-kmSq", "lat-long", "rounded-4", "SVS", "1-UL", "best-clustGeo", "DBSC", "BayesOptClustGeo"))

    return(df)

}


diff_dfs <- calculateDiff(results_sum_cg, cg_benchmark_method)
auprc_diff_df_cg <- diff_dfs$auprc_diff_df
auprc_diff_df_cg$method = factor(auprc_diff_df_cg$method, levels=c("clustGeo-25-60","clustGeo-50-60","clustGeo-75-60","clustGeo-25-70","clustGeo-50-70","clustGeo-75-70","clustGeo-25-80","clustGeo-50-80","clustGeo-75-80","clustGeo-25-90","clustGeo-50-90","clustGeo-75-90"))



# clustGeo AUPRC difference plot
png("results/species_experiments/plotted/plots/auprc_diff_clustGeo.png", units="in", width=5, height=6, res=400)

plot = ggplot(auprc_diff_df_cg, aes(x = reorder(method, auprc_diff), y = auprc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=rev(colors_cg)) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= paste0("AUPRC Improvement over ", cg_benchmark_method), x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


# clustGeo percentage AUPRC difference plot
png("results/species_experiments/plotted/plots/auprc_perc_diff_clustGeo.png", units="in", width=5, height=6, res=400)

plot = ggplot(auprc_diff_df_cg, aes(x = reorder(method, auprc_perc_diff), y = auprc_perc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=rev(colors_cg)) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= paste0("Percentage AUPRC Improvement over ", cg_benchmark_method), x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


# clustGeo percentage AUPRC difference against speciec mean plot
png("results/species_experiments/plotted/plots/auprc_perc_diff_species_mean_clustGeo.png", units="in", width=5, height=6, res=400)

plot = ggplot(auprc_diff_df_cg, aes(x = reorder(method, auprc_perc_diff_from_species_mean), y = auprc_perc_diff_from_species_mean, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=rev(colors_cg)) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= "Percentage Improvement over Species Mean AUPRC", x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


sink("results/species_experiments/summarized/tests_auprc_clustGeo.txt")
cat("\nEffects of clustGeo clustering parameters on percentage AUPRC improvement over ", cg_benchmark_method,"\n\n")

# Kruskal-wallis test
kruskal_test_cg_auprc <- kruskal.test(auprc_perc_diff ~ method, data = auprc_diff_df_cg) 
print(kruskal_test_cg_auprc)

# Dunn's test
dunn_test_cg_auprc <- dunnTest(auprc_perc_diff ~ method, data = auprc_diff_df_cg, method = "bh")
print(dunn_test_cg_auprc)

sink()


# Set tuned clustGeo pars
tuned_cg_df <- tuneDf(results_sum_cg)
write.csv(tuned_cg_df, "results/species_experiments/summarized/clustGeo_tuned_auprc.csv", row.names=FALSE)

# calculate differences between best-clustGeo and BayesOptClustGeo clustering parameters
tuned_cg_df_formatted <- tuned_cg_df

tuned_cg_df_formatted[,c("base","best_alpha","best_lambda")] <- str_split_fixed(tuned_cg_df_formatted$method, "-", 3)

tuned_cg_df_formatted$best_alpha <- as.numeric(tuned_cg_df_formatted$best_alpha)/100
tuned_cg_df_formatted$best_alpha <- as.numeric(tuned_cg_df_formatted$best_alpha)

pars_results <- read.delim("results/species_experiments/summarized/best_pars_summarized.csv", sep = ",", header = T)

tuned_cg_df_formatted$bayesopt_alpha <- ""
tuned_cg_df_formatted$bayesopt_lambda <- ""

for (s in species_all){

    tuned_cg_df_formatted[tuned_cg_df_formatted$species==s,]$bayesopt_alpha <- pars_results[(pars_results$species==s) & (pars_results$par_name=="alpha"),]$best_par
    tuned_cg_df_formatted[tuned_cg_df_formatted$species==s,]$bayesopt_lambda <- pars_results[(pars_results$species==s) & (pars_results$par_name=="lambda"),]$best_par

}

tuned_cg_df_formatted$diff_alpha <- round(as.numeric(tuned_cg_df_formatted$best_alpha) - as.numeric(tuned_cg_df_formatted$bayesopt_alpha), 4)
tuned_cg_df_formatted$diff_lambda <- round(as.numeric(tuned_cg_df_formatted$best_lambda) - as.numeric(tuned_cg_df_formatted$bayesopt_lambda))


col_order <- c("species", "best_alpha", "best_lambda", "diff_alpha", "diff_lambda", "auprc_mean", "auprc_sd")
tuned_cg_df_formatted <- tuned_cg_df_formatted[, col_order]

write.csv(tuned_cg_df_formatted, "results/species_experiments/summarized/clustGeo_tuned_auprc_formatted.csv", row.names=FALSE)




results_sum_t <- results_sum



results_sum_t <- setTunedClustGeo(results_sum_t, tuned_cg_df)


tuned_df <- tuneDf(results_sum_t)


write.csv(tuned_df, "results/species_experiments/summarized/tuned_auprc.csv", row.names=FALSE)



# # All algorithms (tuned clustGeo)
results_sum_t_t <- results_sum_t

species_mean_auprc <- results_sum_t_t %>%
  group_by(Species) %>%
  dplyr::summarize(mean_auprc = mean(auprc, na.rm = TRUE)) %>%
  arrange(desc(mean_auprc))



# Reorder the species factor levels based on mean AUPRC
results_sum_t_t$Species <- factor(results_sum_t_t$Species, levels = rev(species_mean_auprc$Species))


png("results/species_experiments/plotted/plots/auprc.png", units="in", width=12, height=16, res=400)


plot = ggplot(results_sum_t_t, aes(x = Species, y = auprc, fill=method)) +
    theme_classic() +
    geom_boxplot(outlier.alpha=0.25, outlier.size = 0.5, lwd = 0.35, width=0.9) +
    scale_fill_manual(values=colors) + 
    labs(y= "AUPRC", x="Species", fill = "Algorithm") +
    coord_flip() + 
    theme(legend.position='right', legend.text=element_text(size=16), legend.title = element_text(size=20)) +
    theme(axis.text=element_text(size=16), axis.title=element_text(size=22))



show(plot)
dev.off()


diff_dfs <- calculateDiff(results_sum_t, benchmark_method)
auprc_diff_df <- diff_dfs$auprc_diff_df

auprc_diff_df$method = factor(auprc_diff_df$method, levels=c("2to10", "2to10-sameObs", "1-kmSq", "lat-long", "rounded-4", "SVS", "1-UL", "best-clustGeo", "DBSC", "BayesOptClustGeo"))



# AUPRC difference plot (clustGeo = tuned)
png("results/species_experiments/plotted/plots/auprc_diff.png", units="in", width=5, height=6, res=400)


plot = ggplot(auprc_diff_df, aes(x = reorder(method, auprc_diff), y = auprc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=colors) +
    ylim(NA, 0.025) +
    theme(axis.text.x = element_text(angle = 55,vjust = 1, hjust=1)) +
    labs(y= paste0("AUPRC Improvement over ", benchmark_method), x="Algorithm", fill = "Algorithm") +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


# Percentage AUPRC difference plot (clustGeo = tuned)
png("results/species_experiments/plotted/plots/auprc_perc_diff.png", units="in", width=5, height=6, res=400)

plot = ggplot(auprc_diff_df, aes(x = reorder(method, auprc_perc_diff), y = auprc_perc_diff, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=colors) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= paste0("Percentage AUPRC Improvement over ", benchmark_method), x="Algorithm", fill = "Algorithm") +
    scale_y_continuous( breaks = c(-10, -7.5, -5, -2.5, 0, 2.5))  +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()

# Percentage AUPRC difference compared to species mean plot (clustGeo = tuned)
png("results/species_experiments/plotted/plots/auprc_perc_diff_species_mean.png", units="in", width=5, height=6, res=400)

plot = ggplot(auprc_diff_df, aes(x = reorder(method, auprc_perc_diff_from_species_mean), y = auprc_perc_diff_from_species_mean, fill=method)) +
    theme_classic() +
    geom_hline(yintercept=0, linetype='dashed', lwd=1, col = 'darkgrey') +
    geom_boxplot() +
    stat_summary(fun=mean, geom="point", shape=23, size=2, col="white", bg="darkred") +
    scale_fill_manual(values=colors) + 
    theme(axis.text.x = element_text(angle = 55, vjust = 1, hjust=1)) +
    labs(y= "Percentage Improvement over Species Mean AUPRC", x="Algorithm", fill = "Algorithm") +
    scale_y_continuous( breaks = c(-10, -7.5, -5, -2.5, 0, 2.5, 5.0))  +
    theme(axis.text=element_text(size=11), axis.title=element_text(size=12), legend.position="none")

show(plot)
dev.off()


sink("results/species_experiments/summarized/tests_auprc.txt")
cat("\nEffects of clustering algorithms on percentage AUPRC improvement over ", benchmark_method,"\n\n")


# Kruskal-wallis test
kruskal_test_auprc <- kruskal.test(auprc_perc_diff ~ method, data = auprc_diff_df) 
print(kruskal_test_auprc)

# Dunn's test
dunn_test_auprc <- dunnTest(auprc_perc_diff ~ method, data = auprc_diff_df, method = "bh")
print(dunn_test_auprc)
sink()

##############################################################
# # Mixed-effect models - AUPRC
##############################################################

auprc_diff_df_n <- diff_dfs$species_auprc_diff_df



auprc_diff_df_n <- merge(auprc_diff_df_n, results_sum_t, by.x = c("species","method","run"), by.y = c("species","method","run"))



auprc_diff_df_n$method <- as.factor(auprc_diff_df_n$method)
auprc_diff_df_n$prevalence <- as.factor(auprc_diff_df_n$prevalence)
auprc_diff_df_n$habitat <- as.factor(auprc_diff_df_n$habitat)
auprc_diff_df_n$specialist <- as.factor(auprc_diff_df_n$specialist)
auprc_diff_df_n$home_range <- as.factor(auprc_diff_df_n$home_range)


# print(unique(auprc_diff_df_n$method))
auprc_diff_df_n$method <- factor(auprc_diff_df_n$method, levels=c( "2to10-sameObs" , "SVS", "1-UL",  "lat-long", "rounded-4", "1-kmSq",  "DBSC", "BayesOptClustGeo", "best-clustGeo", "2to10"))

auprc_diff_df_n$prevalence <- factor(auprc_diff_df_n$prevalence, levels = c("low", "medium", "high"))
auprc_diff_df_n$habitat <- factor(auprc_diff_df_n$habitat, levels = c("forested", "seral"))
auprc_diff_df_n$specialist <- factor(auprc_diff_df_n$specialist, levels = c("generalist", "specialist"))
auprc_diff_df_n$home_range <- factor(auprc_diff_df_n$home_range, levels = c("small", "medium", "large"))



# Set reference levels of categorical predictors

auprc_diff_df_n$prevalence = relevel(auprc_diff_df_n$prevalence, ref="low")
auprc_diff_df_n$habitat = relevel(auprc_diff_df_n$habitat, ref="forested")
auprc_diff_df_n$specialist = relevel(auprc_diff_df_n$specialist, ref="generalist")
auprc_diff_df_n$home_range = relevel(auprc_diff_df_n$home_range, ref="small")
auprc_diff_df_n$method = relevel(auprc_diff_df_n$method, ref="2to10-sameObs")


auprc_diff_df_n$genspec = auprc_diff_df_n$specialist



lmer_algortihm <- lmer(auprc_perc_diff ~ method * + (1 | species), data = auprc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm <- summary(lmer_algortihm)

lmer_algorithm_prevalence <- lmer(auprc_perc_diff ~ method:prevalence + (1 | species), data = auprc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_prevalence <- summary(lmer_algorithm_prevalence)

lmer_algorithm_habitat <- lmer(auprc_perc_diff ~ method:habitat + (1 | species), data = auprc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_habitat <- summary(lmer_algorithm_habitat)

lmer_algorithm_specialist <- lmer(auprc_perc_diff ~ method:genspec + (1 | species), data = auprc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_specialist <- summary(lmer_algorithm_specialist)

lmer_algorithm_hrange <- lmer(auprc_perc_diff ~ method:home_range + (1 | species), data = auprc_diff_df_n, control=lmerControl(check.rankX="silent.drop.cols"))
summary_output_algorithm_hrange <- summary(lmer_algorithm_hrange)


# Custom function to format and print coefficient matrix with significance stars
get_coefmat <- function(coefmat, digits = 6) {
    # Round the coefficients matrix to the specified number of digits
    rounded_coefmat <- round(coefmat, digits = digits)
    
    # Add a column for significance stars
    p_values <- rounded_coefmat[, "Pr(>|t|)"]
    stars <- ifelse(p_values < 0.001, "***",
            ifelse(p_values < 0.01, "**",
            ifelse(p_values < 0.05, "*",
            ifelse(p_values < 0.1, ".", " "))))
    
    # Combine rounded coefficients and stars into a data frame
    df_coefmat <- as.data.frame(rounded_coefmat)
    df_coefmat$signif_code <- stars
    
    return(df_coefmat)
}

# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auprc_algorithm_mixed_coefficients.csv", row.names=TRUE)




# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_prevalence$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auprc_algorithm_prevalence_mixed_coefficients.csv", row.names=TRUE)


# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_habitat$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auprc_algorithm_habitat_mixed_coefficients.csv", row.names=TRUE)


# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_specialist$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auprc_algorithm_specialist_mixed_coefficients.csv", row.names=TRUE)


# Get the formatted coefficient matrix
df_coefmat <- get_coefmat(summary_output_algorithm_hrange$coefficients)

# Write the summarized results to CSV
write.csv(df_coefmat, "results/species_experiments/summarized/auprc_algorithm_homerange_mixed_coefficients.csv", row.names=TRUE)

# Plot mixed effect model results


plotCoefficients <- function(df_path, output_path, dimensions, title_str, rem_str, title_fontsize) {
    
    
    # Load the data
    coef_data <- read.csv(df_path, row.names = 1)
    
    # Remove the intercept term
    coef_data <- coef_data[-1,]
    
    coef_data$Coefficient <- rownames(coef_data)
    
    # Remove the word "method" from coefficient names
    coef_data$Coefficient <- gsub("method", "", coef_data$Coefficient)
    
    # Remove rem_str from coefficient names
    coef_data$Coefficient <- gsub(rem_str, "", coef_data$Coefficient)
    
    # Rename 'Std. Error' to 'Std_Error' to avoid space issues
    if ("Std..Error" %in% colnames(coef_data)) {
        colnames(coef_data)[colnames(coef_data) == "Std..Error"] <- "Std_Error"
    }
    
    # Reorder the data for plotting
    coef_data <- coef_data[order(coef_data$Estimate), ]
    
    # Define the possible significance levels and their properties
    signif_levels <- c("***", "**", "*", ".", " ")
    signif_colors <- c("***" = "red", "**" = "orange", "*" = "green", "." = "blue", " " = "black")
    signif_labels <- c("***" = "p-value < 0.001", "**" = "p-value < 0.01", "*" = "p-value < 0.05", "." = "p-value < 0.1", " " = "p-value >= 0.1")
    
    # Filter the significance levels present in the data
    present_signif_levels <- intersect(signif_levels, unique(coef_data$signif_code))
    
    # Calculate the number of rows for the legend dynamically
    legend_rows <- ceiling(length(present_signif_levels) / 2)
    
    # Create the coefficient plot
    plot <- ggplot(coef_data, aes(x = Estimate, y = reorder(Coefficient, Estimate))) +
        geom_errorbarh(aes(xmin = Estimate - Std_Error, xmax = Estimate + Std_Error), color = "darkgrey", height = 0.15) +
        geom_point(aes(color = signif_code), size = 2) +
        labs(title = paste0("", title_str),
             x = "Coefficient",
             y = "Effect") +
        scale_color_manual(
            name = "Significance",
            values = signif_colors[present_signif_levels],
            labels = signif_labels[present_signif_levels]
        ) +
        theme_bw() +
        theme(axis.text.y = element_text(size = 10, angle = 0, hjust = 1),
              legend.position = "bottom",
              legend.justification = "center", # Center the legend
              legend.box = "horizontal",
              legend.title = element_text(size = 10, hjust = 0),  # Align legend title to the left
              legend.text = element_text(size = 9, hjust = 0),   # Align legend text to the left
              plot.title = element_text(size = title_fontsize, hjust = 0)) +  # Move the title further to the left
        guides(color = guide_legend(nrow = legend_rows, byrow = TRUE))
    
    # Save the plot to a 300 dpi PNG file with specified dimensions
    ggsave(output_path, plot = plot, dpi = 300, width = dimensions[1], height = dimensions[2])
}




plotCoefficients(
    "results/species_experiments/summarized/auprc_algorithm_mixed_coefficients.csv", 
    "results/species_experiments/plotted/plots/mixed_coefficients_algorithms_auprc.png", 
    c(6,5),
    "% AUPRC Improvement ~ Algorithm",
    "",
    9
)


plotCoefficientsTraits <- function(df_paths, output_path, dimensions, title_strs, rem_strs, title_fontsize) {
    
    plots <- list()
    combined_signif_levels <- c()
    
    for (i in 1:4) {
        # Load the data
        coef_data <- read.csv(df_paths[i], row.names = 1)
        
        # Remove the intercept term
        coef_data <- coef_data[-1,]
        
        coef_data$Coefficient <- rownames(coef_data)
        
        # Remove the word "method" from coefficient names
        coef_data$Coefficient <- gsub("method", "", coef_data$Coefficient)
        
        # Remove rem_str from coefficient names
        coef_data$Coefficient <- gsub(rem_strs[i], "", coef_data$Coefficient)
        
        # Rename 'Std. Error' to 'Std_Error' to avoid space issues
        if ("Std..Error" %in% colnames(coef_data)) {
            colnames(coef_data)[colnames(coef_data) == "Std..Error"] <- "Std_Error"
        }
        
        # Reorder the data for plotting
        coef_data <- coef_data[order(coef_data$Estimate), ]
        
        # Define the possible significance levels and their properties
        signif_levels <- c("***", "**", "*", ".", " ")
        signif_colors <- c("***" = "red", "**" = "orange", "*" = "green", "." = "blue", " " = "black")
        signif_labels <- c("***" = "p-value < 0.001", "**" = "p-value < 0.01", "*" = "p-value < 0.05", "." = "p-value < 0.1", " " = "p-value >= 0.1")
        
        # Filter the significance levels present in the data
        present_signif_levels <- intersect(signif_levels, unique(coef_data$signif_code))
        combined_signif_levels <- unique(c(combined_signif_levels, present_signif_levels))
        
        # Create the coefficient plot
        p <- ggplot(coef_data, aes(x = Estimate, y = reorder(Coefficient, Estimate))) +
            geom_errorbarh(aes(xmin = Estimate - Std_Error, xmax = Estimate + Std_Error), color = "darkgrey", height = 0.15) +
            geom_point(aes(color = signif_code), size = 2) +
            labs(title = title_strs[i],
                 x = "Coefficient",
                 y = "Effect") +
            scale_color_manual(
                name = "Significance",
                values = signif_colors[signif_levels],
                labels = signif_labels[signif_levels]
            ) +
            theme_bw() +
            theme(axis.text.y = element_text(size = 15, angle = 0, hjust = 1),
                  axis.title.y = element_text(size = 15),
                  axis.text.x = element_text(size = 14),
                  axis.title.x = element_text(size = 14), 
                  plot.title = element_text(size = title_fontsize-2, hjust = 0),
                  legend.position = "bottom",
                  legend.box = "horizontal",
                  legend.title = element_text(size = 18, hjust = 0),
                  legend.text = element_text(size = 16, hjust = 0),
                #   legend.spacing.x = unit(0.5, 'in')
                  ) +
            guides(color = guide_legend(nrow = 1, byrow = TRUE))
        
        plots[[i]] <- p
    }
    
    # Extract the combined legend
    combined_legend <- ggplot(data.frame(signif_code = combined_signif_levels), aes(x = 1, y = 1, color = signif_code)) +
        geom_point() +
        scale_color_manual(
            name = "Significance",
            values = signif_colors[combined_signif_levels],
            labels = signif_labels[combined_signif_levels]
        ) +
        theme_void() +
        theme(legend.position = "bottom",
            legend.box = "horizontal",
            legend.title = element_text(size = 18, hjust = 0, margin = margin(r = 0.35, unit='in')),
            legend.text = element_text(size = 16, hjust = 0, margin = margin(t = 0, unit='in')),
            legend.key.spacing.x = unit(0.35, 'in'),
            legend.key.width = unit(0.5, 'in')) +
        guides(color = guide_legend(nrow = 1, byrow = TRUE))
    
    # Extract the legend from the combined legend plot
    g_legend <- function(a.gplot) {
        tmp <- ggplotGrob(a.gplot)
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
    }
    
    legend <- g_legend(combined_legend)
    
    # Remove legends from individual plots
    plots <- lapply(plots, function(p) p + theme(legend.position = "none"))
    
    # Combine the plots into a 2x2 grid
    combined_plot <- arrangeGrob(
        grobs = plots,
        ncol = 2,
        nrow = 2,
        # top = textGrob("Interactions of Species Traits, Clustering Algorithm, and AUPRC\n", gp = gpar(fontsize = title_fontsize)),
        bottom = legend
    )
    
    # Save the combined plot to a 300 dpi PNG file with specified dimensions
    ggsave(output_path, plot = combined_plot, dpi = 300, width = dimensions[1], height = dimensions[2])
}


plotCoefficientsTraits(
    c(
        "results/species_experiments/summarized/auprc_algorithm_prevalence_mixed_coefficients.csv", 
        "results/species_experiments/summarized/auprc_algorithm_homerange_mixed_coefficients.csv",
        "results/species_experiments/summarized/auprc_algorithm_habitat_mixed_coefficients.csv", 
        "results/species_experiments/summarized/auprc_algorithm_specialist_mixed_coefficients.csv"
    ),
    "results/species_experiments/plotted/plots/mixed_coefficients_traits_auprc.png",
    c(15, 14),
    c("Prevalence",
      "Home Range",
      "Habitat",
      "Generalist/Specialist"),
    c("prevalence", "home_range", "habitat", "genspec"),
    20
)

##############################################################


##############################################################
##############################################################
# Aggregated significance test results
##############################################################
##############################################################


# Kruskal-Wallis test and Dunn's test results

kruskal_wallis_df <- data.frame(
    algorithms = character(0),
    metric = character(0), 
    kruskal_wallis_chi_squared_statistic = numeric(0),
    p_value = numeric(0)
)

kruskal_wallis_df <- rbind(
    kruskal_wallis_df, 
    data.frame(
        algorithms = "best-clustGeo parameter combinations",
        metric = "% AUC improvement",
        kruskal_wallis_chi_squared_statistic = kruskal_test_cg_auc$statistic,
        p_value =  kruskal_test_cg_auc$p.value
    )
)

kruskal_wallis_df <- rbind(
    kruskal_wallis_df, 
    data.frame(
        algorithms = "Clustering algorithms",
        metric = "% AUC improvement",
        kruskal_wallis_chi_squared_statistic = kruskal_test_auc$statistic,
        p_value =  kruskal_test_auc$p.value
    )
)

kruskal_wallis_df <- rbind(
    kruskal_wallis_df, 
    data.frame(
        algorithms = "best-clustGeo parameter combinations",
        metric = "% AUPRC improvement",
        kruskal_wallis_chi_squared_statistic = kruskal_test_cg_auprc$statistic,
        p_value =  kruskal_test_cg_auprc$p.value
    )
)

kruskal_wallis_df <- rbind(
    kruskal_wallis_df, 
    data.frame(
        algorithms = "Clustering algorithms",
        metric = "% AUPRC improvement",
        kruskal_wallis_chi_squared_statistic = kruskal_test_auprc$statistic,
        p_value =  kruskal_test_auprc$p.value
    )
)

numeric_cols <- c("kruskal_wallis_chi_squared_statistic")
kruskal_wallis_df[,numeric_cols] <- round(kruskal_wallis_df[,numeric_cols], digits=3)


numeric_cols <- c("p_value")
kruskal_wallis_df[,numeric_cols] <- signif(kruskal_wallis_df[,numeric_cols], digits=3)


write.csv(kruskal_wallis_df, "results/species_experiments/summarized/kruskal_wallis_test_results.csv", row.names=FALSE)


# Dunn's test post-hoc results



dunn_df_cg_auc <- dunn_test_cg_auc$res
dunn_df_auc <- dunn_test_auc$res
dunn_df_cg_auprc <- dunn_test_cg_auprc$res
dunn_df_auprc <- dunn_test_auprc$res

dunn_df_cg_auc[,c("group_1","group_2")] <- str_split_fixed(dunn_df_cg_auc$Comparison, " - ", 2)
dunn_df_auc[,c("group_1","group_2")] <- str_split_fixed(dunn_df_auc$Comparison, " - ", 2)
dunn_df_cg_auprc[,c("group_1","group_2")] <- str_split_fixed(dunn_df_cg_auprc$Comparison, " - ", 2)
dunn_df_auprc[,c("group_1","group_2")] <- str_split_fixed(dunn_df_auprc$Comparison, " - ", 2)



cols_to_save <- c("group_1", "group_2", "Z", "P.unadj", "P.adj")


dunn_df_cg_auc <- dunn_df_cg_auc[,cols_to_save]
dunn_df_auc <- dunn_df_auc[,cols_to_save]
dunn_df_cg_auprc <- dunn_df_cg_auprc[,cols_to_save]
dunn_df_auprc <- dunn_df_auprc[,cols_to_save]

new_colnames <- c("group_1", "group_2", "z_test_statistic", "unadjusted_p_value", "benjamini_hochberg_adjusted_p_value")

colnames (dunn_df_cg_auc) <- new_colnames
colnames (dunn_df_auc) <- new_colnames
colnames (dunn_df_cg_auprc) <- new_colnames
colnames (dunn_df_auprc) <- new_colnames

numeric_cols <- c("z_test_statistic")
dunn_df_cg_auc[,numeric_cols] <- round(dunn_df_cg_auc[,numeric_cols], digits=3)
dunn_df_auc[,numeric_cols] <- round(dunn_df_auc[,numeric_cols], digits=3)
dunn_df_cg_auprc[,numeric_cols] <- round(dunn_df_cg_auprc[,numeric_cols], digits=3)
dunn_df_auprc[,numeric_cols] <- round(dunn_df_auprc[,numeric_cols], digits=3)

numeric_cols <- c("unadjusted_p_value", "benjamini_hochberg_adjusted_p_value")
dunn_df_cg_auc[,numeric_cols] <- signif(dunn_df_cg_auc[,numeric_cols], digits=3)
dunn_df_auc[,numeric_cols] <- signif(dunn_df_auc[,numeric_cols], digits=3)
dunn_df_cg_auprc[,numeric_cols] <- signif(dunn_df_cg_auprc[,numeric_cols], digits=3)
dunn_df_auprc[,numeric_cols] <- signif(dunn_df_auprc[,numeric_cols], digits=3)

alg_order <- c("SVS", "1-UL", "lat-long", "2to10", "2to10-sameObs", "rounded-4", "1-kmSq",  "best-clustGeo", "BayesOptClustGeo", "DBSC")

dunn_df_cg_auc$group_1 <- factor(dunn_df_cg_auc$group_1, levels=alg_order)
dunn_df_cg_auc$group_2 <- factor(dunn_df_cg_auc$group_2, levels=alg_order)
dunn_df_auc$group_1 <- factor(dunn_df_auc$group_1, levels=alg_order)
dunn_df_auc$group_2 <- factor(dunn_df_auc$group_2, levels=alg_order)
dunn_df_cg_auprc$group_1 <- factor(dunn_df_cg_auprc$group_1, levels=alg_order)
dunn_df_cg_auprc$group_2 <- factor(dunn_df_cg_auprc$group_2, levels=alg_order)
dunn_df_auprc$group_1 <- factor(dunn_df_auprc$group_1, levels=alg_order)
dunn_df_auprc$group_2 <- factor(dunn_df_auprc$group_2, levels=alg_order)

dunn_df_cg_auc <- dunn_df_cg_auc[order(dunn_df_cg_auc$group_1, dunn_df_cg_auc$group_2),]
dunn_df_auc <- dunn_df_auc[order(dunn_df_auc$group_1, dunn_df_auc$group_2),]
dunn_df_cg_auprc <- dunn_df_cg_auprc[order(dunn_df_cg_auprc$group_1, dunn_df_cg_auprc$group_2),]
dunn_df_auprc <- dunn_df_auprc[order(dunn_df_auprc$group_1, dunn_df_auprc$group_2),]


write.csv(dunn_df_cg_auc, "results/species_experiments/summarized/dunn_test_results_clustGeo_auc.csv", row.names=FALSE)
write.csv(dunn_df_auc, "results/species_experiments/summarized/dunn_test_results_auc.csv", row.names=FALSE)
write.csv(dunn_df_cg_auprc, "results/species_experiments/summarized/dunn_test_results_clustGeo_auprc.csv", row.names=FALSE)
write.csv(dunn_df_auprc, "results/species_experiments/summarized/dunn_test_results_auprc.csv", row.names=FALSE)


##############################################################



##############################################################
##############################################################
# 
##############################################################
##############################################################

plotClustGroupCG <- function (species, alg_all, extent="all"){

    nrow_c <- 4
    height_p <- 26

    png(paste0("results/species_experiments/plotted/maps/clust_maps_cg.png"), units="in", width=12, height=height_p, res=300)

    plots <- list()
    i <- 1
    for (m in alg_all){

        pts_df <- read.delim(paste0("results/species_experiments/raw/clusterings/", species,"_", m, ".csv"), sep = ",", header = T)
        pts_df$site <- as.factor(pts_df$site)
        pl <- ggplot() +
                geom_point(
                    data = pts_df, 
                    aes(
                        x=longitude, 
                        y=latitude,
                        fill=site, 
                        ), 
                    shape=21, 
                    size=1.5, 
                    color="black", 
                    show.legend = FALSE
                ) +
                theme_bw() +
                theme(plot.margin = margin(t=0.25,r=0.0,b=0,l=0.0, "in")) +
                theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
      

        plots[[i]] <- pl

        i <- i + 1

        
    }
    


    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=3, common.legend=TRUE, legend="bottom", labels=alg_all)
        
    
    show(plot_main)
    dev.off()

        
}


##############################################################





##############################################################
##############################################################
# SDM maps and cluster visualization
##############################################################
##############################################################



plotMapGroupCG <- function (species, pts_df, alg_all){

    nrow_c <- 4
    height_p <- 24

    png(paste0("results/species_experiments/plotted/maps/", species, "_cg_occu_maps.png"), units="in", width=12, height=height_p, res=300)

    plots <- list()
    i <- 1
    for (m in alg_all){

        
        perf_metrics <- read.delim(paste0("results/species_experiments/raw/model_parameters/", species,"_", m, ".csv"), sep = ",", header = T)
            
        par_lists <- get_model_parameter_list(perf_metrics, 1)

        occ.estimates_list <- par_lists$occ.estimates_list

        predicted_map <- predict_sdm_map(occ.estimates_list, OR.normalized)
        
        predicted_occu_prob <- predicted_map$occ_prob


        names(predicted_occu_prob) <- c("occu")

        predicted_occu_prob <- aggregate(predicted_occu_prob, fact=6)

        predicted_occu_prob <- as.data.frame(predicted_occu_prob, xy=TRUE)
        
        
        x_min = min(predicted_occu_prob$x)
        x_max = max(predicted_occu_prob$x)
        y_min = min(predicted_occu_prob$y)
        y_max = max(predicted_occu_prob$y)
        
        pl <- ggplot()+
            geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgrey", show.legend = FALSE)  +
            geom_raster(data = predicted_occu_prob, aes(x = x, y = y, fill = `occu`)) +
            scale_fill_viridis_c(option = "B", limits = c(0.0, 1.0)) +
            theme_void() +
            coord_fixed() +
            theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
            theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt")) +
            theme(
                legend.text = element_text(size = 14), 
                legend.title = element_text(size = 16, margin = margin(r=10, b=20)), 
                legend.key.size = unit(1, 'cm')
                ) +
            labs(fill = "Occupancy Probability")



        plots[[i]] <- pl

        i <- i + 1

        rm(pl)
        rm(predicted_occu_prob)
        
    }
    


    plot_main <- ggarrange(plotlist=plots, nrow=nrow_c, ncol=3, hjust = -1.0, common.legend=TRUE, legend="bottom", labels=alg_all)
        
    
    show(plot_main)
    dev.off()

        
}


plotMapGroupAll <- function (species, pts_df, alg_all){

   

    obs_plot <- NA
    plots <- list()
    i <- 1



    for (m in alg_all){

        
        perf_metrics <- read.delim(paste0("results/species_experiments/raw/model_parameters/", species,"_", m, ".csv"), sep = ",", header = T)
            
        par_lists <- get_model_parameter_list(perf_metrics, 1)

        occ.estimates_list <- par_lists$occ.estimates_list

        predicted_map <- predict_sdm_map(occ.estimates_list, OR.normalized)
        
        predicted_occu_prob <- predicted_map$occ_prob


        names(predicted_occu_prob) <- c("occu")

        predicted_occu_prob <- aggregate(predicted_occu_prob, fact=6)

        predicted_occu_prob <- as.data.frame(predicted_occu_prob, xy=TRUE)
        
     

        x_min = min(predicted_occu_prob$x)
        x_max = max(predicted_occu_prob$x)
        y_min = min(predicted_occu_prob$y)
        y_max = max(predicted_occu_prob$y)

        if (i==1){

            obs_plot <- ggplot()+ 
            
                geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgray", show.legend = FALSE)  +
                geom_raster(data = predicted_occu_prob, aes(x = x, y = y), fill="#E6E6E6", show.legend = FALSE) +
                geom_point(data = pts_df, aes(x=longitude, y=latitude, color=species_observed, shape=species_observed, fill=species_observed, size=species_observed), show.legend = TRUE) +
                scale_color_manual(name = "Observation", values=c("Detection" = "black", "Non-detection" = "black")) +
                scale_fill_manual(name = "Observation", values=c("Detection" = "#39FF14", "Non-detection" = "#83A1CD")) +
                scale_shape_manual(name = "Observation", values=c("Detection" = 24, "Non-detection" = 22)) +
                scale_size_manual(name = "Observation", values=c("Detection" = 1.6, "Non-detection" = 1.5)) +
                labs(title="Species Observations") +
                theme_void() +
                coord_fixed() +
                theme(
                    plot.title = element_text(hjust = 0.5, vjust = -56, face="bold", size = 15)
                ) +
                theme(
                    legend.position = "inside",
                    legend.position.inside = c(0.5, -0.16),
                    legend.direction="vertical",
                    legend.text = element_text(size = 14), 
                    legend.title = element_text(size = 15, margin = margin(l=20, b=10)), 
                    legend.spacing.x = unit(3, "cm"), 
                    legend.key.size = unit(1, 'cm')
                    ) +
                theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) +
                theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt")) 



        }


        pl <- ggplot() +
            geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgray", show.legend = FALSE)  +
            geom_raster(data = predicted_occu_prob, aes(x = x, y = y, fill = `occu`)) +
            scale_fill_viridis_c(option = "B", limits = c(0.0, 1.0)) +
            theme_void() +
            coord_fixed() +
            theme(plot.margin = margin(t = 0.25, r = 0, b = 0, l = 0, unit = "cm")) +
            theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(), axis.ticks.length = unit(0, "pt")) +
            theme(
                legend.text = element_text(size = 14), 
                legend.title = element_text(size = 15, margin = margin(r=10, b=20)), 
                legend.key.size = unit(1, 'cm')
                ) +
            labs(fill = "Occupancy Probability")


        
        plots[[i]] <- pl

        i <- i + 1

        rm(pl)
        rm(predicted_occu_prob)
        
    }
    
    
    alg_all[startsWith(alg_all, "clustGeo")] <- "best-clustGeo"

    label_x_pos <- c(0.0, 0.35, 0.2, 0.375, 0.4, 0.3, 0.3, 0.35, 0.05, -0.1)
    
    nrow_c <- 2
    

    plot_occu <- ggarrange(
        plotlist=plots, 
        nrow=nrow_c, 
        ncol=5, 
        common.legend=TRUE, 
        legend="bottom", 
        labels=alg_all,
        label.x=label_x_pos
    )

    plot_main <- 
        obs_plot + 
        plot_occu  + 
        plot_layout(nrow = 1, byrow = TRUE, widths = c(1, 5))

    height_p <- 9.5
    width_p <- 17
      
    png(paste0("results/species_experiments/plotted/maps/", species, "_occu_maps.png"), units="in", width=width_p, height=height_p, res=300)

    show(plot_main)
    dev.off()

        
}


plotClustGroup <- function (species, alg_all){

    nrow_c <- 1
    

    zoomed_A = list(
        longitude = c(-123.305, -123.29175),
        latitude = c(44.4125, 44.425)
    )
    obs_plot <- NA
    plots <- list()

    valid_boundary <- terra::vect("occupancy_feature_raster/boundary/boundary.shp")
    crs(valid_boundary) <- crs(OR.normalized)
    
    region <- terra::crop(OR.normalized$elevation, valid_boundary, mask = TRUE)

    base_rast <- as.data.frame(region, xy=TRUE)
    
    

    x_min = min(base_rast$x)
    x_max = max(base_rast$x)
    y_min = min(base_rast$y)
    y_max = max(base_rast$y)

    i <- 1
    for (m in alg_all){


        pts_df <- read.delim(paste0("results/species_experiments/raw/clusterings/", species,"_", m, ".csv"), sep = ",", header = T)
        pts_df$site <- as.factor(pts_df$site)
        n_sites <- length(unique(pts_df$site))


        if (i==1){

            obs_plot <- ggplot()+ 
            
                geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max), fill = "darkgray", show.legend = FALSE)  +
                geom_raster(data = base_rast, aes(x = x, y = y), fill="#E6E6E6", show.legend = FALSE) +
                geom_point(data = pts_df, aes(x=longitude, y=latitude), color="black", size= 1, shape= 21, fill = "#FBFAF5", show.legend = FALSE) +
                annotate(
                    "segment", 
                    x = -123.75, 
                    xend = zoomed_A$longitude[1], 
                    y = 44.425, 
                    yend = (zoomed_A$latitude[1] + zoomed_A$latitude[2])/2, 
                    colour = "blue"
                    ) +
                annotate(
                    "label", 
                    x = -124.05, 
                    y = 44.425, 
                    label = "Zoomed-in region",
                    size = 3.25,
                    ) +
                geom_rect(aes(xmin = zoomed_A$longitude[1], xmax = zoomed_A$longitude[2], ymin = zoomed_A$latitude[1], ymax = zoomed_A$latitude[2]), fill = "red", color="red", alpha=0.1, show.legend = FALSE)  +
                labs(
                    title = "Species\nObservations",
                    x = "Longitude",
                    y = "Latitude"
                ) +

                theme_bw() +
                coord_fixed() +
                theme(
                    plot.title = element_text(hjust = 0.5, vjust = -25, face="bold", size = 12)
                ) +
                theme(plot.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) 
                


        }

        pts_df <- pts_df[
            (pts_df$longitude > zoomed_A$longitude[1]) &
            (pts_df$longitude < zoomed_A$longitude[2]) &
            (pts_df$latitude > zoomed_A$latitude[1]) &
            (pts_df$latitude < zoomed_A$latitude[2]) ,]

        pts_df$lat_long <- paste0(pts_df$latitude, "_", pts_df$longitude)
        pts_df <- dplyr::distinct(pts_df, lat_long, .keep_all = T)
        

        # print(nrow(pts_df))
        # print(length(unique(pts_df$site)))

        plots[[i]] <- ggplot() +
                        geom_rect(
                            aes(xmin = zoomed_A$longitude[1], 
                            xmax = zoomed_A$longitude[2], 
                            ymin = zoomed_A$latitude[1], 
                            ymax = zoomed_A$latitude[2]), 
                            fill = "#E6E6E6", 
                            color="#E6E6E6", 
                            show.legend = FALSE
                            ) +
                        geom_point(
                            data = pts_df, 
                            aes(
                                x=longitude, 
                                y=latitude,
                                fill=site 
                                ), 
                            shape = 21,
                            size = 3.5, 
                            color = "black",
                            show.legend = FALSE
                        ) +
                        geom_mark_ellipse(
                            data = pts_df,
                            aes(
                                x = longitude,
                                y = latitude, 
                                fill = site, 
                                color = site ),
                            alpha = 0.1,
                            show.legend = FALSE
                            ) +
                        theme_bw() +
                        coord_fixed() +
                        labs(
                            x = "Longitude",
                            y = "Latitude"
                        ) +
                        theme(plot.margin = margin(t=0.25,r=0.1,b=0.1,l=0.1, "in"))

        


        i <- i + 1

        
    }
    


    alg_all[startsWith(alg_all, "clustGeo")] <- "best-clustGeo"
 
    label_x_pos <- c( 0.35, 0.4, 0.45, 0.25, 0.125)
    

    plot_clust <- ggarrange(
        plotlist=plots, 
        nrow=2, 
        ncol=3, 
        labels=alg_all,
        label.x=label_x_pos,
        font.label = list(size = 12, face = "bold", vjust = -1)
         
    )

    plot_main <- 
        obs_plot + 
        plot_clust  + 
        plot_layout(nrow = 1, byrow = TRUE, widths = c(1, 4.25))

 

    height_p <- 6
    width_p <- 12

    png(paste0("results/species_experiments/plotted/maps/",species,"_clust_maps.png"), units="in", width=width_p, height=height_p, res=300)
    show(plot_main)
    dev.off()


        
}


# Plot all grouped maps and clusters

# species_all <- c("NOFL")

for (s in species_all){

    train.df.og <- read.delim(paste0("checklist_data/species/", s, "/", s, "_zf_filtered_region_2017.csv"), sep = ",", header = T)
        
    train.df.og <- train.df.og[!is.na(train.df.og$duration_minutes),]

    train.df.og <- train.df.og[train.df.og$observation_date >= "2017-05-15" & train.df.og$observation_date <= "2017-07-09",]

    test.df.og <- read.delim(paste0("checklist_data/species/", s, "/", s, "_zf_filtered_region_2018.csv"), sep = ",", header = T)
        
    test.df.og <- test.df.og[!is.na(test.df.og$duration_minutes),]

    test.df.og <- test.df.og[test.df.og$observation_date >= "2018-05-15" & test.df.og$observation_date <= "2018-07-09",]

    test.df.og <- test.df.og[,colnames(train.df.og)]
    pts_df <- rbind(train.df.og, test.df.og)

    pts_df[pts_df$species_observed==TRUE,]$species_observed <-  "Detection"
    pts_df[pts_df$species_observed==FALSE,]$species_observed <-  "Non-detection"
    pts_df$species_observed = factor(pts_df$species_observed, levels=c("Non-detection", "Detection"))


    pts_df <- pts_df %>%
        dplyr::arrange(species_observed)


    alg_all_cg <- alg_all[startsWith(alg_all, "clustGeo")]

    s_cgt <- tuned_cg_df[tuned_cg_df$species==s,]$method

    alg_all_t <- alg_all_o

    alg_all_t[alg_all_t == "clustGeo"] <- s_cgt


    plotMapGroupCG(s, pts_df, alg_all_cg)
    plotMapGroupAll(s, pts_df, alg_all_t)

    alg_clust <- alg_all_t[startsWith(alg_all_t, "clustGeo") | alg_all_t %in% c("rounded-4","1-kmSq","DBSC","BayesOptClustGeo")]
    plotClustGroup(s, alg_clust)

}



##############################################################


