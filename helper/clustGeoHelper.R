######################
# ClustGeo helper

# December 10, 2024
######################

library(ClustGeo) # spatial clustering with clustGeo


clustGeoSites <- function(alpha, checklists, occ_covs, det_covs, num_sites=NULL, ratio=NULL, filter = TRUE){
    
    checklists_filtered <- checklists    
    
    if(is.null(num_sites)){
        num_sites = round(nrow(checklists_filtered)/ratio)
    }
    
    env_data <- dist(subset(checklists_filtered, select = unlist(occ_covs)))
    geo_data <- dist(subset(checklists_filtered, select = c("latitude", "longitude")))
    
    tree <- hclustgeo(env_data, geo_data, alpha = alpha)
    part <- cutree(tree, num_sites)
    checklists_filtered$site <- part
    
    checklists <- data.frame(checklists_filtered)
    return(checklists)
}




