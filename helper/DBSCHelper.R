####################
# DBSC helper

# December 10, 2024
####################


library(sqldf) # sql operation on dataframes
library(RTriangle) # Triangle geometry operations
library(igraph) # graph node and edge operations


# extracts and formats vertex points
formatVert <- function(species_df, occ_covs){
    v_pts <- subset(species_df, select = c("latitude", "longitude", "checklist_id", occ_covs))
 
    v_pts$checklist_id <- as.character(v_pts$checklist_id)
    v_pts$latlong <- paste0(v_pts$latitude, v_pts$longitude)
    
    v_pts <- sqldf("SELECT * from v_pts group by latlong")
    v_pts <- subset(v_pts, select=-c(latlong))
    return(v_pts)
}




# calculates the global mean and standard deviation
# of the edges in a Delaunay Triangulation
# note: distance includes geo-coords
calcGlobalMeanSD <- function(DT, attr_df){
    edge_wts <- c()
    attr_df_no_checklist <- subset(attr_df, select = -c(checklist_id))
    for(i in 1:(length(DT$E)/2)){
        e <- DT$E[i,]
        v1 <- e[1]
        v2 <- e[2]
        
        pt1 <- attr_df[v1,]
        pt2 <- attr_df[v2,]
        
        d <- dist(rbind(attr_df_no_checklist[v1,], attr_df_no_checklist[v2,]))[1]
        edge_wts <- rbind(edge_wts, list(v1_name=pt1$checklist_id,v2_name=pt2$checklist_id,d=d))
    }
    if(nrow(edge_wts) == 1){
        sd <- 0
    } else {
        sd <- sd(as.double(as.data.frame(edge_wts)$d))
    }
    
    return(list(mean=mean(as.double(as.data.frame(edge_wts)$d)), sd=sd, edgeWts=edge_wts))
}

# removes Global Long Edges from a DT to form C-DT
##################################################
# for each vertex v in DT this func calculates the 
# local mean (mean length of edges incident to v)
# and determiens if any edge of v, dist(e_v) > gdc
# where:
#     gdc = globalmean + (globalmean/localmean(v))*globalSD
# if so, e_v is removed from DT
removeGlobalLongEdges <- function(DT, globalMean, globalSD, edgeWts, DT_graph){
    gdt.v <- igraph::as_data_frame(DT_graph, what = "vertices") #%>% setDT()
    edg_att <- edge_attr(DT_graph)
    
    to_delete_edge_ids <- c()
    for(v in V(DT_graph)$name){
        # nghbrs <- neighbors(DT_graph, v)$name
        in_edges_of_v <- edgeWts[edgeWts$v1_name == v,]
        out_edges_of_v <- edgeWts[edgeWts$v2_name == v,]
        
        edges_of_v <- rbind(in_edges_of_v, out_edges_of_v)
        localMean <- mean(as.double(edges_of_v$d))
        
        gdc <- globalMean + (globalMean/localMean)*globalSD
        for(e in 1:nrow(edges_of_v)){
            v1_name <- edges_of_v[e,]$v1_name
            v2_name <- edges_of_v[e,]$v2_name
            d <- edges_of_v[e,]$d
            
            if(d > gdc){
                edge_id <- get_edge_ids(DT_graph, vp = c(as.character(v1_name), as.character(v2_name)), error = FALSE)
                if(!(edge_id %in% to_delete_edge_ids)){
                    to_delete_edge_ids <- append(to_delete_edge_ids, edge_id)
                }
            }
        }
    }
    DT_graph <- DT_graph - edge(to_delete_edge_ids)
    return(list(graph=DT_graph, nRmv=length(to_delete_edge_ids)))
}

# local variation is the sd of the length of all
# edges incident to a vertex v
################################################
# given a vertex name and graph, this func finds
# v's neighbors and adds the distance to a list
# then returns the sd. returns0 if less than 2 
# items
local_variation <- function(g, v_name, gdt.v, edg_att){
    nghbrs <- neighbors(g, v_name)$name
    edge_val <- as.double()
    x <- 0
    for(n1 in nghbrs){
        x <- x + 1
        edge_id <- get_edge_ids(g, vp = c(v_name,n1))
        if(edge_id){
            d <- edg_att$d[[edge_id]]
            edge_val <- append(edge_val, d)
        }
    }
    
    if(x < 2){
        return(0)
    } else {
        return(sd(edge_val))
    }
}

# removes long local edges from each disconnected
# graph in C-DT
#################################################
# for each graph G_i in C-DT, this func iterates 
# through each vertex, v, and creates its 2-order 
# neighborhood. calculate ldc and determine if
# any edge in 2neighborhood is greater than ldc.
# if so, delete it from G_i.
#     ldc = 2-order-mean(v) + BETA*mean_variation(v)
#
# where the 2-order-mean is the mean of all edge
# dist in 2neighborhood and the mean variation is 
# the mean value of all local_variation calculations 
# for each vertex in 2neighborhood
removeLongLocalEdges <- function(de, Beta=2){
    list_of_graphs <- list()
    j <- 0
    for(g in de){
        edg_att <- igraph::edge.attributes(g)
        to_delete_edge_ids <- c()
        j <- j + 1
        gdt.v <- igraph::as_data_frame(g, what = "vertices") %>% setDT()
        if(length(gdt.v$name) > 1){
            paths_in_2 <- make_ego_graph(g, order=2)
            
            for(v_nghbrs in paths_in_2){
                edges_in_2 <- as_edgelist(v_nghbrs)
                list_dist_2_order <- c()    
                for(e in 1:nrow(edges_in_2)){
                    vertices_of_edge <- edges_in_2[e,]
                    edge_id <- get_edge_ids(g, vp = vertices_of_edge)
                    
                    # if we found the edge
                    if(edge_id != 0){
                        edge_val <- edg_att$d[[edge_id]]
                        list_dist_2_order <- append(list_dist_2_order, edge_val)
                    }
                }
                
                near_v <- unique(as.character(edges_in_2))
                local_variation_list <- c()
                for(v in near_v){
                    local_variation_list <- append(local_variation_list, local_variation(g, v, gdt.v, edg_att))
                }
                
                mean_2_order <- mean(list_dist_2_order)
                mean_variation <- mean(local_variation_list)
                
                if(is.na(mean_2_order)){
                    mean_2_order = 100000
                    disp("broken mean 2 order")
                }
                if(is.na(mean_variation)){
                    mean_variation = 100000
                    disp("broken mean variation")
                }
                
                ldc <- mean_2_order + Beta * mean_variation
                
                # if length of a 2-order edge > ldc -----> delete that mf
                for(i in 1:nrow(edges_in_2)){
                    e <- edges_in_2[i,]
                    edge_id <- get_edge_ids(g, vp = e)
                    
                    if(edge_id != 0){
                        e_dist <- edg_att$d[[edge_id]]
                        
                        if(e_dist > ldc){
                            to_delete_edge_ids <- append(to_delete_edge_ids, edge_id)
                        }    
                    }
                }
            }
        }
        g <- g - edge(to_delete_edge_ids) 
        list_of_graphs[[j]] <- g
    }
    return(list_of_graphs)
}

# calculates T1 from C-DT
#################################################
# for each graph G_i in C-DT, for every vertex v
# in G_i, find v's min neighbor in the spatial
# domain.    detect and remove outliers by the rule
# of 3 SDs from the min_edge list. return the 
# average difference
calcT1 <- function(de){
    min_edge_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(min_edge_df) <- c("edge_id", "edge_dist")
    for(g in de){
        
        v_df <- igraph::as_data_frame(g, what = "vertices") %>% setDT()
        edg_att <- igraph::edge.attributes(g)
        
        if(length(edg_att$d) != 0){
            for(v_name in v_df$name){
                # idx <- which(v_df$name == v_name)
                nghbrs <- neighbors(g, v_name)$name
                
                min_edge_wt <- Inf
                min_edge_id <- NULL
                for(n_1 in nghbrs){
                    # idx2 <- which(gdt.v$name == n)
                    e_id <- get_edge_ids(g, vp = c(v_name, n_1))    
                    e_dist <- edg_att$d[[e_id]]
                    if(e_dist < min_edge_wt){
                        min_edge_wt <- e_dist
                        min_edge_id <- e_id
                    }
                }
                if(!is.null(min_edge_id)){
                    min_edge_df <- rbind(min_edge_df, data.frame(edge_id=min_edge_id,edge_dist=min_edge_wt))
                }
            }
            
        }
    }
    
    sd <- sd(as.double(min_edge_df$edge_dist))
    edges_to_drop <- c()
    count <- 0
    for(i in 1:nrow(min_edge_df)){
        ed <- min_edge_df[i,]
        if(ed$edge_dist >= 3*sd){
            count <- count + 1
            edges_to_drop <- append(edges_to_drop, i)
        }
    }     
    # disp("edges dropped: ", as.character(count))
    min_edge_df <- min_edge_df[-c(edges_to_drop),]
    return(mean(min_edge_df$edge_dist))
    
}



isExpandingCore <- function(v1_name, g, T1){
    nghbrs <- neighbors(g, v1_name)$name
    for(n1_name in nghbrs){
        if(isSpatiallyDirectlyReachable(v1_name, n1_name, g, T1)){
            return(TRUE)
        }
    }
    return(FALSE)
}



isSpatiallyDirectlyReachable <- function(v1_name, v2_name, g, T1){
    nghbrs <- neighbors(g, v1_name)$name
    edg_att <- edge_attr(g)
    
    edge_id <- get_edge_ids(g, vp = c(v1_name, v2_name))
    
    d <- edg_att$d[edge_id]
    
    if(v2_name %in% nghbrs){
        if(d <= T1){
            return(TRUE)
        }
    }
    return(FALSE)
}



calcDensityIndicator <- function(v_name, g, v_df, T1){
    # idx <- which(v_df$name == v_name)
    nghbrs <- neighbors(g, v_name)$name
    num_nghbrs <- length(nghbrs)
    
    edg_att <- igraph::edge.attributes(g)
    
    num_sdr <- 0
    for(v2_name in nghbrs){
        # is sdr
        e_id <- get_edge_ids(g, vp = c(v_name, v2_name))    
        d <- edg_att$d[[e_id]]
        if(d <= T1){
            num_sdr <- num_sdr + 1
        }
    }
    
    return(num_sdr + num_sdr/num_nghbrs)
}

# calculate DI for each vertex and sort
# for each graph. return a sorted list on
# density indicator
calcDI.DF <- function(new_graphs, T1){
    DI_df <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(DI_df) <- c("v_name", "density_indicator")
    g_num <- 0
    for(g in new_graphs){
        g_num <- g_num + 1
        # for each vertex
        v_df <- igraph::as_data_frame(g, what = "vertices") %>% setDT()
        for(i in 1:nrow(v_df)){
            di <- calcDensityIndicator(v_df[i,]$name, g, v_df, T1)
            
            if(is.nan(di)){
                di <- 0
            }
            
            DI_df <- rbind(DI_df, data.frame(v_name=v_df[i,]$name, density_indicator=di, graph=g_num))
        }
    }
    return(DI_df[order(-DI_df$density_indicator),])
}


# selecting spatial core
select.SPCC <- function(sorted_DI, new_graphs){
    if(max(sorted_DI$density_indicator) == 0){
        return(NULL)
    }
    max_DI <- sorted_DI[sorted_DI$density_indicator == max(sorted_DI$density_indicator),]
    if(length(max_DI) > 1){
        
        min_avg_diff <- Inf
        min_spcc <- NULL
        for(i in 1:nrow(max_DI)){
            v <- max_DI[i,]
            edge_val <- c()
            v_name <- as.character(v$v_name)
            gra <- new_graphs[[v$graph]]
            
            edg_att <- edge_attr(gra)
            nghbrs <- neighbors(gra, v_name)$name
            
            # TODO: get_edge_ids(g, vp=c(v,n1,v,n2,v,n3 ...))
            for(n1 in nghbrs){
                edge_id <- get_edge_ids(gra, vp = c(v_name, as.character(n1)))
                edge_val <- append(edge_val, edg_att$d[[edge_id]])
            }
            
            if(mean(edge_val) < min_avg_diff){
                min_avg_diff <- mean(edge_val)
                min_spcc <- v
            }
        }
    }
    return(min_spcc)
}


isSpatiallyReachable <- function(v_name, CLU_names, g, v_pts, T1){
    
    v_df <- igraph::as_data_frame(g, what = "vertices") %>% setDT()
    
    nghbrs <- neighbors(g, v_name)$name
    v_pts_no_checklist <- subset(v_pts, select = -c(checklist_id))
    # CLU_idx <- which(v_df$name %in% CLU_v)
    CLU_names[[length(CLU_names)+1]] <- "-1"
    nghbrs[[length(nghbrs)+1]] <- "-1"
    # CLU_idx <- which(v_df$name %in% CLU_v)
    if(length(intersect(CLU_names, nghbrs)) > 1){
        # CLU_df <- as.data.frame(v_df[v_df$name %in% CLU_v])
        CLU_names <- CLU_names[1:(length(CLU_names)-1)]
        nghbrs <- nghbrs[1:(length(nghbrs)-1)]    
        v_info <- v_pts_no_checklist[v_pts$checklist_id == v_name,]
        avg_att <- colMeans(v_pts_no_checklist[v_pts$checklist_id %in% CLU_names,])
        
        
        if(dist(rbind(avg_att, v_info))[1] <= T1){
            return(TRUE)
        }
    } 
    return(FALSE)
}


################
# this is an implementation of Step 3 of the DBSC alogrithm.    
# First, we calculate the density indicator for every vertex.    
################
# (0): We find the spatial clustering core, spcc with largest 
# DI and unlabeled
################
# (i): rank the expanding core neighbors of spcc. create 
# clust and add spcc to it
################
# (ii): add these expanding cores in order of DI; every 
# neighbor must be both sdr from spcc and sr from clust
################
# (iii): look at each subsequent bfs level. the first expanding 
# core to be added is treated as the spcc and we go back to (i)
################
# when no more expanding cores can be added, a cluster is found. 
# go back to Step (0).
################
################
createClusters <- function(de_new_graphs, v_pts, T1){
    sorted_DI <- calcDI.DF(de_new_graphs, T1)
    sorted_DI$clust <- -1
    m_spcc <- select.SPCC(sorted_DI, de_new_graphs)
    clust_count <- 1
    
    s_DI <- sorted_DI
    while(!is.null(m_spcc)){
        spatialCore <- m_spcc
        spatial_core_name <- as.character(spatialCore$v_name)
        gra <- de_new_graphs[[spatialCore$graph]]
        
        clust <- c(spatial_core_name)
        bfs <- bfs(gra, root = spatial_core_name, unreachable = F, dist = T)
        for(k in 1:max(bfs$dist)){
            k_order_nghbrs <- rownames(as.data.frame(bfs$dist[bfs$dist == k]))
            nghbrs_sorted <- as.character(s_DI[s_DI$v_name %in% k_order_nghbrs,]$v_name)
            
            # Step (iii)/(i)
            if(k > 1){
                spatial_core_name <- NULL
                for(n1 in nghbrs_sorted){
                    
                    if(isSpatiallyReachable(n1, clust, gra, v_pts, T1)){
                        if(isExpandingCore(n1, gra, T1)){
                            spatial_core_name <- n1
                            nghbrs_sorted <- nghbrs_sorted[nghbrs_sorted != spatial_core_name]
                            clust <- append(clust, spatial_core_name)
                            break
                        }
                    }
                }
            }
            if(is.null(spatial_core_name)){
                break
            }
            
            # Step (ii)
            # find which in nghbrs_sorted are exp_core; add these if sr & sdr to clust
            for(n1 in nghbrs_sorted){
                if(!(n1 %in% clust)){
                    if(isExpandingCore(n1, gra, T1)){
                        if(isSpatiallyReachable(n1, clust, gra, v_pts, T1)){
                            if(isSpatiallyDirectlyReachable(n1, spatial_core_name, gra, T1)){
                                clust <- append(clust, n1)    
                            }
                        }
                    }
                }
            }
            
            # the issue is that the clust is never removed from sorted_DI 
        }
        sorted_DI[sorted_DI$v_name %in% clust,]$clust <- clust_count
        clust_count <- clust_count + 1
        
        s_DI <- sorted_DI[sorted_DI$clust == -1,]
        m_spcc <- select.SPCC(s_DI, de_new_graphs)
    }
    return(list(sorted_DI, s_DI))
}







runDBSC <- function(species_df, occ_covs, det_covs){
    v_pts <- formatVert(species_df, occ_covs)
    
    DT = RTriangle::triangulate(RTriangle::pslg(v_pts[c("latitude", "longitude")]))
    
    globalMeanDT <- calcGlobalMeanSD(DT, v_pts)
    globalMean <- globalMeanDT$mean
    globalSD <- globalMeanDT$sd
    edgeWts <- as.data.frame(globalMeanDT$edgeWts)
    
    DT_graph <- igraph::graph_from_data_frame(edgeWts, directed = F)
    DT_subgraphs <- removeGlobalLongEdges(DT, globalMean, globalSD, edgeWts, DT_graph)
    de <- decompose(DT_subgraphs$graph)
    
    new_graphs <- removeLongLocalEdges(de)
    de_new_graphs <- list()
    g_num <- 1
    for(n_g in new_graphs){
        decomp <- decompose(n_g)
        for(d in decomp){
            de_new_graphs[[g_num]] <- d
            g_num <- g_num + 1
        }
    }
    
    T1 <- calcT1(de_new_graphs)
    
    # step 3; make clusters & assign noise to separate cluster
    clust_li <- createClusters(de_new_graphs, v_pts, T1)
    sorted_DI <- as.data.frame(clust_li[1])
    max_clust <- max(sorted_DI$clust) + 1
    for(i in 1:nrow(sorted_DI)){
        pt <- sorted_DI[i,]
        if(pt$clust == -1){
            sorted_DI[i,]$clust <- max_clust
            max_clust <- max_clust + 1
        }
    }
    
    
    # assign clusters to checklists with the same lat/long
    species_df$clust <- -1
    for(j in 1:nrow(sorted_DI)){
        c_id <- as.character(sorted_DI[j,]$v_name)
        lat_long <- v_pts[v_pts$checklist_id == c_id,][c("latitude", "longitude")]
        species_df[species_df$latitude == lat_long$latitude & species_df$longitude == lat_long$longitude,]$clust <- sorted_DI[j,]$clust
    }
    
    species_df$site <- species_df$clust
    return(species_df)
}







