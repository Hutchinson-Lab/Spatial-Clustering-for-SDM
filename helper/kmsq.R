####################
# 1 km sq

# December 10, 2024
####################

library(sf) # spatial geometry operations
library(auk) # filter_repat_visits function

MIN_OBS <- 1
MAX_OBS <- 100000



# ##########
# Function for kmSq method
# Overlays a grid of squares each of side length rad_m meters
# Points which fall in the same square/cell are in same sites
# Takes species_df, filters it, overlays grid cell polygons 
#   to check where each row in species_df falls in and assigns 
#   a new column "site" to checklists_filtered based on id of polygon it falls in 
# ##########

kmsq.Sites.Core <- function(species_df, rad_m, filter = TRUE) {
   
    if (filter) {
        checklists_filtered <- filter_repeat_visits(
            species_df,
            min_obs = MIN_OBS,
            max_obs = MAX_OBS,
            annual_closure = TRUE,
            date_var = "formatted_date",
            site_vars = c("locality_id")
        )
        checklists_filtered <- subset(checklists_filtered, select = -c(site, n_observations, closure_id))
    } else {
        checklists_filtered <- species_df
    }

    
    og.crs <- st_crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

    # Convert to sf object
    species.sf <- st_as_sf(checklists_filtered, coords = c("longitude", "latitude"), crs = og.crs, remove = FALSE)

    # Transform to Albers CRS
    Albers.crs <- st_crs("+proj=aea +lat_1=42 +lat_2=48 +lon_0=-122 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    species.albers <- st_transform(species.sf, crs = Albers.crs)
    
    # Get bounding box coordinates
    bbox <- st_bbox(species.albers)
    
    # Define boundary and resolution of grid
    x <- seq(from = bbox["xmin"] - (rad_m * 5), to = bbox["xmax"] + (rad_m * 5), by = rad_m)
    y <- seq(from = bbox["ymin"] - (rad_m * 5), to = bbox["ymax"] + (rad_m * 5), by = rad_m)
    
    # Create a grid of points
    grid_points <- expand.grid(x = x, y = y)
    grid_sf <- st_as_sf(grid_points, coords = c("x", "y"), crs = Albers.crs)
    
    # Create grid polygons
    grid <- st_make_grid(grid_sf, cellsize = rad_m, square = TRUE)
    grid_sf <- st_as_sf(grid)
    grid_sf$id <- seq_len(nrow(grid_sf))
    
    # Perform spatial join
    joined <- sapply(st_intersects(species.albers, grid_sf), function(z) if (length(z)==0) NA_integer_ else z[1])
    
    checklists_df <- st_drop_geometry(species.albers)
    checklists_df$site <- as.character(unlist(joined))
    
    return(checklists_df)
}







kmsq.Sites <- function(species_df, rad_m){
    checklists <- kmsq.Sites.Core(species_df, rad_m, filter = TRUE)
    return(checklists)
}
