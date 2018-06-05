read_hdata <- function(filename) {
  
  ext_tracks_widths <- c(7, 10, 2, 2, 3, 5, 5, 6, 4, 5, 4, 4, 5, 3, 4, 3, 3, 3,
                         4, 3, 3, 3, 4, 3, 3, 3, 2, 6, 1)
  ext_tracks_colnames <- c("storm_id", "storm_name", "month", "day",
                           "hour", "year", "latitude", "longitude",
                           "max_wind", "min_pressure", "rad_max_wind",
                           "eye_diameter", "pressure_1", "pressure_2",
                           paste("radius_34", c("ne", "se", "sw", "nw"), sep = "_"),
                           paste("radius_50", c("ne", "se", "sw", "nw"), sep = "_"),
                           paste("radius_64", c("ne", "se", "sw", "nw"), sep = "_"),
                           "storm_type", "distance_to_land", "final")
  
  ext_tracks <- read_fwf("ebtrk_atlc_1988_2015.txt", 
                         fwf_widths(ext_tracks_widths, ext_tracks_colnames),
                         na = "-99")
  
  return(ext_tracks)
  
}

tidy_hdata <- function(data) {
  data %>% 
    
    
    dplyr::mutate_(storm_id = ~stringr::str_c(stringr::str_to_title(storm_name), year, sep = '-'),
                   date = ~stringr::str_c(year, '-', month, '-', day, ' ', hour, ':', '00', ':', '00'),
                   longitude = ~-longitude
    ) %>% 
    # Select only the relevant columns
    dplyr::select_(.dots = c('storm_id', 'date', 'longitude', 'latitude', 
                             'radius_34_ne', 'radius_34_se', 'radius_34_sw', 'radius_34_nw',
                             'radius_50_ne', 'radius_50_se', 'radius_50_sw', 'radius_50_nw',
                             'radius_64_ne', 'radius_64_se', 'radius_64_sw', 'radius_64_nw')
    ) %>%
    
    #There is a better way to do this part, this is the wide to long transfmration
    tidyr::gather(variable, value, -storm_id, -date,-latitude, -longitude, -storm_id, -date) %>% mutate_(wind_speed = ~str_extract(variable, "(34|50|64)"),
                                                                                                         variable = ~str_extract(variable, "(ne|nw|se|sw)")) %>% tidyr::spread(variable, value) %>% select_(.dots = c('storm_id', 'date', 'latitude', 'longitude', 'wind_speed', 'ne', 'nw', 'se', 'sw'))
}

filter_hdata <- function(data, hurricane, observation) {
  data <- filter_(data, ~storm_id == hurricane & date == observation)
  
}

geom_hurricane <- function(mapping = NULL, data = NULL, stat = 'identity',
                           position = 'identity', na.rm = FALSE,
                           show.legend = NA, inherit.aes = TRUE, ...) {
  ggplot2::layer(
    geom = geom_hurricane_proto, mapping = mapping,
    data = data, stat = stat, position = position, 
    show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
    
  )
}

geom_hurricane_proto <- ggplot2::ggproto("geom_hurricane_proto", Geom,
                                         required_aes = c("x", "y",
                                                          "r_ne", "r_se", "r_nw", "r_sw"
                                         ),
                                         default_aes = aes(fill = 1, colour = 1, alpha = 1, scale_radii = 1),
                                         draw_key = draw_key_polygon,
                                         draw_group = function(data, panel_scales, coord) {
                                           
                                           ## Transform the data first
                                           coords <- coord$transform(data, panel_scales)
                                           
                                           # Convert nautical miles to meters and multiply by scale factor
                                           data <- data %>% mutate_(r_ne = ~r_ne*1852*scale_radii,
                                                                    r_se = ~r_se*1852*scale_radii,
                                                                    r_sw = ~r_sw*1852*scale_radii,
                                                                    r_nw = ~r_nw*1852*scale_radii
                                           )
                                           
                                           
                                           # Loop over the data and create the points for each quandrant
                                           for (i in 1:nrow(data)) {
                                             
                                             # Create the Northwest Quandrant
                                             df_nw <- base::data.frame(colour = data[i,]$colour,
                                                                       fill = data[i,]$fill,
                                                                       geosphere::destPoint(p = c(data[i,]$x, data[i,]$y),
                                                                                            b = 270:360,
                                                                                            d = data[i,]$r_nw),
                                                                       group = data[i,]$group,
                                                                       PANEL = data[i,]$PANEL,
                                                                       alpha = data[i,]$alpha
                                             )
                                             
                                             # Create the Northeast Quandrant
                                             df_ne <- base::data.frame(colour = data[i,]$colour,
                                                                       fill = data[i,]$fill,
                                                                       geosphere::destPoint(p = c(data[i,]$x, data[i,]$y),
                                                                                            b = 1:90,
                                                                                            d = data[i,]$r_ne),
                                                                       group = data[i,]$group,
                                                                       PANEL = data[i,]$PANEL,
                                                                       alpha = data[i,]$alpha
                                             )
                                             
                                             # Create the Southeast Quandrant
                                             df_se <- base::data.frame(colour = data[i,]$colour,
                                                                       fill = data[i,]$fill,
                                                                       geosphere::destPoint(p = c(data[i,]$x, data[i,]$y),
                                                                                            b = 90:180,
                                                                                            d = data[i,]$r_se),
                                                                       group = data[i,]$group,
                                                                       PANEL = data[i,]$PANEL,
                                                                       alpha = data[i,]$alpha
                                             )
                                             
                                             # Create the Southwest Quandrant
                                             df_sw <- data.frame(colour = data[i,]$colour,
                                                                 fill = data[i,]$fill,
                                                                 geosphere::destPoint(p = c(data[i,]$x, data[i,]$y),
                                                                                      b = 180:270,
                                                                                      d = data[i,]$r_sw),
                                                                 group = data[i,]$group,
                                                                 PANEL = data[i,]$PANEL,
                                                                 alpha = data[i,]$alpha
                                             )
                                             
                                             # bind all the rows into a dataframe
                                             df_points <- dplyr::bind_rows(list(df_nw, df_ne, df_se, df_sw))
                                             
                                           }
                                           
                                           
                                           # Rename columns x and y from lon and lat repectively
                                           df_points <- df_points %>% dplyr::rename_('x' = 'lon',
                                                                                     'y' = 'lat'
                                           )
                                           
                                           # Convert to character
                                           df_points$colour <- base::as.character(df_points$colour)
                                           df_points$fill <- base::as.character(df_points$fill)
                                           
                                           
                                           ## transform data points
                                           coords_df <- coord$transform(df_points, panel_scales)
                                           
                                           ## Construct grid polygon
                                           grid::polygonGrob(
                                             x= coords_df$x,
                                             y = coords_df$y,
                                             gp = grid::gpar(col = coords_df$colour, fill = coords_df$fill, alpha = coords_df$alpha)
                                           )
                                           
                                         }
                                         
)

