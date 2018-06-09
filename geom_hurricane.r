read_EBT_data <- function(filename) {
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


geom_hurricane <- function(mapping = NULL, data = NULL, stat = 'identity',
                           position = 'identity', 
                           show.legend = NA, ...) {
  ggplot2::layer(
    geom = geom_hurricane_proto, mapping = mapping,
    data = data, stat = stat, position = position, 
    show.legend = show.legend, inherit.aes = TRUE,
    params = list(na.rm = FALSE, ...))
}

geom_hurricane_proto <- ggplot2::ggproto("geom_hurricane_proto", Geom,
                                         required_aes = c("x", "y","r_ne", "r_se", "r_nw", "r_sw"),
                                         default_aes = aes(fill=1, 
                                                           colour=1, 
                                                           size=0.5,
                                                           linetype=1,
                                                           alpha=.5,
                                                           arc_step=1,
                                                           scale_radii=1),
                                         draw_key = draw_key_polygon,
                                         draw_group = function(data, panel_scales, coord) {
                                           coords <- coord$transform(data, panel_scales)
                                           
                                           data <- data %>% mutate_(r_ne = ~r_ne*1852*scale_radii,
                                                                    r_se = ~r_se*1852*scale_radii,
                                                                    r_sw = ~r_sw*1852*scale_radii,
                                                                    r_nw = ~r_nw*1852*scale_radii)
                                           
                                           
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

