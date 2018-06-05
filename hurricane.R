library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggmap)


df<-read_hdata("ebtrk_atlc_1988_2015.txt")

katrina <- df %>% tidy_hdata() %>% filter(storm_id == "Katrina-2005")

katrina <- df %>% tidy_hdata() %>% filter_hdata(h = "Katrina-2005", observation = "2005-08-28 18:00:00")


ggplot(data = katrina) +
  geom_hurricane(aes(x = longitude, y = latitude,
                     r_ne = ne, r_se = se, r_nw = nw, r_sw = sw,
                     fill = wind_speed, color = wind_speed)) +
  scale_color_manual(name = "Wind speed (kts)",
                     values = c("red", "orange", "yellow")) +
  scale_fill_manual(name = "Wind speed (kts)",
                    values = c("red", "orange", "yellow")) 

#
#map_data <- get_map("Louisiana", zoom = 6, maptype = "toner-background")
#base_map <- ggmap(map_data, extent = "device")

base_map +
  geom_hurricane(data = katrina, aes(x = longitude, y = latitude,
                                     r_ne = ne, r_se = se,
                                     r_nw = nw, r_sw = sw,
                                     fill = wind_speed,
                                     color = wind_speed)) +
  scale_color_manual(name = "Wind speed (kts)",
                     values = c("red", "orange", "yellow")) +
  scale_fill_manual(name = "Wind speed (kts)",
                    values = alpha(c("red", "orange", "yellow"), 0.5))