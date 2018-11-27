# Package
library(maps)
library(mapproj)
library(ggplot2)
library(gridExtra)

# Data
states <- map_data("state")
county <-  map_data("county")
tx_county <- subset(county, region == "texas")
brazos <- subset(tx_county, subregion == "brazos")

# Plot
p_US<-ggplot(states, aes(long, lat, group = group)) +
  geom_polygon(fill = I("grey85")) +
  geom_polygon(data = tx_county, fill = NA, color = "white") +
  geom_polygon(data = brazos, fill = "maroon") +
  geom_path(color = "gray") +
  coord_map(project="globular") +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw()

p_TX<-ggplot(tx_county, aes(long, lat, group = group)) +
  geom_polygon(fill = I("grey85")) +
  geom_polygon(data = brazos, fill = "maroon") +
  geom_path(color = "gray") +
  coord_map(project="globular") +
  labs(x=NULL, y=NULL) +
  theme_bw()

p_BCS <- ggplot(data = brazos, mapping = aes(x = long, y = lat)) + 
  geom_polygon(fill = I("maroon")) + 
  geom_text(aes(label = "TAMU", x = -96.336, y = 30.618), col = "white") +
  coord_map(project="globular") +
  labs(x=NULL, y=NULL) +
  theme_bw()
