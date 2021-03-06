if(!require(tidyverse)) {
  install.packages("tidyverse"); require(tidyverse)
}

###############

df<-read.csv("EPA_NATA_2014_HCG.csv")
names(df)

df1 <- df %>% group_by(TRACT, COUNTY, FIPS, population) %>% 
  summarise_at(names(df)[6:40], sum, na.rm = TRUE)

write.csv(df1, file = "TRACT_flood_sum.csv", row.names=FALSE)

################

df<-read.csv("Calculatedarea_HarrisCounty_Flood_Intersect_CSTlevel.csv")
df$GEOID <- as.factor(df$GEOID)

df1 <- df %>% group_by(GEOID, FloodYear) %>% 
  summarise(Sum_Shape_Area = sum(Shape_Area)) %>% 
  as.data.frame() %>%
  mutate(Sum_Shape_Area = replace_na(Sum_Shape_Area, 0)) %>%
  spread(key = FloodYear, value = Sum_Shape_Area) %>%
  replace_na(list(`0` = 0, `100` = 0, `500` = 0)) %>%
  mutate(Total_area = `0`+`100`+`500`) %>%
  mutate(Total_area = replace_na(Total_area, 0)) %>%
  mutate(pct_000_area = `0`/Total_area) %>%
  mutate(pct_100_area = `100`/Total_area) %>%
  mutate(pct_500_area = `500`/Total_area)

df1$`<NA>` <- NULL

write.csv(df1, file = "flood_average_area_(sum).csv", row.names=FALSE)

