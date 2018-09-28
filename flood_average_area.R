if(!require(dplyr)) {
  install.packages("dplyr"); require(dplyr)
}

df<-read.csv("Calculatedarea_HarrisCounty_Flood_Intersect_CSTlevel.csv")
df$NAME <- as.factor(df$NAME)

df1 <- df %>% group_by(NAME, FloodYear) %>% 
  summarise(Sum_Shape_Area = sum(Shape_Area))

df1$average_area <- df1$FloodYear / df1$Sum_Shape_Area 

unique(df$NAME)
unique(df$FloodYear)

x<-rep(NA, nrow(df))

for(i in 1:nrow(df)){
  for(j in unique(df$NAME)){
    for(k in unique(df$FloodYear)){
      if(df$NAME[i] == j & df$FloodYear[i] == k){
        x[i] <- as.numeric(df1[ which(df1$NAME==j & df1$FloodYear == k), "average_area"])
      }
    }
  }
}

df$average_area <- x

write.csv(df, file = "flood_average_area.csv", row.names=FALSE)
write.csv(df1, file = "flood_average_area_(sum).csv", row.names=FALSE)

