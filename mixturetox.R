library(ggplot2)
library(dplyr)
library(tidyr) #seperate
library(scales)

df <- read.csv("mixture.csv")
colnames(df)<-c("chemical", 
                "1_r1","1_r2","1_r3","1_r4","1_r5","1_r6",
                "2_r1","2_r2","2_r3","2_r4","2_r5","2_r6",
                "3_r1","3_r2","3_r3","3_r4","3_r5","3_r6", 
                "4_r1","4_r2","4_r3","4_r4","4_r5","4_r6",
                "5_r1","5_r2","5_r3","5_r4","5_r5","5_r6")
DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
DF$dosen <- as.numeric(DF$dose)
DF1 <- DF %>% mutate(dose = 10^(dosen-5))
colnames(DF1)[4]<-"response"

ggplot(DF1, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 4) + theme_bw() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))