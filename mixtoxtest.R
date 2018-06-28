# devtools::install_github("lahothorn/SiTuR")
library(tukeytrend)
library(mixtox)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr) #seperate

df <- read.csv("mixtoxtest.csv")
colnames(df)<-c("chemical", "1_r1","1_r2", "2_r1","2_r2", 
            "3_r1","3_r2", "4_r1","4_r2", "5_r1","5_r2")
DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
DF$dosen <- as.numeric(DF$dose)
DF1 <- DF %>% mutate(dose = 10^(dosen-3))
colnames(DF1)[4]<-"response"

ggplot(DF1, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 7) + 
  scale_x_log10() + theme_bw()

df1 <- filter(DF1, chemical == "ALDRIN")
fitw <- lm(response~dose, data = df1)
ttw<-tukeytrendfit(fitw, dose="dose", scaling = c("ari", "ord", "log"))
s<-summary(asglht(ttw))
s$test$tstat
s$test$pvalues[1]
  
df2 <- filter(DF1, chemical == "BENZIDINE")
fitw <- lm(response~dose, data = df2)
ttw<-tukeytrendfit(fitw, dose="dose", scaling = c("ari", "ord", "log"))
summary(asglht(ttw))

df3 <- filter(DF1, chemical == "LEAD (Nitrate)")
fitw <- lm(response~dose, data = df3)
ttw<-tukeytrendfit(fitw, dose="dose", scaling = c("ari", "ord", "log"))
summary(asglht(ttw))


