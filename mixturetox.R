library(ggplot2)
library(dplyr)
library(tidyr) #seperate
library(scales)
library(readxl)
library(gridExtra)

# df <- read.csv("mixture.csv")

sheets <- excel_sheets("Mixture_Neuron.xlsx")

df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
RName<-as.matrix(df[,2])

df <- data.frame(df[,c(6:13)])
row.names(df) <- RName
names(df) <- c("AC50 min","AC50 Max","Expo min","Expo max","POD min","POD max","RFD min","RFD max")

DF <- df %>% as.matrix() %>%
  reshape::melt() %>% 
  magrittr::set_colnames(c("chemical", "response", "dose"))

df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])

df$pct.MW.AC50.L <- df$`Molecular weight` * df$`Min. AC50`/sum(df$`Min. AC50`)
df$pct.MW.AC50.H <- df$`Molecular weight` * df$`Max. AC50`/sum(df$`Max. AC50`)
df$pct.MW.Css.L <- df$`Molecular weight` * df$Css.medExpos_medRTK.plasma.uM/sum(df$Css.medExpos_medRTK.plasma.uM)
df$pct.MW.Css.H <- df$`Molecular weight` * df$Css.95percExpos_95RTK.plasma.uM/sum(df$Css.95percExpos_95RTK.plasma.uM)
df$pct.MW.RFD.L <- df$`Molecular weight` * df$RFD.Low/sum(df$RFD.Low)
df$pct.MW.RFD.H <- df$`Molecular weight` * df$RFD.High/sum(df$RFD.High)
df$pct.MW.POD.L <- df$`Molecular weight` * df$POD.Lowest/sum(df$POD.Lowest)
df$pct.MW.POD.H <- df$`Molecular weight` * df$POD.Highest/sum(df$POD.Highest)

df$ugl.AC50.L <- df$`Molecular weight` * df$`Min. AC50`
df$ugl.AC50.H <- df$`Molecular weight` * df$`Max. AC50`
df$ugl.Css.L <- df$`Molecular weight` * df$Css.medExpos_medRTK.plasma.uM
df$ugl.Css.H <- df$`Molecular weight` * df$Css.95percExpos_95RTK.plasma.uM
df$ugl.RFD.L <- df$`Molecular weight` * df$RFD.Low
df$ugl.RFD.H <- df$`Molecular weight` * df$RFD.High
df$ugl.POD.L <- df$`Molecular weight` * df$POD.Lowest
df$ugl.POD.H <- df$`Molecular weight` * df$POD.Highest

#
mix.MW.AC50.L <- sum(df$pct.MW.AC50.L)
mix.ugl.AC50.L <- mean(df$ugl.AC50.L)
mixuM.AC50.L <- mix.ugl.AC50.L / mix.MW.AC50.L

mix.MW.AC50.H <- sum(df$pct.MW.AC50.H)
mix.ugl.AC50.H <- mean(df$ugl.AC50.H)
mixuM.AC50.H <- mix.ugl.AC50.H / mix.MW.AC50.H

mix.MW.Css.L <- sum(df$pct.MW.Css.L)
mix.ugl.Css.L <- mean(df$ugl.Css.L)
mixuM.Css.L <- mix.ugl.Css.L / mix.MW.Css.L

mix.MW.Css.H <- sum(df$pct.MW.Css.H)
mix.ugl.Css.H <- mean(df$ugl.Css.H)
mixuM.Css.H <- mix.ugl.Css.H / mix.MW.Css.H

mix.MW.RFD.L <- sum(df$pct.MW.RFD.L)
mix.ugl.RFD.L <- mean(df$ugl.RFD.L)
mixuM.RFD.L <- mix.ugl.RFD.L / mix.MW.RFD.L

mix.MW.RFD.H <- sum(df$pct.MW.RFD.H)
mix.ugl.RFD.H <- mean(df$ugl.RFD.H)
mixuM.RFD.H <- mix.ugl.RFD.H / mix.MW.RFD.H

mix.MW.POD.L <- sum(df$pct.MW.POD.L)
mix.ugl.POD.L <- mean(df$ugl.POD.L)
mixuM.POD.L <- mix.ugl.POD.L / mix.MW.POD.L

mix.MW.POD.H <- sum(df$pct.MW.POD.H)
mix.ugl.POD.H <- mean(df$ugl.POD.H)
mixuM.POD.H <- mix.ugl.POD.H / mix.MW.POD.H

############

png(file="index-DR.png",width=6600,height=4800,res=300)
#pdf("index-DR.pdf", 22, 16)
ggplot(DF, aes(x = reorder(chemical, dose), y=dose, label = response)) +
  geom_path(linetype = 1, color = "grey40") +
  geom_label(size = 4, aes(fill = response)) +
  xlab(expression('chemical'))+
  ylab("reponse conc.") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x, n= 4),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme(legend.position = "none") +
  coord_flip()
dev.off()

DR <- function(i){
  sheets <- readxl::excel_sheets("Mixture_Neuron.xlsx")

  df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[i])
  df <- as.data.frame(df)
  
  df[,1] <- c("AC50 mn", "POD mn", "RfD mx", "Expo mx", "Expo mn", "AC50 mx", "POD mx", "RfD mn")
  colnames(df)<-c("chemical", 
                  "1_r1","1_r2","1_r3","1_r4","1_r5","1_r6",
                  "2_r1","2_r2","2_r3","2_r4","2_r5","2_r6",
                  "3_r1","3_r2","3_r3","3_r4","3_r5","3_r6", 
                  "4_r1","4_r2","4_r3","4_r4","4_r5","4_r6",
                  "5_r1","5_r2","5_r3","5_r4","5_r5","5_r6")
  DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
  DF$dosen <- as.numeric(DF$dose)
  DF1 <- DF %>% mutate(dilution = 10^(dosen-5)) %>% mutate(conc = 10^(dosen-5) * mixuM.POD.H)
  colnames(DF1)[4]<-"response"
  
  ggplot(DF1, aes(x = dilution, y = response)) +
#   ggplot(DF1, aes(x = conc, y = response)) +
    geom_point(aes(colour = round))+
    geom_path(aes(colour = round), size = 0.1) + 
    facet_wrap( ~ chemical, ncol = 4) + theme_bw() +
    ggtitle(sheets[i])+ theme(legend.position = "none") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
}

pdf("mix-DR.pdf", 11, 6)
DR(2);DR(3);DR(4);DR(5);DR(6);DR(7);DR(8);DR(9);DR(10);DR(11) 
dev.off()

png(file="mix-DR.png",width=6800,height=3200,res=300)
grid.arrange(DR(2), DR(3), DR(4), DR(5), 
             DR(6), DR(7), DR(8), DR(9),
             DR(10), DR(11), ncol=5)
dev.off()

png(file="mix-DR.png",width=3200,height=1800,res=300)
DR(3)
dev.off()

#############

df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[2])
df <- as.data.frame(df)
df[,1] <- c("AC50 mn", "POD mn", "RfD mx", "Expo mx", "Expo mn", "AC50 mx", "POD mx", "RfD mn")
colnames(df)<-c("chemical", 
                "1_r1","1_r2","1_r3","1_r4","1_r5","1_r6",
                "2_r1","2_r2","2_r3","2_r4","2_r5","2_r6",
                "3_r1","3_r2","3_r3","3_r4","3_r5","3_r6", 
                "4_r1","4_r2","4_r3","4_r4","4_r5","4_r6",
                "5_r1","5_r2","5_r3","5_r4","5_r5","5_r6")
DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
DF$effect <- sheets[2]
DF$dosen <- as.numeric(DF$dose)
colnames(DF)[4]<-"response"
DF1 <- DF %>% 
  mutate(conc = ifelse(chemical == "AC50 mn", 10^(dosen-5) * mixuM.AC50.L,
                       ifelse(chemical == "POD mn", 10^(dosen-5) * mixuM.POD.L,
                              ifelse(chemical == "RFD mx", 10^(dosen-5) * mixuM.RFD.H,
                                     ifelse(chemical == "Expo mx", 10^(dosen-5) * mixuM.Css.H,
                                            ifelse(chemical == "Expo mn", 10^(dosen-5) * mixuM.Css.L,
                                                   ifelse(chemical == "AC50 mx", 10^(dosen-5) * mixuM.AC50.H,
                                                          ifelse(chemical == "POD mx", 10^(dosen-5) * mixuM.POD.H,
                                                                 10^(dosen-5) * mixuM.RFD.L))))))))
for (i in 3:11) {
  df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[i])
  df <- as.data.frame(df)
  
  df[,1] <- c("AC50 mn", "POD mn", "RfD mx", "Expo mx", "Expo mn", "AC50 mx", "POD mx", "RfD mn")
  colnames(df)<-c("chemical", 
                  "1_r1","1_r2","1_r3","1_r4","1_r5","1_r6",
                  "2_r1","2_r2","2_r3","2_r4","2_r5","2_r6",
                  "3_r1","3_r2","3_r3","3_r4","3_r5","3_r6", 
                  "4_r1","4_r2","4_r3","4_r4","4_r5","4_r6",
                  "5_r1","5_r2","5_r3","5_r4","5_r5","5_r6")
  DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
  DF$effect <- sheets[i]
  DF$dosen <- as.numeric(DF$dose)
  colnames(DF)[4]<-"response"
  DF2 <- DF %>% 
    mutate(conc = ifelse(chemical == "AC50 mn", 10^(dosen-5) * mixuM.AC50.L,
                         ifelse(chemical == "POD mn", 10^(dosen-5) * mixuM.POD.L,
                                ifelse(chemical == "RFD mx", 10^(dosen-5) * mixuM.RFD.H,
                                       ifelse(chemical == "Expo mx", 10^(dosen-5) * mixuM.Css.H,
                                              ifelse(chemical == "Expo mn", 10^(dosen-5) * mixuM.Css.L,
                                                     ifelse(chemical == "AC50 mx", 10^(dosen-5) * mixuM.AC50.H,
                                                            ifelse(chemical == "POD mx", 10^(dosen-5) * mixuM.POD.H,
                                                                   10^(dosen-5) * mixuM.RFD.L))))))))
  DF1 <- rbind(DF1, DF2) 
}

p <- ggplot(DF1, aes(x = conc, y = response)) +
  geom_point(aes(colour = round))+
  #ggtitle(sheets[i]) +
  theme(legend.position = "none") +
  #geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ effect, ncol = 4) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

png(file="mix-mixDR.png",width=4800,height=2800,res=300)
p
dev.off()
