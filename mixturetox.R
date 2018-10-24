library(ggplot2)
library(dplyr)
library(tidyr) #seperate
library(scales)
library(readxl)

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
  c.df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
  
  c.df$pct.MW.AC50.L <- c.df$`Molecular weight` * c.df$`Min. AC50`/sum(c.df$`Min. AC50`)
  c.df$pct.MW.AC50.H <- c.df$`Molecular weight` * c.df$`Max. AC50`/sum(c.df$`Max. AC50`)
  c.df$pct.MW.Css.L <- c.df$`Molecular weight` * c.df$Css.medExpos_medRTK.plasma.uM/sum(c.df$Css.medExpos_medRTK.plasma.uM)
  c.df$pct.MW.Css.H <- c.df$`Molecular weight` * c.df$Css.95percExpos_95RTK.plasma.uM/sum(c.df$Css.95percExpos_95RTK.plasma.uM)
  c.df$pct.MW.RFD.L <- c.df$`Molecular weight` * c.df$RFD.Low/sum(c.df$RFD.Low)
  c.df$pct.MW.RFD.H <- c.df$`Molecular weight` * c.df$RFD.High/sum(c.df$RFD.High)
  c.df$pct.MW.POD.L <- c.df$`Molecular weight` * c.df$POD.Lowest/sum(c.df$POD.Lowest)
  c.df$pct.MW.POD.H <- c.df$`Molecular weight` * c.df$POD.Highest/sum(c.df$POD.Highest)
  
  c.df$ugl.AC50.L <- c.df$`Molecular weight` * c.df$`Min. AC50`
  c.df$ugl.AC50.H <- c.df$`Molecular weight` * c.df$`Max. AC50`
  c.df$ugl.Css.L <- c.df$`Molecular weight` * c.df$Css.medExpos_medRTK.plasma.uM
  c.df$ugl.Css.H <- c.df$`Molecular weight` * c.df$Css.95percExpos_95RTK.plasma.uM
  c.df$ugl.RFD.L <- c.df$`Molecular weight` * c.df$RFD.Low
  c.df$ugl.RFD.H <- c.df$`Molecular weight` * c.df$RFD.High
  c.df$ugl.POD.L <- c.df$`Molecular weight` * c.df$POD.Lowest
  c.df$ugl.POD.H <- c.df$`Molecular weight` * c.df$POD.Highest

  #
  mix.MW.AC50.L <- sum(c.df$pct.MW.AC50.L)
  mix.ugl.AC50.L <- mean(c.df$ugl.AC50.L)
  mixuM.AC50.L <- mix.ugl.AC50.L / mix.MW.AC50.L
  
  mix.MW.AC50.H <- sum(c.df$pct.MW.AC50.H)
  mix.ugl.AC50.H <- mean(c.df$ugl.AC50.H)
  mixuM.AC50.H <- mix.ugl.AC50.H / mix.MW.AC50.H
  
  mix.MW.Css.L <- sum(c.df$pct.MW.Css.L)
  mix.ugl.Css.L <- mean(c.df$ugl.Css.L)
  mixuM.Css.L <- mix.ugl.Css.L / mix.MW.Css.L
  
  mix.MW.Css.H <- sum(c.df$pct.MW.Css.H)
  mix.ugl.Css.H <- mean(c.df$ugl.Css.H)
  mixuM.Css.H <- mix.ugl.Css.H / mix.MW.Css.H
  
  mix.MW.RFD.L <- sum(c.df$pct.MW.RFD.L)
  mix.ugl.RFD.L <- mean(c.df$ugl.RFD.L)
  mixuM.RFD.L <- mix.ugl.RFD.L / mix.MW.RFD.L
  
  mix.MW.RFD.H <- sum(c.df$pct.MW.RFD.H)
  mix.ugl.RFD.H <- mean(c.df$ugl.RFD.H)
  mixuM.RFD.H <- mix.ugl.RFD.H / mix.MW.RFD.H
  
  mix.MW.POD.L <- sum(c.df$pct.MW.POD.L)
  mix.ugl.POD.L <- mean(c.df$ugl.POD.L)
  mixuM.POD.L <- mix.ugl.POD.L / mix.MW.POD.L
  
  mix.MW.POD.L <- sum(c.df$pct.MW.POD.L)
  mix.ugl.POD.L <- mean(c.df$ugl.POD.L)
  mixuM.POD.L <- mix.ugl.POD.L / mix.MW.POD.L
  
  mix.MW.POD.H <- sum(c.df$pct.MW.POD.H)
  mix.ugl.POD.H <- mean(c.df$ugl.POD.H)
  mixuM.POD.H <- mix.ugl.POD.H / mix.MW.POD.H

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

sheets <- readxl::excel_sheets("Mixture_Neuron.xlsx")
df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
df$ugl <- df$`Molecular weight`* df$POD.Highest
total_weight <- sum(df$ugl)
df$pct <- df$ugl/total_weight
df$pct.MW.POD.H <- df$`Molecular weight` * df$pct

mix.MW <- sum(df$pct.MW.POD.H)
mix.ugl <- mean(df$ugl)
mixuM <- mix.ugl / mix.MW

