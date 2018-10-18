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
names(df) <- c("AC50 min","AC50 Max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")

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
  df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[i])
  df <- as.data.frame(df)
  
  colnames(df)<-c("chemical", 
                  "1_r1","1_r2","1_r3","1_r4","1_r5","1_r6",
                  "2_r1","2_r2","2_r3","2_r4","2_r5","2_r6",
                  "3_r1","3_r2","3_r3","3_r4","3_r5","3_r6", 
                  "4_r1","4_r2","4_r3","4_r4","4_r5","4_r6",
                  "5_r1","5_r2","5_r3","5_r4","5_r5","5_r6")
  DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
  DF$dosen <- as.numeric(DF$dose)
  DF1 <- DF %>% mutate(dilution = 10^(dosen-5))
  colnames(DF1)[4]<-"response"
  
  ggplot(DF1, aes(x = dilution, y = response)) +
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

png(file="mix-DR.png",width=6400,height=3600,res=300)
grid.arrange(DR(2), DR(3), DR(4), DR(5), 
             DR(6), DR(7), DR(8), DR(9),
             DR(10), DR(11), ncol=5)
dev.off()