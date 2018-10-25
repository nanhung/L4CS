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

df$ugl.AC50.L <- df$`Molecular weight` * df$`Min. AC50`
df$ugl.AC50.H <- df$`Molecular weight` * df$`Max. AC50`
df$ugl.Css.L <- df$`Molecular weight` * df$Css.medExpos_medRTK.plasma.uM
df$ugl.Css.H <- df$`Molecular weight` * df$Css.95percExpos_95RTK.plasma.uM
df$ugl.RFD.L <- df$`Molecular weight` * df$RFD.Low
df$ugl.RFD.H <- df$`Molecular weight` * df$RFD.High
df$ugl.POD.L <- df$`Molecular weight` * df$POD.Lowest
df$ugl.POD.H <- df$`Molecular weight` * df$POD.Highest

df$pct.MW.AC50.L <- df$`Molecular weight` * df$ugl.AC50.L / sum(df$ugl.AC50.L)
df$pct.MW.AC50.H <- df$`Molecular weight` * df$ugl.AC50.H / sum(df$ugl.AC50.H)
df$pct.MW.Css.L <- df$`Molecular weight` * df$ugl.Css.L /sum(df$ugl.Css.L)
df$pct.MW.Css.H <- df$`Molecular weight` * df$ugl.Css.H /sum(df$ugl.Css.H)
df$pct.MW.RFD.L <- df$`Molecular weight` * df$ugl.RFD.L /sum(df$ugl.RFD.L)
df$pct.MW.RFD.H <- df$`Molecular weight` * df$ugl.RFD.H /sum(df$ugl.RFD.H)
df$pct.MW.POD.L <- df$`Molecular weight` * df$ugl.POD.L /sum(df$ugl.POD.L)
df$pct.MW.POD.H <- df$`Molecular weight` * df$ugl.POD.H /sum(df$ugl.POD.H)

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

##################### PredTox

sheets <- readxl::excel_sheets("Mixture_Neuron.xlsx")
df1 <- readxl::read_xlsx("42_Chem_Neuron.xlsx", sheet = sheets[2])
df2 <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
df2[,1] <- df1[,2]

colnames(df1)<-c("chemical", "1_r1","1_r2", "2_r1","2_r2", 
                 "3_r1","3_r2", "4_r1","4_r2", "5_r1","5_r2")
DF <- df1 %>% as.data.frame () %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
DF$dosen <- as.numeric(DF$dose)
DF2 <- DF %>% mutate(dose = 10^(dosen-3))
colnames(DF2)[4]<-"response"

df2 <- data.frame(df2[,c(2,6:13)])
names(df2) <- c("chemical","AC50 min","AC50 max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")

chem <- as.data.frame(df1[,1])

source("mixECx.R")

effv_est <- function(i, init_n = 1){
  X <- Fit(i, init_n)
  C <- df2[i,'POD max']
  y <- X$p[3] / (1 + (C/X$p[1])^X$p[2])
  return(y)  
}


effv <- df2$`POD max` / sum(df2$`POD max`)


effPoints <- rev((c(0.025, 0.03, 0.05, 0.1, 0.15, 0.2, 
                    0.25, 0.3, 0.35, 0.4, 0.45, 0.47, 0.5, 0.52, 
                    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.999)))
pctEcx <- t(t(effv/sum(effv)))


model <- rep('Hill_three_rev', 42)
param <- matrix(c(c(df2$`AC50 max`[1], 1, 1),
                  c(df2$`AC50 max`[2], 1, 1),
                  c(df2$`AC50 max`[3], 1, 1),
                  Fit(4, effect = sheets[3])$p,
                  Fit(5, effect = sheets[3])$p,
                  Fit(6, effect = sheets[3])$p,
                  Fit(7, effect = sheets[3])$p,
                  c(df2$`AC50 max`[8], 1, 1),
                  c(df2$`AC50 max`[9], 1, 1),
                  c(df2$`AC50 max`[10], 1, 1),
                  c(df2$`AC50 max`[11], 1, 1),
                  Fit(12, effect = sheets[3])$p,
                  c(df2$`AC50 max`[13], 1, 1),
                  c(df2$`AC50 max`[14], 1, 1),
                  Fit(15, effect = sheets[3])$p,
                  Fit(16, effect = sheets[3])$p,
                  c(df2$`AC50 max`[17], 1, 1),
                  c(df2$`AC50 max`[18], 1, 1),
                  Fit(18, effect = sheets[3])$p,
                  Fit(19, effect = sheets[3])$p,
                  Fit(20, effect = sheets[3])$p,
                  c(df2$`AC50 max`[21], 1, 1),
                  Fit(22, effect = sheets[3])$p,
                  Fit(23, effect = sheets[3])$p,
                  c(df2$`AC50 max`[24], 1, 1),
                  Fit(25, effect = sheets[3])$p,
                  c(df2$`AC50 max`[26], 1, 1),
                  c(43.922406, 14.189828,  1.088035),
                  Fit(28, effect = sheets[3])$p,
                  c(df2$`AC50 max`[29], 1, 1),
                  c(df2$`AC50 max`[30], 1, 1),
                  c(43.907095, 14.421544,  1.066666),
                  c(df2$`AC50 max`[32], 1, 1),
                  Fit(33, effect = sheets[3])$p,
                  c(df2$`AC50 max`[34], 1, 1),
                  c(df2$`AC50 max`[35], 1, 1),
                  Fit(36, effect = sheets[3])$p,
                  Fit(37, effect = sheets[3])$p,
                  c(df2$`AC50 max`[38], 1, 1),
                  Fit(39, effect = sheets[3])$p,
                  Fit(40, effect = sheets[3])$p,
                  Fit(41, effect = sheets[3])$p,
                  Fit(42, effect = sheets[3])$p), byrow = T, ncol =3)

ECx(model, param, effPoints)

concAdd <- function(pctEcx, effPoints) {
  ecPoints <- ECx(model, param, effPoints)
  ca <- 1/(t(pctEcx) %*% (1/ecPoints))
  return(ca)
}

ca <- concAdd(pctEcx, rev(effPoints))

ca <- c(10e-5, as.numeric(ca))
effPoints<- c(0.999, effPoints)

x<-filter(DF1, effect == "cellnum")["conc"] %>% as.matrix()
y<-filter(DF1, effect == "cellnum")["response"]/100


png(file="mix-ca.png",width=3600,height=2800,res=300)
plot(x,as.matrix(y), log = "x", pch = 19, 
     main = "Concentration addition",
     xlab = "log(c) mol/L", ylab = "Inhibition (%)")
lines(ca, effPoints, lwd = 2, col = 2)
dev.off()

##########################

effv_POD_min <- df2$`POD min` / sum(df2$`POD min`)
effv_POD_max <- df2$`POD max` / sum(df2$`POD max`)
effv_AC50_min <- df2$`AC50 min` / sum(df2$`AC50 min`)
effv_AC50_max <- df2$`AC50 max` / sum(df2$`AC50 max`)
effv_Expo_max <- df2$`Expo max` / sum(df2$`Expo max`)
effv_Expo_min <- df2$`Expo min` / sum(df2$`Expo min`)
effv_RFD_min <- df2$`RFD min` / sum(df2$`RFD min`)
effv_RFD_max <- df2$`RFD max` / sum(df2$`RFD max`)

X <- data.frame(effv_AC50_min, effv_AC50_max,
                effv_Expo_min, effv_Expo_max, effv_POD_min, effv_POD_max, effv_RFD_min, effv_RFD_max)

row.names(X) <- df2[,1]
colnames(X) <- c("AC50 min","AC50 max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")

pct_df <- X %>% as.matrix() %>% reshape2::melt()
names(pct_df) <- c("chemical", "EC", "percentage")
pct_df <- ddply(pct_df, .(EC), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

png(file="mixtox-3.png",width=4800,height=2800,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = EC, fill = chemical), data = pct_df, stat="identity")+
  ggtitle("Percentage of individual chemicals in the mixtures")+
  #geom_text(data=pct_df, aes(x = EC, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+ scale_fill_viridis(discrete=TRUE) +
  ylab("Percentage (%)")
dev.off()
