#install.packages(c("ggplot2","dplyr", "tidyr", "scales", "readxl", "gridExtra"))

library(plyr)
library(ggplot2)
library(dplyr)
library(tidyr) #seperate
library(scales)
library(readxl)
library(gridExtra)
library(reshape2)

# df <- read.csv("mixture.csv")

sheets <- excel_sheets("Mixture_Neuron.xlsx")

df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])
RName<-as.matrix(df[,2])

df <- data.frame(df[,c(6:13)])
row.names(df) <- RName
names(df) <- c("AC50 min","AC50 Max","Expo min","Expo max","POD min","POD max","RFD min","RFD max")

X <- df %>% as.matrix() %>%
  reshape::melt() %>% 
  magrittr::set_colnames(c("chemical", "response", "dose"))

df <- readxl::read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1])


df$pct.MW.AC50.L <- df$`Min. AC50` / sum(df$`Min. AC50`)
df$pct.MW.AC50.H <- df$`Max. AC50` / sum(df$`Max. AC50`)
df$pct.MW.Css.L <- df$Css.medExpos_medRTK.plasma.uM /sum(df$Css.medExpos_medRTK.plasma.uM)
df$pct.MW.Css.H <- df$Css.95percExpos_95RTK.plasma.uM /sum(df$Css.95percExpos_95RTK.plasma.uM)
df$pct.MW.RFD.L <- df$RFD.Low /sum(df$RFD.Low)
df$pct.MW.RFD.H <- df$RFD.High /sum(df$RFD.High)
df$pct.MW.POD.L <- df$POD.Lowest /sum(df$POD.Lowest)
df$pct.MW.POD.H <- df$POD.Highest /sum(df$POD.Highest)

#
mix.MW.AC50.L <- sum(df$pct.MW.AC50.L)
mix.AC50.L <- mean(df$`Min. AC50`)
mixuM.AC50.L <- mix.AC50.L / mix.MW.AC50.L

mix.MW.AC50.H <- sum(df$pct.MW.AC50.H)
mix.AC50.H <- mean(df$`Max. AC50`)
mixuM.AC50.H <- mix.AC50.H / mix.MW.AC50.H

mix.MW.Css.L <- sum(df$pct.MW.Css.L)
mix.Css.L <- mean(df$Css.medExpos_medRTK.plasma.uM)
mixuM.Css.L <- mix.Css.L / mix.MW.Css.L

mix.MW.Css.H <- sum(df$pct.MW.Css.H)
mix.Css.H <- mean(df$Css.95percExpos_95RTK.plasma.uM)
mixuM.Css.H <- mix.Css.H / mix.MW.Css.H

mix.MW.RFD.L <- sum(df$pct.MW.RFD.L)
mix.RFD.L <- mean(df$RFD.Low)
mixuM.RFD.L <- mix.RFD.L / mix.MW.RFD.L

mix.MW.RFD.H <- sum(df$pct.MW.RFD.H)
mix.RFD.H <- mean(df$RFD.High)
mixuM.RFD.H <- mix.RFD.H / mix.MW.RFD.H

mix.MW.POD.L <- sum(df$pct.MW.POD.L)
mix.POD.L <- mean(df$POD.Lowest)
mixuM.POD.L <- mix.POD.L / mix.MW.POD.L

mix.MW.POD.H <- sum(df$pct.MW.POD.H)
mix.POD.H <- mean(df$POD.Highest)
mixuM.POD.H <- mix.POD.H / mix.MW.POD.H

effv_POD_min <- df$POD.Lowest / sum(df$POD.Lowest)
effv_POD_max <- df$POD.Highest / sum(df$POD.Highest)
effv_AC50_min <- df$`Min. AC50` / sum(df$`Min. AC50`)
effv_AC50_max <- df$`Max. AC50` / sum(df$`Max. AC50`)
effv_Expo_max <- df$Css.medExpos_medRTK.plasma.uM / sum(df$Css.medExpos_medRTK.plasma.uM)
effv_Expo_min <- df$Css.95percExpos_95RTK.plasma.uM / sum(df$Css.95percExpos_95RTK.plasma.uM)
effv_RFD_min <- df$RFD.Low / sum(df$RFD.Low)
effv_RFD_max <- df$RFD.High / sum(df$RFD.High)

X <- data.frame(effv_AC50_min, effv_AC50_max,
                effv_Expo_min, effv_Expo_max, effv_POD_min, effv_POD_max, effv_RFD_min, effv_RFD_max)

row.names(X) <- as.matrix(df[,2])
colnames(X) <- c("AC50 min","AC50 max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")

pct_df <- X %>% as.matrix() %>% reshape2::melt()
names(pct_df) <- c("chemical", "EC", "percentage")
pct_df <- plyr::ddply(pct_df, .(EC), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

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

png(file="mixtox-3.png",width=4800,height=2800,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = EC, fill = chemical), data = pct_df, stat="identity")+
  ggtitle("Percentage of individual chemicals in the mixtures")+
  #geom_text(data=pct_df, aes(x = EC, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+ viridis::scale_fill_viridis(discrete=TRUE) +
  ylab("Percentage (%)")
dev.off()

png(file="index-DR.png",width=6600,height=4800,res=300)
#pdf("index-DR.pdf", 22, 16)
ggplot(X, aes(x = reorder(chemical, dose), y=dose, label = response)) +
  geom_path(linetype = 1, color = "grey40") +
  geom_label(size = 4, aes(fill = response)) +
  xlab('chemical')+
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
  geom_point(aes(colour = chemical))+
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

df3 <- df2
df2 <- data.frame(df2[,c(2,6:13)])
names(df2) <- c("chemical","AC50 min","AC50 max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")

chem <- as.data.frame(df1[,1])

#effv_est <- function(i, init_n = 1){
#  X <- Fit(i, init_n)
#  C <- df2[i,'POD max']
#  y <- X$p[3] / (1 + (C/X$p[1])^X$p[2])
#  return(y)  
#}

source("mixECx.R")

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


effPoints <- rev((c(0.05, 0.1, 0.15, 0.2, 
                    0.25, 0.3, 0.35, 0.4, 0.45, 0.47, 0.5, 0.52, 
                    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.999)))


#eff <- df2$`AC50 max`
#effv <- eff / sum(eff)
#pctEcx <- t(t(effv/sum(effv)))


ca_estimate <- function(data){
  eff <- data
  effv <- eff / sum(eff)
  pctEcx <- t(t(effv/sum(effv)))
  
  concAdd <- function(pctEcx, effPoints) {
    ecPoints <- ECx(model, param, effPoints)
    ca <- 1/(t(pctEcx) %*% (1/ecPoints))
    return(ca)
  }
  
  ca <- concAdd(pctEcx, rev(effPoints))
  ca <- c(1e-4, as.numeric(ca))
  
  return(ca)
}

ia_estimate <- function(data){
  eff <- data
  effv <- eff / sum(eff)
  pctEcx <- t(t(effv/sum(effv)))
  
  ia <- indAct(model, param, pctEcx, effPoints)
  ia <- c(1e-4, as.numeric(ia))
  return(ia)
}

ia_AC5mn <- ia_estimate(df2$`AC50 min`)
ia_AC5mx <- ia_estimate(df2$`AC50 max`)
ia_Expmx <- ia_estimate(df2$`Expo min`)
ia_Expmn <- ia_estimate(df2$`Expo max`)
ia_PODmn <- ia_estimate(df2$`POD min`)
ia_PODmx <- ia_estimate(df2$`POD max`)
ia_RfDmn <- ia_estimate(df2$`RFD min`)
ia_RfDmx <- ia_estimate(df2$`RFD max`)

ca_AC5mx <- ca_estimate(df2$`AC50 max`)
ca_AC5mn <- ca_estimate(df2$`AC50 min`)
ca_Expmx <- ca_estimate(df2$`Expo min`)
ca_Expmn <- ca_estimate(df2$`Expo max`)
ca_PODmn <- ca_estimate(df2$`POD min`)
ca_PODmx <- ca_estimate(df2$`POD max`)
ca_RfDmn <- ca_estimate(df2$`RFD min`)
ca_RfDmx <- ca_estimate(df2$`RFD max`)

EffPoints <- c(0.999, effPoints)
rsp_rng <- range(filter(DF1, effect == "cellnum")["response"]/100)

x_AC5mx <-filter(DF1, effect == "cellnum" & chemical == "AC50 mx")["conc"] %>% as.matrix()
y_AC5mx <-filter(DF1, effect == "cellnum" & chemical == "AC50 mx")["response"]/100
x_AC5mn <-filter(DF1, effect == "cellnum" & chemical == "AC50 mn")["conc"] %>% as.matrix()
y_AC5mn <-filter(DF1, effect == "cellnum" & chemical == "AC50 mn")["response"]/100
x_Expmx <-filter(DF1, effect == "cellnum" & chemical == "Expo mx")["conc"] %>% as.matrix()
y_Expmx <-filter(DF1, effect == "cellnum" & chemical == "Expo mx")["response"]/100
x_Expmn <-filter(DF1, effect == "cellnum" & chemical == "Expo mn")["conc"] %>% as.matrix()
y_Expmn <-filter(DF1, effect == "cellnum" & chemical == "Expo mn")["response"]/100
x_PODmx <-filter(DF1, effect == "cellnum" & chemical == "POD mx")["conc"] %>% as.matrix()
y_PODmx <-filter(DF1, effect == "cellnum" & chemical == "POD mx")["response"]/100
x_PODmn <-filter(DF1, effect == "cellnum" & chemical == "POD mn")["conc"] %>% as.matrix()
y_PODmn <-filter(DF1, effect == "cellnum" & chemical == "POD mn")["response"]/100
x_RfDmx <-filter(DF1, effect == "cellnum" & chemical == "RfD mx")["conc"] %>% as.matrix()
y_RfDmx <-filter(DF1, effect == "cellnum" & chemical == "RfD mx")["response"]/100
x_RfDmn <-filter(DF1, effect == "cellnum" & chemical == "RfD mn")["conc"] %>% as.matrix()
y_RfDmn <-filter(DF1, effect == "cellnum" & chemical == "RfD mn")["response"]/100

png(file="mix-caia.png",width=3600,height=1800,res=300)
par(mfrow = c(2,4))
plot(x_AC5mn,as.matrix(y_AC5mn), log = "x", pch = 19, 
     main = "AC50 min", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_AC5mn, EffPoints, lwd = 2, col = 2)
lines(ia_AC5mn, EffPoints, lwd = 2, col = 3)

plot(x_AC5mx,as.matrix(y_AC5mx), log = "x", pch = 19, 
     main = "AC50 max", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_AC5mx, EffPoints, lwd = 2, col = 2)
lines(ia_AC5mx, EffPoints, lwd = 2, col = 3)

plot(x_Expmn, as.matrix(y_Expmn), log = "x", pch = 19, 
     main = "Expo min", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_Expmn, EffPoints, lwd = 2, col = 2)
lines(ia_Expmn, EffPoints, lwd = 2, col = 3)

plot(x_Expmx, as.matrix(y_Expmx), log = "x", pch = 19, 
     main = "Expo max", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_Expmx, EffPoints, lwd = 2, col = 2)
lines(ia_Expmx, EffPoints, lwd = 2, col = 3)

plot(x_PODmn,as.matrix(y_PODmn), log = "x", pch = 19, 
     main = "POD min", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_PODmn, EffPoints, lwd = 2, col = 2)
lines(ia_PODmn, EffPoints, lwd = 2, col = 3)

plot(x_PODmx,as.matrix(y_PODmx), log = "x", pch = 19, 
     main = "POD max", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_PODmx, EffPoints, lwd = 2, col = 2)
lines(ia_PODmx, EffPoints, lwd = 2, col = 3)

plot(x_RfDmn,as.matrix(y_RfDmn), log = "x", pch = 19, 
     main = "RfD min", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_RfDmn, EffPoints, lwd = 2, col = 2)
lines(ia_RfDmn, EffPoints, lwd = 2, col = 3)

plot(x_RfDmx,as.matrix(y_RfDmx), log = "x", pch = 19, 
     main = "RfD max", ylim = rsp_rng,
     xlab = expression(paste("Concentration (", mu,"M)")), ylab = "Inhibition (%)")
lines(ca_RfDmx, EffPoints, lwd = 2, col = 2)
lines(ia_RfDmx, EffPoints, lwd = 2, col = 3)

dev.off()




