## Load package and data ####
library(tidyverse)
library(readxl)
library(reshape)
library(scales)
source("mixECx.R") # defined function (revised from mixtox package)

sheets <- excel_sheets("Mixture_Neuron.xlsx")
x1 <- read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[1]) %>% as.data.frame() # mixture info
x2 <- read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[2]) %>% as.data.frame() # mixture dose-response
X <- read_xlsx("42_Chem_Neuron.xlsx", sheet = sheets[2]) # Use cell number change as response

## Estimate mixture concentration ####
mixuM.AC50.L <-  mean(x1$`Min. AC50`)
mixuM.AC50.H <-  mean(x1$`Max. AC50`)
mixuM.POD.L <-  mean(x1$POD.Highest)
mixuM.POD.H <-  mean(x1$POD.Highest)
mixuM.Css.L <-  mean(x1$Css.medExpos_medRTK.plasma.uM)
mixuM.Css.H <-  mean(x1$Css.95percExpos_95RTK.plasma.uM)
mixuM.RfD.L <-  mean(x1$RFD.Low)
mixuM.RfD.H <-  mean(x1$RFD.High)

## Generate Tidy Data
x2[,1] <- c("AC50 mn", "POD mn", "RfD mx", "Expo mx", "Expo mn", "AC50 mx", "POD mx", "RfD mn")
colnames(x2)<-c("chemical", 
                "1_r1","1_r2","1_r3","1_r4","1_r5","1_r6",
                "2_r1","2_r2","2_r3","2_r4","2_r5","2_r6",
                "3_r1","3_r2","3_r3","3_r4","3_r5","3_r6", 
                "4_r1","4_r2","4_r3","4_r4","4_r5","4_r6",
                "5_r1","5_r2","5_r3","5_r4","5_r5","5_r6")
X2 <- x2 %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
X2$effect <- sheets[2]
X2$dosen <- as.numeric(X2$dose)
colnames(X2)[4]<-"response"

DF <- X2 %>% 
  mutate(conc = ifelse(chemical == "AC50 mn", 10^(dosen-5) * mixuM.AC50.L,
                       ifelse(chemical == "POD mn", 10^(dosen-5) * mixuM.POD.L,
                              ifelse(chemical == "RFD mx", 10^(dosen-5) * mixuM.RfD.H,
                                     ifelse(chemical == "Expo mx", 10^(dosen-5) * mixuM.Css.H,
                                            ifelse(chemical == "Expo mn", 10^(dosen-5) * mixuM.Css.L,
                                                   ifelse(chemical == "AC50 mx", 10^(dosen-5) * mixuM.AC50.H,
                                                          ifelse(chemical == "POD mx", 10^(dosen-5) * mixuM.POD.H,
                                                                 10^(dosen-5) * mixuM.RfD.L))))))))
for (i in 3:11) {
  x2 <- read_xlsx("Mixture_Neuron.xlsx", sheet = sheets[i]) %>% as.data.frame() 
  
  x2[,1] <- c("AC50 mn", "POD mn", "RfD mx", "Expo mx", "Expo mn", "AC50 mx", "POD mx", "RfD mn")
  colnames(x2)<-c("chemical", 
                  "1_r1","1_r2","1_r3","1_r4","1_r5","1_r6",
                  "2_r1","2_r2","2_r3","2_r4","2_r5","2_r6",
                  "3_r1","3_r2","3_r3","3_r4","3_r5","3_r6", 
                  "4_r1","4_r2","4_r3","4_r4","4_r5","4_r6",
                  "5_r1","5_r2","5_r3","5_r4","5_r5","5_r6")
  X2 <- x2 %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
  X2$effect <- sheets[i]
  X2$dosen <- as.numeric(X2$dose)
  colnames(X2)[4]<-"response"
  DF2 <- X2 %>% 
    mutate(conc = ifelse(chemical == "AC50 mn", 10^(dosen-5) * mixuM.AC50.L,
                         ifelse(chemical == "POD mn", 10^(dosen-5) * mixuM.POD.L,
                                ifelse(chemical == "RFD mx", 10^(dosen-5) * mixuM.RfD.H,
                                       ifelse(chemical == "Expo mx", 10^(dosen-5) * mixuM.Css.H,
                                              ifelse(chemical == "Expo mn", 10^(dosen-5) * mixuM.Css.L,
                                                     ifelse(chemical == "AC50 mx", 10^(dosen-5) * mixuM.AC50.H,
                                                            ifelse(chemical == "POD mx", 10^(dosen-5) * mixuM.POD.H,
                                                                   10^(dosen-5) * mixuM.RfD.L))))))))
  DF <- rbind(DF, DF2) 
}

# plot
ggplot(DF, aes(x = conc, y = response)) +
  geom_point(aes(colour = chemical))+
  theme(legend.position = "none") +
  facet_wrap( ~ effect, ncol = 5) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

## Curve Fitting ####
colnames(X)<-c("chemical", "1_r1","1_r2", "2_r1","2_r2", "3_r1","3_r2", "4_r1","4_r2", "5_r1","5_r2")
XX <- X %>% as.data.frame () %>% reshape::melt() %>% separate(variable, c("dose", "round")) %>% 
  mutate (Dose = 10^(as.numeric(dose)-3))
colnames(XX)[4]<-"response"
X1 <- data.frame(x1[,c(2,6:13)])
names(X1) <- c("chemical","AC50 min","AC50 max","Expo max","Expo min","POD min","POD max","RFD min","RFD max")
chem <- as.data.frame(X[,1])
model <- rep('Hill_three_rev', 42)
param <- matrix(c(c(X1$`AC50 max`[1], 1, 1),
                  c(X1$`AC50 max`[2], 1, 1),
                  c(X1$`AC50 max`[3], 1, 1),
                  Fit(4, effect = sheets[3])$p,
                  Fit(5, effect = sheets[3])$p,
                  Fit(6, effect = sheets[3])$p,
                  Fit(7, effect = sheets[3])$p,
                  c(X1$`AC50 max`[8], 1, 1),
                  c(X1$`AC50 max`[9], 1, 1),
                  c(X1$`AC50 max`[10], 1, 1),
                  c(X1$`AC50 max`[11], 1, 1),
                  Fit(12, effect = sheets[3])$p,
                  c(X1$`AC50 max`[13], 1, 1),
                  c(X1$`AC50 max`[14], 1, 1),
                  Fit(15, effect = sheets[3])$p,
                  Fit(16, effect = sheets[3])$p,
                  c(X1$`AC50 max`[17], 1, 1),
                  c(X1$`AC50 max`[18], 1, 1),
                  Fit(18, effect = sheets[3])$p,
                  Fit(19, effect = sheets[3])$p,
                  Fit(20, effect = sheets[3])$p,
                  c(X1$`AC50 max`[21], 1, 1),
                  Fit(22, effect = sheets[3])$p,
                  Fit(23, effect = sheets[3])$p,
                  c(X1$`AC50 max`[24], 1, 1),
                  Fit(25, effect = sheets[3])$p,
                  c(X1$`AC50 max`[26], 1, 1),
                  c(43.922406, 14.189828,  1.088035),
                  Fit(28, effect = sheets[3])$p,
                  c(X1$`AC50 max`[29], 1, 1),
                  c(X1$`AC50 max`[30], 1, 1),
                  c(43.907095, 14.421544,  1.066666),
                  c(X1$`AC50 max`[32], 1, 1),
                  Fit(33, effect = sheets[3])$p,
                  c(X1$`AC50 max`[34], 1, 1),
                  c(X1$`AC50 max`[35], 1, 1),
                  Fit(36, effect = sheets[3])$p,
                  Fit(37, effect = sheets[3])$p,
                  c(X1$`AC50 max`[38], 1, 1),
                  Fit(39, effect = sheets[3])$p,
                  Fit(40, effect = sheets[3])$p,
                  Fit(41, effect = sheets[3])$p,
                  Fit(42, effect = sheets[3])$p), byrow = T, ncol =3)

effPoints <- rev((c(0.05, 0.1, 0.15, 0.2, 
                    0.25, 0.3, 0.35, 0.4, 0.45, 0.47, 0.5, 0.52, 
                    0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.999)))
EffPoints <- c(0.999, effPoints)

## Mixture toxicity prediction ####
ia_AC5mn <- ia_estimate(X1$`AC50 min`)
ia_AC5mx <- ia_estimate(X1$`AC50 max`)
ia_Expmx <- ia_estimate(X1$`Expo min`)
ia_Expmn <- ia_estimate(X1$`Expo max`)
ia_PODmn <- ia_estimate(X1$`POD min`)
ia_PODmx <- ia_estimate(X1$`POD max`)
ia_RfDmn <- ia_estimate(X1$`RFD min`)
ia_RfDmx <- ia_estimate(X1$`RFD max`)

ca_AC5mx <- ca_estimate(X1$`AC50 max`)
ca_AC5mn <- ca_estimate(X1$`AC50 min`)
ca_Expmx <- ca_estimate(X1$`Expo min`)
ca_Expmn <- ca_estimate(X1$`Expo max`)
ca_PODmn <- ca_estimate(X1$`POD min`)
ca_PODmx <- ca_estimate(X1$`POD max`)
ca_RfDmn <- ca_estimate(X1$`RFD min`)
ca_RfDmx <- ca_estimate(X1$`RFD max`)


x_AC5mx <-filter(DF, effect == "cellnum" & chemical == "AC50 mx")["conc"] %>% as.matrix()
y_AC5mx <-filter(DF, effect == "cellnum" & chemical == "AC50 mx")["response"]/100
x_AC5mn <-filter(DF, effect == "cellnum" & chemical == "AC50 mn")["conc"] %>% as.matrix()
y_AC5mn <-filter(DF, effect == "cellnum" & chemical == "AC50 mn")["response"]/100
x_Expmx <-filter(DF, effect == "cellnum" & chemical == "Expo mx")["conc"] %>% as.matrix()
y_Expmx <-filter(DF, effect == "cellnum" & chemical == "Expo mx")["response"]/100
x_Expmn <-filter(DF, effect == "cellnum" & chemical == "Expo mn")["conc"] %>% as.matrix()
y_Expmn <-filter(DF, effect == "cellnum" & chemical == "Expo mn")["response"]/100
x_PODmx <-filter(DF, effect == "cellnum" & chemical == "POD mx")["conc"] %>% as.matrix()
y_PODmx <-filter(DF, effect == "cellnum" & chemical == "POD mx")["response"]/100
x_PODmn <-filter(DF, effect == "cellnum" & chemical == "POD mn")["conc"] %>% as.matrix()
y_PODmn <-filter(DF, effect == "cellnum" & chemical == "POD mn")["response"]/100
x_RfDmx <-filter(DF, effect == "cellnum" & chemical == "RfD mx")["conc"] %>% as.matrix()
y_RfDmx <-filter(DF, effect == "cellnum" & chemical == "RfD mx")["response"]/100
x_RfDmn <-filter(DF, effect == "cellnum" & chemical == "RfD mn")["conc"] %>% as.matrix()
y_RfDmn <-filter(DF, effect == "cellnum" & chemical == "RfD mn")["response"]/100


# Plot
rsp_rng <- range(filter(DF, effect == "cellnum")["response"]/100)

par(mfrow = c(2,4))
mixtoxPlot(x_AC5mn, y_AC5mn, ca_AC5mn, ia_AC5mn, EffPoints, main = "AC50 min", ylim = rsp_rng)
mixtoxPlot(x_AC5mx, y_AC5mx, ca_AC5mx, ia_AC5mx, EffPoints, main = "AC50 max", ylim = rsp_rng)
mixtoxPlot(x_Expmn, y_Expmn, ca_Expmn, ia_Expmn, EffPoints, main = "Expo min", ylim = rsp_rng)
mixtoxPlot(x_Expmx, y_Expmx, ca_Expmx, ia_Expmx, EffPoints, main = "Expo max", ylim = rsp_rng)
mixtoxPlot(x_PODmn, y_PODmn, ca_PODmn, ia_PODmn, EffPoints, main = "POD min", ylim = rsp_rng)
mixtoxPlot(x_PODmx, y_PODmx, ca_PODmx, ia_PODmx, EffPoints, main = "POD max", ylim = rsp_rng)
mixtoxPlot(x_RfDmn, y_RfDmn, ca_RfDmn, ia_RfDmn, EffPoints, main = "RfD min", ylim = rsp_rng)
mixtoxPlot(x_RfDmx, y_RfDmx, ca_RfDmx, ia_RfDmx, EffPoints, main = "RfD max", ylim = rsp_rng)




