# devtools::install_github("lahothorn/SiTuR")
#library(tukeytrend)
library(mixtox)
library(ggplot2)
library(dplyr)
#library(plyr)
library(tidyr) #seperate
library(broom)


df <- read.csv("mixtoxtest.csv")
colnames(df)<-c("chemical", "1_r1","1_r2", "2_r1","2_r2", 
            "3_r1","3_r2", "4_r1","4_r2", "5_r1","5_r2")
DF <- df %>% reshape::melt() %>% separate(variable, c("dose", "round")) 
DF$dosen <- as.numeric(DF$dose)
DF1 <- DF %>% mutate(dose = 10^(dosen-3))
colnames(DF1)[4]<-"response"

dfChem <- DF1 %>% group_by(chemical) %>%
  do(fitw = lm(response~dose, data = .))

dfChemCoef <- tidy(dfChem, fitw)

dfChemdose <- dfChemCoef %>% filter(term %in% "dose") %>% 
  mutate(screen = ifelse(p.value < 0.001 & estimate < 0, "Tox", "Non-Tox"))


png(file="dose_response.png",width=5200,height=2800,res=300)
ggplot(DF1, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 7) + 
  scale_x_log10() + theme_bw()
dev.off()

png(file="significant_test.png",width=3600,height=2800,res=300)
ggplot(dfChemdose, aes(x = estimate, y = -log(p.value), label = chemical)) +
  xlab(expression('slope'))+
  ylab("-log10 (p-value)") +
  geom_hline(yintercept=-log(0.001), col = "grey60", linetype = "dotted")+
  geom_hline(yintercept=-log(0.05), col = "grey60", linetype = "dashed")+
  annotate("text", x = -1.2, y = -log(0.001)+0.2, label = "p = 0.001", col = "grey60")+
  annotate("text", x = -1.2, y = -log(0.05)+0.2, label = "p = 0.05", col = "grey60")+
  geom_label(size = 2, aes(fill = screen), colour = "white", fontface = "bold") +
  theme(legend.position='none')
dev.off()

#####################

# tukeytrendfit

chem <- "MERCURIC CHLORIDE"
p.df<-matrix(NA, 42,3)
row.names(p.df) <- unique(DF1$chemical)
colnames(p.df) <- c("Non", "Rank", "Log")

for (chem in unique(DF1$chemical)){
  df1 <- filter(DF1, chemical %in% chem)
  fitw <- lm(response~dose, data = df1)
  ttw<-tukeytrendfit(fitw, dose="dose", scaling = c("ari", "ord", "log"))
  s<-summary(asglht(ttw))
  p.df[chem,] <- s$test$pvalues  
}

P.DF <- p.df %>% reshape::melt() %>%
  magrittr::set_colnames(c("chemical", "scaling", "p.value")) %>%
  mutate(logp = -log(p.value))

P.DF <- P.DF[!is.infinite(P.DF$logp),]

png(file="significant_test2.png",width=4800,height=2800,res=300)
ggplot(P.DF, aes(x = reorder(chemical, logp), y=logp, label = scaling)) +
  geom_label(size = 3, aes(fill = scaling), colour = "white", fontface = "bold") +
  geom_path(linetype = 2, color = "grey40") +
  geom_hline(yintercept=-log(0.05), col = "grey60", linetype = "dashed")+
  xlab(expression('chemical'))+
  ylab("-log10 (p-value)") +
  annotate("text", x = 1, y = -log(0.05)+0.2, label = "p = 0.05", col = "grey60")+
  coord_flip()
dev.off()

######################

chem <- c("DDT, O,P'-", "ALDRIN", "DIELDRIN", "CADMIUM(Chloride)", "HEPTACHLOR",
          "DDD, P,P'-", "MERCURIC CHLORIDE", "ENDOSULFAN", "DICOFOL", "DI(2-ETHYLHEXYL)PHTHALATE")
          #"HEPTACHLOR EPOXIDE", "CHLORPYRIFOS", "METHOXYCHLOR", 
          #"ENDRIN", "Potassium Chromate",
          #"DDT, P,P'-", "NICKEL",
          #"COBALT","BENZO(B)FLUORANTHEN")

par(mfrow = c(2,5))
param2fitplot("DDT, O,P'-", "Logit")
param2fitplot("ALDRIN", "Logit")
param2fitplot("DIELDRIN", "Logit")
param2fitplot("CADMIUM(Chloride)", "Logit")
param2fitplot("HEPTACHLOR", "Logit")
param2fitplot("DDD, P,P'-", "Logit")
#param2fitplot("MERCURIC CHLORIDE", "Weibull")
param2fitplot("ENDOSULFAN", "Logit")
param2fitplot("DICOFOL", "Logit")
#param2fitplot("DI(2-ETHYLHEXYL)PHTHALATE", "Logit")
param2fitplot("HEPTACHLOR EPOXIDE", "Logit")
param2fitplot("CHLORPYRIFOS", "Logit")


#detach(package:plyr)
library(scales)

Tox_chem <- dfChemdose %>% filter(screen == "Tox") 
Tox_chem$chemical

DF <- DF1 %>% filter(chemical == Tox_chem$chemical[i]) %>% 
  group_by(round) %>% mutate(normalize.response = response/max(response))
DR <- as.data.frame(DF %>% group_by(dose) %>% 
                      summarise(N.R = mean(normalize.response)))

ggplot(DF, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() + theme(legend.position='none')

ggplot(DF, aes(x = dose, y = normalize.response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_bw() + theme(legend.position='none')


FIT <- function(i){
  DF <- DF1 %>% filter(chemical == Tox_chem$chemical[i]) %>% 
    group_by(round) %>% mutate(normalize.response = response/max(response))
  DR <- as.data.frame(DF %>% group_by(dose) %>% 
                        summarise(N.R = mean(normalize.response)))
  
  print(Tox_chem$chemical[i])
  #print("Hill")
  #print(tuneFit(DR$dose, DR$N.R, eq = "Hill"))
  print("Logit")
  #print(tuneFit(DR$dose, DR$N.R, eq = "Logit"))
  #print("Weibull")
  #print(tuneFit(DR$dose, DR$N.R, eq = "Weibull"))
  #print("BCL")
  #print(tuneFit(DR$dose, DR$N.R, eq = "BCL"))
  #print("BCW")
  #print(tuneFit(DR$dose, DR$N.R, eq = "BCW"))
  #print("GL")
  #print(tuneFit(DR$dose, DR$N.R, eq = "GL"))
}

FIT(7)

#######################

Tox_chem <- dfChemdose %>% filter(screen == "Tox") 
param2fitplot <- function (chem, model){
  chem.no <- which(Tox_chem$chemical == chem)
  DF <- DF1 %>% filter(chemical == Tox_chem$chemical[chem.no]) %>% 
    group_by(round) %>% mutate(normalize.response = response/max(response))
  DR <- as.data.frame(DF %>% group_by(dose) %>% 
                        summarise(N.R = mean(normalize.response)))
  tuneFit1 <- tuneFit(DR$dose, DR$N.R, eq = model)
  fit1 <- curveFit(DR$dose, DR$N.R, eq = model, 
                   param = c(tuneFit1$sta[1], tuneFit1$sta[2]))
  figPlot(fit1, xlab=expression(paste("Log", "Conc. (", mu,"M)")), ylab="Inhibition (%)")
  mtext(paste(chem, ", R2 = ", round(tuneFit1$sta[4],3)))
}

par(mfrow = c(3,5))
param2fitplot("ALDRIN", "Logit")
param2fitplot("CADMIUM(Chloride)", "Logit")
param2fitplot("CHLORPYRIFOS", "Logit")
param2fitplot("DDD, P,P'-", "Logit")
param2fitplot("DDT, O,P'-", "Logit")
param2fitplot("DDT, P,P'-", "Logit")
param2fitplot("DICOFOL", "Logit")
param2fitplot("DIELDRIN", "Logit")
param2fitplot("ENDOSULFAN", "Logit")
param2fitplot("HEPTACHLOR", "Logit")
param2fitplot("HEPTACHLOR EPOXIDE", "Logit")
param2fitplot("METHOXYCHLOR", "Logit")
param2fitplot("NICKEL", "Logit")
param2fitplot("Potassium Chromate", "Logit")

########################

chem <- c("DDT, O,P'-", "ALDRIN", "DIELDRIN", "CADMIUM(Chloride)", "HEPTACHLOR",
          "DDD, P,P'-", "MERCURIC CHLORIDE", "ENDOSULFAN", "DICOFOL", "DI(2-ETHYLHEXYL)PHTHALATE")

Screen_df <- DF1 %>% filter(chemical %in% chem) %>% 
  group_by(round, chemical) %>% mutate(normalize.response = response/max(response))

ggplot(Screen_df, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 5) + 
  scale_x_log10() + theme_bw() + theme(legend.position = "top")

ggplot(Screen_df, aes(x = dose, y = normalize.response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 5) + 
  scale_x_log10() + theme_bw() + theme(legend.position = "top")



