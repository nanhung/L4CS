# devtools::install_github("lahothorn/SiTuR")
# devtools::install_github("nanhung/mixtox")

library(tukeytrend)

library(mixtox)
library(ggplot2)
library(dplyr)
#library(plyr)
library(tidyr) #seperate
library(broom)
library(scales)

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


png(file="dose_response.png",width=5200,height=3600,res=300)
ggplot(DF1, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 6) + theme_bw() +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))
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

###
p3<-ggplot(DF1, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 6) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme(legend.position = "none", rect = element_rect(fill = "transparent"))
g<-ggplotGrob(p3)

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

png(file="significant_test2.png",width=5600,height=4000,res=300)
ggplot(P.DF, aes(x = reorder(chemical, logp), y=logp, label = scaling)) +
  geom_label(size = 3, aes(fill = scaling), colour = "white", fontface = "bold") +
  geom_path(linetype = 2, color = "grey40") +
  geom_hline(yintercept=-log(0.05), col = "grey60", linetype = "dashed")+
  geom_hline(yintercept=-log(0.001), col = "grey60", linetype = "dotted")+
  xlab(expression('chemical'))+
  ylab("-log10 (p-value)") +
  annotate("text", x = 1, y = -log(0.05)+0.2, label = "p = 0.05", col = "grey60")+
  annotate("text", x = 1, y = -log(0.001)+0.2, label = "p = 0.001", col = "grey60")+
  annotation_custom(grob=g, xmin = 1, xmax = 26,
                    ymin=7, ymax=34)+
  coord_flip()
dev.off()

######################

chem <- c("DDT, O,P'-", "ALDRIN", "DIELDRIN", "CADMIUM(Chloride)", "HEPTACHLOR",
          "DDD, P,P'-", "MERCURIC CHLORIDE", "ENDOSULFAN", "DICOFOL", "CHLORPYRIFOS")
          #"HEPTACHLOR EPOXIDE", "DI(2-ETHYLHEXYL)PHTHALATE", "METHOXYCHLOR", 
          #"ENDRIN", "Potassium Chromate",
          #"DDT, P,P'-", "NICKEL",
          #"COBALT","BENZO(B)FLUORANTHEN")

Screen_df <- DF1 %>% filter(chemical %in% chem) %>% 
  group_by(round, chemical) %>% mutate(normalize.response = response/max(response))

png(file="chem10_non-norm.png",width=2400,height=3600,res=300)
ggplot(Screen_df, aes(x = dose, y = response)) +
  geom_point(aes(colour = round))+
  geom_path(aes(colour = round), size = 0.1) + 
  facet_wrap( ~ chemical, ncol = 2) + 
  theme_bw() + theme(legend.position = "top") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) 
dev.off()

##########################

Tox_chem <- dfChemdose %>% filter(screen == "Tox") 
param2fitplot <- function (chem, model){
  chem.no <- which(Tox_chem$chemical == chem)
  DF <- DF1 %>% filter(chemical == Tox_chem$chemical[chem.no]) %>% 
    #group_by(round) %>% mutate(normalize.response = response/100)
    group_by(round) %>% mutate(normalize.response = response/max(response))
  DR <- as.data.frame(DF %>% group_by(dose) %>% 
                        summarise(N.R = mean(normalize.response)))
  tuneFit1 <- tuneFit(DR$dose, DR$N.R, eq = model)
  fit1 <- curveFit(DR$dose, DR$N.R, eq = model, 
                   param = c(tuneFit1$sta[1], tuneFit1$sta[2]))
  figPlot(fit1, xlab=expression(paste("Log", "Conc. (", mu,"M)")), ylab="Response (%)")
  mtext(paste(chem, ", R2 = ", round(tuneFit1$sta[4],3)))
}

png(file="Fit.png",width=4800,height=2200,res=300)
#par(mfrow = c(2,3), oma=c(0,0,2,0))
par(mfrow = c(2,5), oma=c(0,0,2,0))
param2fitplot("DDT, O,P'-", "Hill") # H
param2fitplot("ALDRIN", "Hill") # L
param2fitplot("DIELDRIN", "Hill")
param2fitplot("CADMIUM(Chloride)", "Hill") # H
param2fitplot("HEPTACHLOR", "Hill")
param2fitplot("DDD, P,P'-", "Hill")
#param2fitplot("MERCURIC CHLORIDE", "Weibull")
param2fitplot("ENDOSULFAN", "Hill") # H
param2fitplot("DICOFOL", "Hill") #L
param2fitplot("HEPTACHLOR EPOXIDE", "Hill")
param2fitplot("CHLORPYRIFOS", "Hill") # H
dev.off()

######################
#detach(package:plyr)

chem <- c("DDT, O,P'-", "ALDRIN", "DIELDRIN", "CADMIUM(Chloride)", "HEPTACHLOR",
          "DDD, P,P'-", "MERCURIC CHLORIDE", "ENDOSULFAN", "DICOFOL", "HEPTACHLOR EPOXIDE",
          "CHLORPYRIFOS")

FIT <- function(i){
  DF <- DF1 %>% filter(chemical == chem[i]) %>% 
    group_by(round) %>% mutate(normalize.response = response/max(response))
  DR <- as.data.frame(DF %>% group_by(dose) %>% 
                        summarise(N.R = mean(normalize.response)))
  
  print(chem[i])
  print("Hill")
  print(tuneFit(DR$dose, DR$N.R, eq = "Hill"))
  print("Logit")
  print(tuneFit(DR$dose, DR$N.R, eq = "Logit"))
  print("Weibull")
  print(tuneFit(DR$dose, DR$N.R, eq = "Weibull"))
  x<-tuneFit(DR$dose, DR$N.R, eq = "Hill", effv = 0.5)
  return(x)
}

x01<-FIT(1)
x02<-FIT(2)
x03<-FIT(3)
x04<-FIT(4)
x05<-FIT(5)
x06<-FIT(6)
x07<-FIT(8)
x08<-FIT(9)
x09<-FIT(10)
x10<-FIT(11)

#######################
library(plyr)

model <- rep("Hill",10)
a<-c(x01$sta[1],x01$sta[2],0,
     x02$sta[1],x02$sta[2],0,
     x03$sta[1],x03$sta[2],0,
     x04$sta[1],x04$sta[2],0,
     x05$sta[1],x05$sta[2],0,
     x06$sta[1],x06$sta[2],0,
     x07$sta[1],x07$sta[2],0,
     x08$sta[1],x08$sta[2],0,
     x09$sta[1],x09$sta[2],0,
     x10$sta[1],x10$sta[2],0)
param <- matrix(a, nrow = 10, ncol = 3, byrow = T)
colnames(param) <- c("alpha", "beta", "gamma")
row.names(param) <- chem[c(1:6,8:11)]

ECX <- ECx(model, param, effv = c(0.95, 0.90, 0.80, 0.70, 0.50))

ECDF <- reshape2::melt(ECX)
names(ECDF) <- c("chemical","ECx","conc.")

png(file="EC.png",width=4800,height=2200,res=300)
ggplot(data=ECDF) + 
  geom_text(aes(x = ECx, y = conc.+3, label = paste0(round(conc., 2)))) +
  geom_point(aes(x=ECx,y=conc.)) +
  geom_line(aes(x=as.numeric(ECx),y=conc.)) +
  facet_wrap(~reorder(chemical, conc.), nrow = 2) +
  ylab(expression(paste("Concentration (", mu,"M)")))
dev.off()


aca <- caPred(model, param, mixType = "acr", effv = c(rep(0.5, 10)))
aia <- iaPred(model, param, mixType = "acr", effv = c(rep(0.5, 10)))
eeca <- caPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))
eeia <- iaPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))
udca <- caPred(model, param, mixType = "udcr", effv = rep(c(0.95, 0.90, 0.80, 0.70, 0.50), 2))
udia <- iaPred(model, param, mixType = "udcr", effv = rep(c(0.95, 0.90, 0.80, 0.70, 0.50), 2))


df0 <- reshape2::melt(aca$pct)
names(df0) <- c("Exp","chemical","percentage")
df0[,2] <- chem[c(1:6,8:11)]
df0[,1] <- rep("ECx", 10)

df1 <- reshape2::melt(eeca$pct)
names(df1) <- c("EC","chemical","percentage")
df1 <- ddply(df1, .(EC), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

png(file="ACR.png",width=3600,height=2400,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = Exp, fill = chemical), data = df0, stat="identity")+
  ggtitle("Percentage of individual chemicals in the ACR mixtures")+
  xlab("Design")+
  ylab("Percentage (%)")
dev.off()

png(file="EECR.png",width=3600,height=2400,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = EC, fill = chemical), data = df1, stat="identity")+
  ggtitle("Percentage of individual chemicals in the EECR mixtures")+
  geom_text(data=df1, aes(x = EC, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+
  ylab("Percentage (%)")
dev.off()

rownames(udca$unitab) <- c("u1","u2","u3","u4","u5","u6","u7","u8","u9","u10")
#colnames(udca$unitab) <- c("l1","l2","l3","l4","l5","l6","l7","l8","l9","l10")
colnames(udca$unitab) <- chem[c(1:6,8:11)]



dat <- reshape2::melt(udca$unitab) %>%
  mutate(EC = ifelse(value %in% 1:2, "EC95",
                     ifelse(value %in% 3:4, "EC90",
                            ifelse(value %in% 5:6, "EC80",
                                   ifelse(value %in% 7:8, "EC70",
                                          ifelse(value %in% 9:10, "EC50", NA))))))

for(i in 1:nrow(dat)){
  for(j in colnames(udca$unitab)){
    for(k in unique(dat$EC)){
      if(dat$Var2[i] == j & dat$EC[i] == k){
        dat$ECC[i] <- ECX[j,k]
      }
    }
  }
}


png(file="UDCR.png",width=3600,height=2400,res=300)
ggplot(data =  dat, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = EC), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Uniform Design Table (ECx)") +
  xlab("number of runs (levels or pseudo-levels)") +
  ylab("chemical") +
  guides(fill=FALSE)
dev.off()

png(file="UDCR-C.png",width=3600,height=2400,res=300)
ggplot(data =  dat, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = round(ECC,2)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Uniform Design Table (concentration)") +
  xlab("number of runs (levels or pseudo-levels)") +
  ylab("chemical") +
  guides(fill=FALSE)
dev.off()

row.names(udca$pct) <- paste("u", 1:10, sep = "")

df2 <- reshape2::melt(udca$pct)
names(df2) <- c("U","chemical","percentage")
df2 <- ddply(df2, .(U), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

png(file="UDCR_pct.png",width=3600,height=2400,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = U, fill = chemical), data = df2, stat="identity")+
  ggtitle("Percentage of individual chemicals in the UDCR mixtures")+
  geom_text(data=df2, aes(x = U, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+
  ylab("Percentage (%)") 
dev.off()

########################

png(file="MixTox.png",width=3600,height=1800,res=300)
par(mfrow = c(1,2))
plot(eeca$ca[1, ], eeca$e * 100, type = "l", xlab = expression(paste("Concentration (", mu,"M)")), 
     ylab = "Percentage (%)", main = "Equal-effect design")
lines(eeca$ca[2, ], eeca$e * 100, col = i)
legend("topright", legend=c("EC05", "EC50"), col=c(1, 2), lty = 1, cex=0.8)

plot(udca$ca[1, ], eeia$e * 100, type = "l", xlab = expression(paste("Concentration (", mu,"M)")), 
     ylab = "Percentage (%)", main = "Uniform design")
for(i in 2:10){
  lines(udca$ca[i, ], eeia$e * 100, col = i)
}
legend("topright", legend=paste("U", 1 : 10, sep = ""), col=c(1:10), lty = 1, cex=0.8)
dev.off()

