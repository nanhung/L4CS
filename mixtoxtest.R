# devtools::install_github("lahothorn/SiTuR")
library(tukeytrend)
library(mixtox)
library(ggplot2)
library(dplyr)
library(plyr)
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

dfChemdose <- dfChemCoef %>% filter(term %in% "dose")

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
  geom_hline(yintercept=-log(0.05), col = "#800000", linetype = "dotdash")+
  annotate("text", x = -1.2, y = 3.5, label = "p = 0.05", col = "#800000")+
  geom_label(size = 2)
dev.off()

#tukeytrendfit
#df1 <- filter(DF1, chemical == "ALDRIN")
#fitw <- lm(response~dose, data = df1)
#ttw<-tukeytrendfit(fitw, dose="dose", scaling = c("ari", "ord", "log"))
#s<-summary(asglht(ttw))
#s$test$pvalues

df1 <- DF1 %>% filter(chemical == "ALDRIN") %>% filter(round == "r1")
df2 <- DF1 %>% filter(chemical == "ALDRIN") %>% filter(round == "r2")
conc <- df1$dose
response <- cbind(df1$response, df2$response)
response <- apply(response, 2, function(x) x/max(x))

tuneFit1 <- tuneFit(conc, rowMeans(response), eq = "Hill")
tuneFit1 <- tuneFit(conc, rowMeans(response), eq = "Logit")
tuneFit1 <- tuneFit(conc, rowMeans(response), eq = "Weibull")
tuneFit1

fit1 <- curveFit(conc, rowMeans(response), eq = "Logit", 
                 param = c(tuneFit1$sta[1], tuneFit1$sta[2]))

figPlot(fit1, xlab=expression(paste("Conc. (", mu,"M)")), ylab="Inhibition (%)")
mtext(paste0("ALDRIN (Logit);"," r2 = ",round(tuneFit1$sta[3],3)))


