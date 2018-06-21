library(mixtox)
library(ggplot2)
library(dplyr)
library(plyr)

compound <- c(rep("PAR", 12), rep("SPE", 12), rep("KAN", 12), rep("STR", 12), rep("DIH", 12), rep("GEN", 12), rep("NEO", 12),
              rep("EE05 mixture", 12), rep("EE50 mixture", 12),
              rep("U1 mixture", 12), rep("U2 mixture", 12),rep("U3 mixture", 12),rep("U4 mixture", 12),rep("U5 mixture", 12),
              rep("U6 mixture", 12), rep("U7 mixture", 12),rep("U8 mixture", 12),rep("U9 mixture", 12),rep("U10 mixture", 12))

conc <- c(antibiotox$PAR$x,
          antibiotox$SPE$x,
          antibiotox$KAN$x,
          antibiotox$STR$x,
          antibiotox$DIH$x,
          antibiotox$GEN$x,
          antibiotox$NEO$x,
          antibiotox$ee05$x,
          antibiotox$ee50$x,
          antibiotox$u1$x,
          antibiotox$u2$x,
          antibiotox$u3$x,
          antibiotox$u4$x,
          antibiotox$u5$x,
          antibiotox$u6$x,
          antibiotox$u7$x,
          antibiotox$u8$x,
          antibiotox$u9$x,
          antibiotox$u10$x)

expr <- c(rowMeans(antibiotox$PAR$y),
          rowMeans(antibiotox$SPE$y),
          rowMeans(antibiotox$KAN$y),
          rowMeans(antibiotox$STR$y),
          rowMeans(antibiotox$DIH$y),
          rowMeans(antibiotox$GEN$y),
          rowMeans(antibiotox$NEO$y),
          rowMeans(antibiotox$ee05$y),
          rowMeans(antibiotox$ee50$y),
          rowMeans(antibiotox$u1$y),
          rowMeans(antibiotox$u2$y),
          rowMeans(antibiotox$u3$y),
          rowMeans(antibiotox$u4$y),
          rowMeans(antibiotox$u5$y),
          rowMeans(antibiotox$u6$y),
          rowMeans(antibiotox$u7$y),
          rowMeans(antibiotox$u8$y),
          rowMeans(antibiotox$u9$y),
          rowMeans(antibiotox$u10$y))

type <- c(rep("Single compound",84), rep("Equal effect design",24), rep("Uniform design", 120))

df<-data.frame(compound, conc, expr, type)
df$type <- factor(df$type, level=c("Uniform design", "Equal effect design", "Single compound"))
df$compound<-factor(df$compound, 
                    level=c("U1 mixture", "U2 mixture", "U3 mixture", "U4 mixture", "U5 mixture",
                            "U6 mixture", "U7 mixture", "U8 mixture", "U9 mixture", "U10 mixture",
                            "EE50 mixture", "EE05 mixture",
                            "PAR","SPE","KAN","STR","DIH","GEN","NEO")) 

png(file="mixtox-1.png",width=3600,height=2800,res=300)
ggplot(df, aes(x=conc, y = expr)) +
  geom_line(aes(color=compound)) + 
  geom_point(aes(color=compound)) +
  facet_grid(type~.) +
  ggtitle("Experiment Data")+
  xlab("Concentration, mol/L")+
  ylab("Inhibition (%)")+
  ylim(min(df$expr),1)+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x, n= 4),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) + guides(size=FALSE)
dev.off()

#
x1 <- antibiotox$PAR$x
y1 <- antibiotox$PAR$y
tuneFit1 <- tuneFit(x1, rowMeans(y1), eq = "Logit")

x2 <- antibiotox$SPE$x
y2 <- antibiotox$SPE$y
tuneFit2 <- tuneFit(x1, rowMeans(y2), eq = "Weibull")

x3 <- antibiotox$KAN$x
y3 <- antibiotox$KAN$y
tuneFit3 <- tuneFit(x3, rowMeans(y3), eq = "Weibull")

x4 <- antibiotox$STR$x
y4 <- antibiotox$STR$y
tuneFit4 <- tuneFit(x4, rowMeans(y4), eq = "Weibull")

x5 <- antibiotox$DIH$x
y5 <- antibiotox$DIH$y
tuneFit5 <- tuneFit(x5, rowMeans(y5), eq = "Logit")

x6 <- antibiotox$GEN$x
y6 <- antibiotox$GEN$y
tuneFit6 <- tuneFit(x6, rowMeans(y6), eq = "Weibull")

x7 <- antibiotox$NEO$x
y7 <- antibiotox$NEO$y
tuneFit7 <- tuneFit(x7, rowMeans(y7), eq = "Weibull")

#

fit1 <- curveFit(x1, y1, eq = "Logit", param = c(tuneFit1$sta[1], tuneFit1$sta[2]), effv = c(0.05, 0.5))
fit2 <- curveFit(x2, y2, eq = "Weibull", param = c(tuneFit2$sta[1], tuneFit2$sta[2]), effv = c(0.05, 0.5))
fit3 <- curveFit(x3, y3, eq = "Weibull", param = c(tuneFit3$sta[1], tuneFit3$sta[2]), effv = c(0.05, 0.5))
fit4 <- curveFit(x4, y4, eq = "Weibull", param = c(tuneFit4$sta[1], tuneFit4$sta[2]), effv = c(0.05, 0.5))
fit5 <- curveFit(x5, y5, eq = "Logit", param = c(tuneFit5$sta[1], tuneFit5$sta[2]), effv = c(0.05, 0.5))
fit6 <- curveFit(x6, y6, eq = "Weibull", param = c(tuneFit6$sta[1], tuneFit6$sta[2]), effv = c(0.05, 0.5))
fit7 <- curveFit(x7, y7, eq = "Weibull", param = c(tuneFit7$sta[1], tuneFit7$sta[2]), effv = c(0.05, 0.5))

png(file="mixtox-2.png",width=3600,height=2800,res=300)
par(mfrow = c(2,4), oma = c(1,1,2,1))
figPlot(fit1, xlab="log(c) mol/L", ylab="Inhibition (%)"); 
mtext(paste0("PAR (Logit);"," r2 = ",round(tuneFit1$sta[3],3)))
figPlot(fit2, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("SPE (Weibull);"," r2 = ",round(tuneFit2$sta[3],3)))
figPlot(fit3, xlab="log(c) mol/L", ylab="Inhibition (%)"); 
mtext(paste0("KAN (Weibull);"," r2 = ",round(tuneFit3$sta[3],3)))
figPlot(fit4, xlab="log(c) mol/L", ylab="Inhibition (%)"); 
mtext(paste0("STR (Weibull);"," r2 = ",round(tuneFit4$sta[3],3)))
figPlot(fit5, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("DIH (Logit);"," r2 = ",round(tuneFit5$sta[3],3)))
figPlot(fit6, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("Gen (Weibull);"," r2 = ",round(tuneFit6$sta[3],3)))
figPlot(fit7, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("NEO (Weibull);"," r2 = ",round(tuneFit7$sta[3],3)))
plot.new()
legend('top', legend = c('95% CI', '95% PI'), col = c('red','blue'),
       lty = 1, lwd = 1, pch = NA, bty = 'n', cex = 2,
       text.col = 'black')
dev.off()

# acr: arbitrary concentration ratio; 
# eecr: equal effect concentration ratio; 
# udcr: uniform design concentration ratio.

model <- antibiotox$sgl$model
param <- antibiotox$sgl$param
eeca <- caPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))
eeia <- iaPred(model, param, mixType = "eecr", effv = c(0.05, 0.5))
udca <- caPred(model, param, mixType = "udcr", effv = rep(c(0.05, 0.1, 0.2, 0.3, 0.5), 2))
udia <- iaPred(model, param, mixType = "udcr", effv = rep(c(0.05, 0.1, 0.2, 0.3, 0.5), 2))

# http://t-redactyl.io/blog/2016/01/creating-plots-in-r-using-ggplot2-part-4-stacked-bar-plots.html
df <- reshape2::melt(eeca$pct)
names(df) <- c("EC","aminoglycoside","percentage")
df <- ddply(df, .(EC), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

png(file="mixtox-3.png",width=3600,height=2800,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = EC, fill = aminoglycoside), data = df, stat="identity")+
  ggtitle("Percentage of individual chemicals in the EECR mixtures")+
  geom_text(data=df, aes(x = EC, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+
  ylab("Percentage (%)")
dev.off()

#
rownames(udca$unitab) <- c("u1","u2","u3","u4","u5","u6","u7","u8","u9","u10")
colnames(udca$unitab) <- c("l1","l2","l3","l4","l5","l6","l7")
dat <- reshape2::melt(udca$unitab) %>%
  mutate(EC = ifelse(value %in% 1:2, "EC05",
                     ifelse(value %in% 3:4, "EC10",
                            ifelse(value %in% 5:6, "EC20",
                                   ifelse(value %in% 7:8, "EC30",
                                          ifelse(value %in% 9:10, "EC50", NA))))))


png(file="mixtox-4.png",width=3600,height=2800,res=300)
ggplot(data =  dat, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = value), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Uniform Design Table") +
  xlab("number of runs (levels or pseudo-levels)") +
  ylab("number of factors (compounds)") +
  guides(fill=FALSE)
dev.off()

png(file="mixtox-4-1.png",width=3600,height=2800,res=300)
ggplot(data =  dat, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), colour = "white") +
  geom_text(aes(label = EC), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ggtitle("Uniform Design Table") +
  xlab("number of runs (levels or pseudo-levels)") +
  ylab("number of factors (compounds)") +
  guides(fill=FALSE)
dev.off()



df <- reshape2::melt(antibiotox$udcr.pct)
names(df) <- c("U","aminoglycoside","percentage")
df <- ddply(df, .(U), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

png(file="mixtox-5.png",width=3600,height=2800,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = U, fill = aminoglycoside), data = df, stat="identity")+
  ggtitle("Percentage of individual chemicals in the UDCR mixtures")+
  geom_text(data=df, aes(x = U, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+
  ylab("Percentage (%)") 
dev.off()

#

png(file="mixtox-6.png",width=3600,height=2800,res=300)
par(mfrow = c(1, 2))
x <- antibiotox$ee05$x
expr <- antibiotox$ee05$y
ee05fit <- curveFit(x, expr, eq = antibiotox$eecr.mix$model[1],
                    param = antibiotox$eecr.mix$param[1, ])
plot(rep(log10(x), ncol(expr)), expr * 100, pch = 20, ylim = c(-10, 110),
     xlab = "log(c) mol/L", ylab = "Inhibition (%)", cex = 1.8, cex.lab = 1.8,
     cex.axis = 1.8)
#lines(log10(x), ee05fit$crcInfo[, 2] * 100, col = 1, lwd = 2)
#lines(log10(x), ee05fit$crcInfo[, 6] * 100, col = "blue", lwd = 1.5, lty = 3)
#lines(log10(x), ee05fit$crcInfo[, 7] * 100, col = "blue", lwd = 1.5,  lty = 3)
lines(log10(eeia$ia[1, ]), eeia$e * 100, col = "red", lwd = 2.5, lty = 2)
lines(log10(eeca$ca[1, ]), eeca$e * 100, col = "green", lwd = 2.5, lty = 2)
mtext("EE05", cex = 2)
legend('topleft', legend = c('Independent action (IA) prediction', 'Concentration action (CA) prediction'), 
       col = c('red','green'),
       lty = 2, lwd = 2.5, pch = NA, bty = 'n')

x <- antibiotox$ee50$x
expr <- antibiotox$ee50$y
ee50fit <- curveFit(x, expr, eq = antibiotox$eecr.mix$model[1],
                    param = antibiotox$eecr.mix$param[1, ])
plot(rep(log10(x), ncol(expr)), expr * 100, pch = 20, ylim = c(-10, 110),
     xlab = "log(c) mol/L", ylab = "Inhibition (%)", cex = 1.8, cex.lab = 1.8,
     cex.axis = 1.8)
#lines(log10(x), ee50fit$crcInfo[, 2] * 100, col = 1, lwd = 2)
#lines(log10(x), ee50fit$crcInfo[, 6] * 100, col = "blue", lwd = 1.5, lty = 3)
#lines(log10(x), ee50fit$crcInfo[, 7] * 100, col = "blue", lwd = 1.5,  lty = 3)
lines(log10(eeia$ia[2, ]), eeia$e * 100, col = "red", lwd = 2.5, lty = 2)
lines(log10(eeca$ca[2, ]), eeca$e * 100, col = "green", lwd = 2.5, lty = 2)
mtext("EE50", cex = 2)
dev.off()

#

Uplot <-function(i, X, title){
  x <- X$x
  expr <- X$y
  fit <- curveFit(x, expr, eq = antibiotox$udcr.mix$model[i], param = antibiotox$udcr.mix$param[i, ])
  plot(rep(log10(x), ncol(expr)), expr * 100, pch = 20, ylim = c(-10, 110),
       xlab = "log(c) mol/L", ylab = "Inhibition (%)", cex = 1.8, cex.lab = 1.8,
       cex.axis = 1.8)
  #lines(log10(x), fit$crcInfo[, 2] * 100, col = 1, lwd = 2)
  #lines(log10(x), fit$crcInfo[, 6] * 100, col = "blue", lwd = 1.5, lty = 3)
  #lines(log10(x), fit$crcInfo[, 7] * 100, col = "blue", lwd = 1.5,  lty = 3)
  lines(log10(udia$ia[1, ]), eeia$e * 100, col = "red", lwd = 2.5, lty = 2)
  lines(log10(udca$ca[1, ]), eeca$e * 100, col = "green", lwd = 2.5, lty = 2)
  mtext(title, cex = 1.5)
}

png(file="mixtox-7.png",width=3600,height=2800,res=300)
par(mfrow = c(2, 5))
Uplot(1, antibiotox$u1, "U1")
legend('topleft', legend = c('IA prediction', 'CA prediction'), col = c('red','green'),
       lty = c(3,2,2), lwd = c(1.5,2.5,2.5), 
       pch = NA, bty = 'n', cex = 1.2,
       text.col = 'black')
Uplot(2, antibiotox$u2, "U2")
Uplot(3, antibiotox$u3, "U3")
Uplot(4, antibiotox$u4, "U4")
Uplot(5, antibiotox$u5, "U5")
Uplot(6, antibiotox$u6, "U6")
Uplot(7, antibiotox$u7, "U7")
Uplot(8, antibiotox$u8, "U8")
Uplot(9, antibiotox$u9, "U9")
Uplot(10, antibiotox$u10, "U10")
dev.off()
