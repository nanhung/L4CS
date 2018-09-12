library(mixtox)
data(cytotox)

compound <- c(rep("Ni", 12), rep("Zn", 12), rep("Cu", 12), rep("Mn", 12), 
              rep("Omim", 12), rep("Dmim", 12), rep("Emim", 12), rep("=Hmim", 12),
              rep("EE10 mixture", 12), rep("EE50 mixture", 12),
              rep("U1 mixture", 12), rep("U2 mixture", 12),rep("U3 mixture", 12),rep("U4 mixture", 12),rep("U5 mixture", 12),
              rep("U6 mixture", 12), rep("U7 mixture", 12),rep("U8 mixture", 12),rep("U9 mixture", 12),rep("U10 mixture", 12))


conc <- c(cytotox$Ni$x,
          cytotox$Zn$x,
          cytotox$Cu$x,
          cytotox$Mn$x,
          cytotox$Omim$x,
          cytotox$Dmim$x,
          cytotox$Emim$x,
          cytotox$Hmim$x,
          cytotox$ee10$x,
          cytotox$ee50$x,
          cytotox$u1$x,
          cytotox$u2$x,
          cytotox$u3$x,
          cytotox$u4$x,
          cytotox$u5$x,
          cytotox$u6$x,
          cytotox$u7$x,
          cytotox$u8$x,
          cytotox$u9$x,
          cytotox$u10$x)

expr <- c(rowMeans(cytotox$Ni$y),
          rowMeans(cytotox$Zn$y),
          rowMeans(cytotox$Cu$y),
          rowMeans(cytotox$Mn$y),
          rowMeans(cytotox$Omim$y),
          rowMeans(cytotox$Dmim$y),
          rowMeans(cytotox$Emim$y),
          rowMeans(cytotox$Hmim$y),
          rowMeans(cytotox$ee10$y),
          rowMeans(cytotox$ee50$y),
          rowMeans(cytotox$u1$y),
          rowMeans(cytotox$u2$y),
          rowMeans(cytotox$u3$y),
          rowMeans(cytotox$u4$y),
          rowMeans(cytotox$u5$y),
          rowMeans(cytotox$u6$y),
          rowMeans(cytotox$u7$y),
          rowMeans(cytotox$u8$y),
          rowMeans(cytotox$u9$y),
          rowMeans(cytotox$u10$y))

type <- c(rep("Single compound",96), rep("Equal effect design",24), rep("Uniform design", 120))
df<-data.frame(compound, conc, expr, type)
df$type <- factor(df$type, level=c("Uniform design", "Equal effect design", "Single compound"))
df$compound<-factor(df$compound, 
                    level=c("U1 mixture", "U2 mixture", "U3 mixture", "U4 mixture", "U5 mixture",
                            "U6 mixture", "U7 mixture", "U8 mixture", "U9 mixture", "U10 mixture",
                            "EE50 mixture", "EE10 mixture",
                            "Ni","Zn","Cu","Mn","Omim","Dmim","Emim","Hmim")) 

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
x1 <- cytotox$Ni$x
y1 <- cytotox$Ni$y
tuneFit1 <- tuneFit(x1, rowMeans(y1), eq = "Logit")

x2 <- cytotox$Zn$x
y2 <- cytotox$Zn$y
tuneFit2 <- tuneFit(x2, rowMeans(y2), eq = "Logit")

x3 <- cytotox$Cu$x
y3 <- cytotox$Cu$y
tuneFit3 <- tuneFit(x3, rowMeans(y3), eq = "Logit")

x4 <- cytotox$Mn$x
y4 <- cytotox$Mn$y
tuneFit4 <- tuneFit(x4, rowMeans(y4), eq = "Logit")

x5 <- cytotox$Omim$x
y5 <- cytotox$Omim$y
tuneFit5 <- tuneFit(x5, rowMeans(y5), eq = "Weibull")

x6 <- cytotox$Dmim$x
y6 <- cytotox$Dmim$y
tuneFit6 <- tuneFit(x6, rowMeans(y6), eq = "Weibull")

x7 <- cytotox$Emim$x
y7 <- cytotox$Emim$y
tuneFit7 <- tuneFit(x7, rowMeans(y7), eq = "Logit")

x8 <- cytotox$Hmim$x
y8 <- cytotox$Hmim$y
tuneFit8 <- tuneFit(x8, rowMeans(y8), eq = "Weibull")

fit1 <- curveFit(x1, y1, eq = "Logit", param = c(tuneFit1$sta[1], tuneFit1$sta[2]))
fit2 <- curveFit(x2, y2, eq = "Logit", param = c(tuneFit2$sta[1], tuneFit2$sta[2]))
fit3 <- curveFit(x3, y3, eq = "Logit", param = c(tuneFit3$sta[1], tuneFit3$sta[2]))
fit4 <- curveFit(x4, y4, eq = "Logit", param = c(tuneFit4$sta[1], tuneFit4$sta[2]))
fit5 <- curveFit(x5, y5, eq = "Weibull", param = c(tuneFit5$sta[1], tuneFit5$sta[2]))
fit6 <- curveFit(x6, y6, eq = "Weibull", param = c(tuneFit6$sta[1], tuneFit6$sta[2]))
fit7 <- curveFit(x7, y7, eq = "Logit", param = c(tuneFit7$sta[1], tuneFit7$sta[2]))
fit8 <- curveFit(x8, y8, eq = "Weibull", param = c(tuneFit8$sta[1], tuneFit8$sta[2]))

png(file="mixtox-2.png",width=3600,height=2000,res=300)
par(mfrow = c(2,4), oma = c(1,1,2,1))
figPlot(fit1, xlab="log(c) mol/L", ylab="Inhibition (%)"); 
mtext(paste0("Ni (Logit);"," r2 = ",round(tuneFit1$sta[3],3)))
figPlot(fit2, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("Zn (Logit);"," r2 = ",round(tuneFit2$sta[3],3)))
figPlot(fit3, xlab="log(c) mol/L", ylab="Inhibition (%)"); 
mtext(paste0("Cu (Logit);"," r2 = ",round(tuneFit3$sta[3],3)))
figPlot(fit4, xlab="log(c) mol/L", ylab="Inhibition (%)"); 
mtext(paste0("Mn (Logit);"," r2 = ",round(tuneFit4$sta[3],3)))
figPlot(fit5, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("Omim (Weibull);"," r2 = ",round(tuneFit5$sta[3],3)))
figPlot(fit6, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("Dmim (Weibull);"," r2 = ",round(tuneFit6$sta[3],3)))
figPlot(fit7, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("Emim (Logit);"," r2 = ",round(tuneFit7$sta[3],3)))
figPlot(fit8, xlab="log(c) mol/L", ylab="Inhibition (%)");
mtext(paste0("Hmim (Weibull);"," r2 = ",round(tuneFit8$sta[3],3)))
dev.off()

#####
model <- cytotox$sgl$model
param <- cytotox$sgl$param
eeca <- caPred(model, param, mixType = "eecr", effv = c(0.1, 0.5))
eeia <- iaPred(model, param, mixType = "eecr", effv = c(0.1, 0.5))
udca <- caPred(model, param, mixType = "udcr", effv = rep(c(0.05, 0.1, 0.2, 0.3, 0.5), 2))
udia <- iaPred(model, param, mixType = "udcr", effv = rep(c(0.05, 0.1, 0.2, 0.3, 0.5), 2))

df <- reshape2::melt(eeca$pct)
names(df) <- c("EC","chemical","percentage")
df <- ddply(df, .(EC), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

png(file="mixtox-3.png",width=3600,height=2800,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = EC, fill = chemical), data = df, stat="identity")+
  ggtitle("Percentage of individual chemicals in the EECR mixtures")+
  geom_text(data=df, aes(x = EC, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+
  ylab("Percentage (%)")
dev.off()

rownames(udca$unitab) <- c("u1","u2","u3","u4","u5","u6","u7","u8","u9","u10")
colnames(udca$unitab) <- c("l1","l2","l3","l4","l5","l6","l7","l8")
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

df <- reshape2::melt(cytotox$udcr.pct)
names(df) <- c("U","chemical","percentage")
df <- ddply(df, .(U), transform, pos = 1- (cumsum(percentage) - (0.5 * percentage)))

png(file="mixtox-5.png",width=3600,height=2800,res=300)
ggplot() + geom_bar(aes(y = percentage*100, x = U, fill = chemical), data = df, stat="identity")+
  ggtitle("Percentage of individual chemicals in the UDCR mixtures")+
  geom_text(data=df, aes(x = U, y = pos*100, label = paste0(round(percentage*100, 2),"%")), size=4) +
  xlab("Design")+
  ylab("Percentage (%)") 
dev.off()


########

png(file="mixtox-6.png",width=3600,height=2000,res=300)
par(mfrow = c(1, 2))
x <- cytotox$ee10$x
expr <- cytotox$ee10$y
#ee10fit <- curveFit(x, expr, eq = cytotox$eecr.mix$model[1], param = cytotox$eecr.mix$param[1, ])
plot(rep(log10(x), ncol(expr)), expr * 100, pch = 20, ylim = c(-10, 110),
     xlab = "log(c) mol/L", ylab = "Inhibition (%)", cex = 1.8, cex.lab = 1.8,
     cex.axis = 1.8)
#lines(log10(x), ee05fit$crcInfo[, 2] * 100, col = 1, lwd = 2)
#lines(log10(x), ee05fit$crcInfo[, 6] * 100, col = "blue", lwd = 1.5, lty = 3)
#lines(log10(x), ee05fit$crcInfo[, 7] * 100, col = "blue", lwd = 1.5,  lty = 3)
lines(log10(eeia$ia[1, ]), eeia$e * 100, col = "red", lwd = 2.5, lty = 2)
lines(log10(eeca$ca[1, ]), eeca$e * 100, col = "green", lwd = 2.5, lty = 2)
mtext("EE10", cex = 2)
legend('topleft', legend = c('Independent action (IA) prediction', 'Concentration action (CA) prediction'), 
       col = c('red','green'),
       lty = 2, lwd = 2.5, pch = NA, bty = 'n')

x <- cytotox$ee50$x
expr <- cytotox$ee50$y
#ee50fit <- curveFit(x, expr, eq = cytotox$eecr.mix$model[1], param = cytotox$eecr.mix$param[1, ])
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

Uplot <-function(i, X, title){
  x <- X$x
  expr <- X$y
  #fit <- curveFit(x, expr, eq = cytotox$udcr.mix$model[i], param = cytotox$udcr.mix$param[i, ])
  plot(rep(log10(x), ncol(expr)), expr * 100, pch = 20, ylim = c(-10, 110),
       xlab = "log(c) mol/L", ylab = "Inhibition (%)", cex = 1.8, cex.lab = 1.8,
       cex.axis = 1.8)
  #lines(log10(x), fit$crcInfo[, 2] * 100, col = 1, lwd = 2)
  #lines(log10(x), fit$crcInfo[, 6] * 100, col = "blue", lwd = 1.5, lty = 3)
  #lines(log10(x), fit$crcInfo[, 7] * 100, col = "blue", lwd = 1.5,  lty = 3)
  lines(log10(udia$ia[i, ]), eeia$e * 100, col = "red", lwd = 2.5, lty = 2)
  lines(log10(udca$ca[i, ]), eeca$e * 100, col = "green", lwd = 2.5, lty = 2)
  mtext(title, cex = 1.5)
}

png(file="mixtox-7.png",width=3600,height=2000,res=300)
par(mfrow = c(2, 5))
Uplot(1, cytotox$u1, "U1")
Uplot(2, cytotox$u2, "U2")
Uplot(3, cytotox$u3, "U3")
Uplot(4, cytotox$u4, "U4")
Uplot(5, cytotox$u5, "U5")
Uplot(6, cytotox$u6, "U6")
Uplot(7, cytotox$u7, "U7")
Uplot(8, cytotox$u8, "U8")
Uplot(9, cytotox$u9, "U9")
Uplot(10, cytotox$u10, "U10")
dev.off()
