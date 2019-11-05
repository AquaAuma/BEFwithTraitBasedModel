rm(list=ls())
library(vegan)
library(stringr)
library(readr)
library(viridis)
library(ggplot2)
library(egg)
library(FactoMineR)


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S1: the Trait space
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

data <- read.csv('SimulationsSmallAlpha.csv')
datA <- data.frame()
for(i in 1:40){
  sub.dat <- subset(data, data$sr==i)
  s <- sample(c(1:nrow(sub.dat)), size=50)
  sub.sub <- sub.dat[s,]
  datA <- rbind(datA, sub.sub)
}
max.resource <- 510*101
datA$biomass <- datA$biomass/max.resource
datA$gainr <- datA$gainr/max.resource
datA$metab <- datA$metab/max.resource

library(ggplot2)
ggNS <- ggplot(datA, aes(x=sr, y=NS01)) + geom_point(color=adjustcolor('grey50', alpha=0.1)) +
  theme(panel.grid.minor = element_blank())+ theme_bw() + 
  xlab('Species Richness') + ylab('Niche Space') + 
  theme(axis.text=element_text(size=18),axis.title=element_text(size=18))+
  scale_x_continuous(limits = c(0, 40))+scale_y_continuous(limits = c(0, 1)) +
  geom_smooth(aes(x=sr, y=NS01), method='loess', span=1, lwd=1, color='black', fill=adjustcolor('grey', alpha=0.3))+
  
  geom_point(aes(x=sr, y=NS001), color=adjustcolor('blue', alpha=0.1)) + 
  geom_smooth(aes(x=sr, y=NS001), method='loess', span=1, lwd=1, color='blue', fill=adjustcolor('blue', alpha=0.3))+
  
  geom_point(aes(x=sr, y=NS0001), color=adjustcolor('indianred2', alpha=0.1)) + 
  geom_smooth(aes(x=sr, y=NS0001), method='loess', span=1, lwd=1, color='indianred', fill=adjustcolor('indianred2', alpha=0.3))+
  
  geom_point(aes(x=sr, y=NS00001), color=adjustcolor('mediumseagreen', alpha=0.1)) + 
  geom_smooth(aes(x=sr, y=NS00001), method='loess', span=1, lwd=1, color='mediumseagreen', fill=adjustcolor('mediumseagreen', alpha=0.3))

windows(80,60)
ggNS


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S2: Standardization on food web simulations
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

data <- read.csv('SimulationsSmallAlpha.csv')
datc <- read.csv('SimulationsHighAlpha.csv')

histA <- ggplot(data, aes(x=sr)) + geom_histogram(color='black', fill='grey', binwidth=1)+
  theme(panel.grid.minor = element_blank())+ theme_bw() + 
  xlab('Species Richness') + ylab('Simulations') + theme(axis.text=element_text(size=18),axis.title=element_text(size=18))+
  scale_x_continuous(limits = c(0, 50))+scale_y_continuous(limits = c(0, 7500))

histC <- ggplot(datc, aes(x=sr)) + geom_histogram(color='black', fill='grey', binwidth=1)+
  theme(panel.grid.minor = element_blank())+ theme_bw() + 
  xlab('Species Richness') + ylab('Simulations') + theme(axis.text=element_text(size=18),axis.title=element_text(size=18))+
  scale_x_continuous(limits = c(0, 50))+scale_y_continuous(limits = c(0, 7500))

windows(100,40)
egg::ggarrange(histA, histC, ncol = 2, nrow = 1)



###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S3: Shapes of richness-ecosystme functioning relationships
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

dat2 <- read.csv('SimulationsSmallAlpha.csv')

# With 40 species max, we can randomly select at least 100 species per richness level to build the boxplots
dat3 <- data.frame()
for(i in 1:40){
  sub.dat <- subset(dat2, dat2$sr==i)
  s <- sample(c(1:nrow(sub.dat)), size=50)
  sub.sub <- sub.dat[s,]
  dat3 <- rbind(dat3, sub.sub)
}

max.resource <- 101*510
dat3$biomass <- dat3$biomass/max.resource
dat3$gainr <- dat3$gainr/max.resource
dat3$metab <- dat3$metab/max.resource


### Biomass ~ Richness
nlsB.mm <- nls(biomass ~ a*sr/(b+sr), data=dat3, start=list(a=1.1, b=2), algorithm = "default", trace = T)
nlsB.p <- nls(biomass ~ a*sr^b, data=dat3, start=list(a=0.3, b=0.35), algorithm = "default", trace = T)
nlsB.linear <- nls(biomass ~ a*sr+b, data=dat3, start=list(a=0.5, b=0.15), algorithm = "default", trace = T)
nlsB.sigm <- nls(biomass ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=0.2, b=2, c=0.2), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsB.mm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$biomass-mean(dat3$biomass))^2)
var.resid.M0 <- (1/n)*sum((dat3$biomass - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:40))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsB.mm)[1]*pred.models$SR/(coef(nlsB.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsB.p)[1]*(pred.models$SR^(coef(nlsB.p)[2]))
pred.models$L <- coef(nlsB.linear)[1]*pred.models$SR + coef(nlsB.linear)[2]
pred.models$S <- coef(nlsB.sigm)[1]/(1+coef(nlsB.sigm)[2]*exp(-coef(nlsB.sigm)[3]*pred.models$SR))

ggfit.B <- ggplot(dat3, aes(x=as.factor(sr), y=biomass)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + xlab('Species Richness')+ylab('Biomass')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25)

### Metabolism ~ Richness
nlsM.mm <- nls(metab ~ a*sr/(b+sr), data=dat3, start=list(a=0.7, b=0.5), algorithm = "default", trace = T)
nlsM.p <- nls(metab ~ a*sr^b, data=dat3, start=list(a=0.5, b=0.08), algorithm = "default", trace = T)
nlsM.linear <- nls(metab ~ a*sr+b, data=dat3, start=list(a=0.1, b=0.1), algorithm = "default", trace = T)
nlsM.sigm <- nls(metab ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=0.15, b=0.5, c=0.5), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsM.sigm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$metab-mean(dat3$metab))^2)
var.resid.M0 <- (1/n)*sum((dat3$metab - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:40))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsM.mm)[1]*pred.models$SR/(coef(nlsM.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsM.p)[1]*(pred.models$SR^(coef(nlsM.p)[2]))
pred.models$L <- coef(nlsM.linear)[1]*pred.models$SR + coef(nlsM.linear)[2]
pred.models$S <- coef(nlsM.sigm)[1]/(1+coef(nlsM.sigm)[2]*exp(-coef(nlsM.sigm)[3]*pred.models$SR))

ggfit.M <- ggplot(dat3, aes(x=as.factor(sr), y=metab)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + xlab('Species Richness')+ylab('Metabolism')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25)

### Production ~ Richness
nlsP.mm <- nls(gainr ~ a*sr/(b+sr), data=dat3, start=list(a=2.2, b=5), algorithm = "default", trace = T)
nlsP.p <- nls(gainr ~ a*sr^b, data=dat3, start=list(a=0.2, b=0.6), algorithm = "default", trace = T)
nlsP.linear <- nls(gainr ~ a*sr+b, data=dat3, start=list(a=0.1, b=0.1), algorithm = "default", trace = T)
nlsP.sigm <- nls(gainr ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=0.5, b=4,c=0.25), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsP.sigm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$gainr-mean(dat3$gainr))^2)
var.resid.M0 <- (1/n)*sum((dat3$gainr - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:40))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsP.mm)[1]*pred.models$SR/(coef(nlsP.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsP.p)[1]*(pred.models$SR^(coef(nlsP.p)[2]))
pred.models$L <- coef(nlsP.linear)[1]*pred.models$SR + coef(nlsP.linear)[2]
pred.models$S <- coef(nlsP.sigm)[1]/(1+coef(nlsP.sigm)[2]*exp(-coef(nlsP.sigm)[3]*pred.models$SR))

ggfit.P <- ggplot(dat3, aes(x=as.factor(sr), y=gainr)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + xlab('Species Richness')+ylab('Production')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25)


### Productivity ~ Richness
nlsPB.mm <- nls(PB1 ~ a*sr/(b+sr), data=dat3, start=list(a=2, b=1), algorithm = "default", trace = T)
nlsPB.p <- nls(PB1 ~ a*sr^b, data=dat3, start=list(a=1, b=0.2), algorithm = "default", trace = T)
nlsPB.linear <- nls(PB1 ~ a*sr+b, data=dat3, start=list(a=0.5, b=1), algorithm = "default", trace = T)
nlsPB.sigm <- nls(PB1 ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=2, b=2, c=0.25), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsPB.sigm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$PB1-mean(dat3$PB1))^2)
var.resid.M0 <- (1/n)*sum((dat3$PB1 - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:40))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsPB.mm)[1]*pred.models$SR/(coef(nlsPB.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsPB.p)[1]*(pred.models$SR^(coef(nlsPB.p)[2]))
pred.models$L <- coef(nlsPB.linear)[1]*pred.models$SR + coef(nlsPB.linear)[2]
pred.models$S <- coef(nlsPB.sigm)[1]/(1+coef(nlsPB.sigm)[2]*exp(-coef(nlsPB.sigm)[3]*pred.models$SR))

ggfit.PB <- ggplot(dat3, aes(x=as.factor(sr), y=PB1)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  xlab('Species Richness')+ylab('Productivity')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25) 

library(egg)
windows()
ggarrange(ggfit.B, ggfit.M, ggfit.P, ggfit.PB, labels=c('A','B','C','D'), nrow=2, ncol=2)


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S3b: Shapes of richness-ecosystme functioning relationships
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

datc <- read.csv('SimulationsHighAlpha.csv')
datC <- data.frame()
for(i in 1:40){
  sub.dat <- subset(datc, datc$sr==i)
  s <- sample(c(1:nrow(sub.dat)), size=50)
  sub.sub <- sub.dat[s,]
  datC <- rbind(datC, sub.sub)
}
datC$biomass <- datC$biomass/10000
datC$gainr <- datC$gainr/10000
datC$metab <- datC$metab/10000
dat3 <- datC

### Biomass ~ Richness
nlsB.mm <- nls(biomass ~ a*sr/(b+sr), data=dat3, start=list(a=1.1, b=2), algorithm = "default", trace = T)
nlsB.p <- nls(biomass ~ a*sr^b, data=dat3, start=list(a=3, b=0.35), algorithm = "default", trace = T)
nlsB.linear <- nls(biomass ~ a*sr+b, data=dat3, start=list(a=0.5, b=0.15), algorithm = "default", trace = T)
nlsB.sigm <- nls(biomass ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=0.2, b=2, c=0.2), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsB.sigm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$biomass-mean(dat3$biomass))^2)
var.resid.M0 <- (1/n)*sum((dat3$biomass - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
#confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:22))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsB.mm)[1]*pred.models$SR/(coef(nlsB.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsB.p)[1]*(pred.models$SR^(coef(nlsB.p)[2]))
pred.models$L <- coef(nlsB.linear)[1]*pred.models$SR + coef(nlsB.linear)[2]
pred.models$S <- coef(nlsB.sigm)[1]/(1+coef(nlsB.sigm)[2]*exp(-coef(nlsB.sigm)[3]*pred.models$SR))

ggfit.B <- ggplot(dat3, aes(x=as.factor(sr), y=biomass)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20), labels=c(1,10,20)) + xlab('Species Richness')+ylab('Biomass')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25)


### Metabolism ~ Richness
nlsM.mm <- nls(metab ~ a*sr/(b+sr), data=dat3, start=list(a=0.7, b=0.5), algorithm = "default", trace = T)
nlsM.p <- nls(metab ~ a*sr^b, data=dat3, start=list(a=0.5, b=0.08), algorithm = "default", trace = T)
nlsM.linear <- nls(metab ~ a*sr+b, data=dat3, start=list(a=0.1, b=0.1), algorithm = "default", trace = T)
nlsM.sigm <- nls(metab ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=0.15, b=0.5, c=0.5), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsM.sigm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$metab-mean(dat3$metab))^2)
var.resid.M0 <- (1/n)*sum((dat3$metab - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
#confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:22))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsM.mm)[1]*pred.models$SR/(coef(nlsM.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsM.p)[1]*(pred.models$SR^(coef(nlsM.p)[2]))
pred.models$L <- coef(nlsM.linear)[1]*pred.models$SR + coef(nlsM.linear)[2]
pred.models$S <- coef(nlsM.sigm)[1]/(1+coef(nlsM.sigm)[2]*exp(-coef(nlsM.sigm)[3]*pred.models$SR))

ggfit.M <- ggplot(dat3, aes(x=as.factor(sr), y=metab)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20), labels=c(1,10,20)) + xlab('Species Richness')+ylab('Metabolism')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25)


### Production ~ Richness
nlsP.mm <- nls(gainr ~ a*sr/(b+sr), data=dat3, start=list(a=2.2, b=5), algorithm = "default", trace = T)
nlsP.p <- nls(gainr ~ a*sr^b, data=dat3, start=list(a=0.2, b=0.6), algorithm = "default", trace = T)
nlsP.linear <- nls(gainr ~ a*sr+b, data=dat3, start=list(a=0.1, b=0.1), algorithm = "default", trace = T)
nlsP.sigm <- nls(gainr ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=2, b=4,c=0.25), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsP.sigm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$gainr-mean(dat3$gainr))^2)
var.resid.M0 <- (1/n)*sum((dat3$gainr - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
#confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:22))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsP.mm)[1]*pred.models$SR/(coef(nlsP.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsP.p)[1]*(pred.models$SR^(coef(nlsP.p)[2]))
pred.models$L <- coef(nlsP.linear)[1]*pred.models$SR + coef(nlsP.linear)[2]
pred.models$S <- coef(nlsP.sigm)[1]/(1+coef(nlsP.sigm)[2]*exp(-coef(nlsP.sigm)[3]*pred.models$SR))

ggfit.P <- ggplot(dat3, aes(x=as.factor(sr), y=gainr)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20), labels=c(1,10,20)) + xlab('Species Richness')+ylab('Production')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25)


### Productivity ~ Richness
nlsPB.mm <- nls(PB1 ~ a*sr/(b+sr), data=dat3, start=list(a=2, b=1), algorithm = "default", trace = T)
nlsPB.p <- nls(PB1 ~ a*sr^b, data=dat3, start=list(a=1, b=0.2), algorithm = "default", trace = T)
nlsPB.linear <- nls(PB1 ~ a*sr+b, data=dat3, start=list(a=0.5, b=1), algorithm = "default", trace = T)
nlsPB.sigm <- nls(PB1 ~ a/(1+b*exp(-c*sr)), data=dat3, start=list(a=2, b=2, c=0.25), algorithm='default', trace=T)

# calcul du R2, AIC and coefs
model <- nlsPB.sigm
n <- nrow(dat3)
var.totale <- (1/n)*sum((dat3$PB1-mean(dat3$PB1))^2)
var.resid.M0 <- (1/n)*sum((dat3$PB1 - fitted(model))^2)
var.expliq.M0 <- var.totale - var.resid.M0
R2.M0 <- var.expliq.M0/var.totale
R2.M0
AIC(model)
coef(model)
#confint(model)

# data frame wth predicted fit curves
pred.models <- data.frame(c(1:22))
names(pred.models) <- 'SR'
pred.models$MM <- coef(nlsPB.mm)[1]*pred.models$SR/(coef(nlsPB.mm)[2]+pred.models$SR)
pred.models$P <- coef(nlsPB.p)[1]*(pred.models$SR^(coef(nlsPB.p)[2]))
pred.models$L <- coef(nlsPB.linear)[1]*pred.models$SR + coef(nlsPB.linear)[2]
pred.models$S <- coef(nlsPB.sigm)[1]/(1+coef(nlsPB.sigm)[2]*exp(-coef(nlsPB.sigm)[3]*pred.models$SR))

ggfit.PB <- ggplot(dat3, aes(x=as.factor(sr), y=PB1)) + geom_point(col=adjustcolor('black', alpha=0.1), cex=2) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20), labels=c(1,10,20)) + 
  xlab('Species Richness')+ylab('Productivity')+
  geom_line(data=pred.models, aes(x=SR, y=MM), col='red', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=P), col='darkolivegreen3', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=L), col='darkgoldenrod1', lwd=1.25)+
  geom_line(data=pred.models, aes(x=SR, y=S), col='blue', lwd=1.25) 

library(egg)
windows()
ggarrange(ggfit.B, ggfit.M, ggfit.P, ggfit.PB, labels=c('A','B','C','D'), nrow=2, ncol=2)


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S4a and b
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

data <- read.csv('SimulationsSmallAlpha.csv')
datA <- data.frame()
for(i in 1:40){
  sub.dat <- subset(data, data$sr==i)
  s <- sample(c(1:nrow(sub.dat)), size=50)
  sub.sub <- sub.dat[s,]
  datA <- rbind(datA, sub.sub)
}
max.resource <- 510*101
datA$biomass <- datA$biomass/max.resource
datA$gainr <- datA$gainr/max.resource
datA$metab <- datA$metab/max.resource

datA <- data.frame(cbind(datA$biomass, datA$metab, datA$gainr, datA$PB1))
names(datA) <- c('biomass', 'metabolism','production','productivity')

library(PerformanceAnalytics)
library(corrplot)
chart.Correlation(datA, method='pearson')

### PCA on several ecosystem functions
EF <- data.frame(cbind(dat3$biomass, dat3$metab, dat3$gainr, dat3$PB1))
names(EF) <- c('Biomass','Metabolism','Production','Productivity')
res.pca = PCA(EF, scale.unit=TRUE, ncp=5, graph=T)
plot.PCA(res.pca, axes=c(1,2), choix='var')


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S4c
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

load(ss1, file='Species.Metrics1.RData')
load(ss2, file='Species.Metrics2.RData')
sss <- rbind(ss1, ss2)

### Species biomass against predation loss, consumption in color
sss$time <- 1/sss$PB
sss$ratio <- sss$PredLoss/(sss$PredLoss+sss$metab)
sss <- subset(sss, sss$time<5e10)

ggloss.sp <- ggplot(sss, aes(x=ratio, y=PB, colour=PredLoss))+geom_point(size=0.5)+
  scale_color_viridis(option='inferno')+theme_minimal()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  xlab('Predation Loss/Sum of Losses')+ylab('Productivity')+theme(panel.border = element_rect(fill=NA))+
  labs(color = "Predation Loss")

datA <- read.csv('SimulationsSmallAlpha.csv')
datA$ratio <- (datA$loss-datA$metab)/(datA$loss)
datA$PredLoss <- datA$loss-datA$metab
ggloss.fw <- ggplot(datA, aes(x=ratio, y=PB1, colour=PredLoss))+geom_point(size=0.5)+
  scale_color_viridis(option='inferno')+theme_minimal()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  xlab('Predation Loss/Sum of Losses')+ylab('Productivity')+theme(panel.border = element_rect(fill=NA))+
  labs(color = "Predation Loss")

windows(90,100)
egg::ggarrange(ggloss.sp, ggloss.fw, labels=c('a', 'b'), ncol = 1, nrow = 2)




###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S6: Explain why dominance occurs
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

load(ss1, file='Species.Metrics1.RData')
load(ss2, file='Species.Metrics2.RData')
sss <- rbind(ss1, ss2)

### Species biomass against predation loss, consumption in color
windows(90,90)
ggplot(sss, aes(x=PredLoss, y=biomass, colour=Consumption))+geom_point(size=0.5)+
  scale_color_viridis(option='inferno')+theme_minimal()+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  xlab('Predation Loss')+ylab('Biomass')+theme(panel.border = element_rect(fill=NA))


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### Appendix S7: Differences in biodiversity, network metrics, and ecosystem functioning 
### in communities with low and high home range coefficients
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

data <- read.csv('SimulationsSmallAlpha.csv')
datA <- data.frame()
for(i in 1:40){
  sub.dat <- subset(data, data$sr==i)
  s <- sample(c(1:nrow(sub.dat)), size=50)
  sub.sub <- sub.dat[s,]
  datA <- rbind(datA, sub.sub)
}
max.resource <- 510*101
datA$biomass <- datA$biomass/max.resource
datA$gainr <- datA$gainr/max.resource
datA$metab <- datA$metab/max.resource

datc <- read.csv('SimulationsHighAlpha.csv')
datC <- data.frame()
for(i in 1:40){
  sub.dat <- subset(datc, datc$sr==i)
  s <- sample(c(1:nrow(sub.dat)), size=50)
  sub.sub <- sub.dat[s,]
  datC <- rbind(datC, sub.sub)
}
max.resource <- 510*101
datC$biomass <- datC$biomass/max.resource
datC$gainr <- datC$gainr/max.resource
datC$metab <- datC$metab/max.resource


datAC <- data.frame(rbind(cbind(as.character(datA$homerange), datA$sr, datA$evesimpson, datA$biomass, datA$metab, datA$gainr, datA$PB1, datA$connectance, datA$meanTL, datA$NS01, datA$NS0001, datA$maxTL, datA$max.prop),
                          cbind(as.character(datC$homerange), datC$sr, datC$evesimpson, datC$biomass, datC$metab, datC$gainr, datC$PB1, datC$connectance, datC$meanTL, datC$NS01, datC$NS0001, datC$maxTL, datC$max.prop)))
names(datAC) <- c('homerange','sr','evenness','biomass','metab','gainr','PB','C','meanTL','NS1','NS2','maxTL', 'max.prop')
datAC$biomass <- as.numeric(as.vector(datAC$biomass))
datAC$metab <- as.numeric(as.vector(datAC$metab))
datAC$gainr <- as.numeric(as.vector(datAC$gainr))
datAC$PB <- as.numeric(as.vector(datAC$PB))
datAC$C <- as.numeric(as.vector(datAC$C))
datAC$meanTL <- as.numeric(as.vector(datAC$meanTL))
datAC$sr <- as.numeric(as.vector(datAC$sr))
datAC$NS1 <- as.numeric(as.vector(datAC$NS1))
datAC$NS2 <- as.numeric(as.vector(datAC$NS2))
datAC$maxTL <- as.numeric(as.vector(datAC$maxTL))
datAC$evenness <- as.numeric(as.vector(datAC$evenness))
datAC$max.prop <- as.numeric(as.vector(datAC$max.prop))
datAC <- datAC[order(datAC$sr),]


ggb.sr <- ggplot(datAC, aes(x=homerange, y=sr)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Species Richness') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.C <- ggplot(datAC, aes(x=homerange, y=C)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Connectance') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.TL <- ggplot(datAC, aes(x=homerange, y=meanTL)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Mean Trophic Level') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.maxTL <- ggplot(datAC, aes(x=homerange, y=maxTL)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Maximum Trophic Level') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.evenness <- ggplot(datAC, aes(x=homerange, y=evenness)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Evenness') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.prop <- ggplot(datAC, aes(x=homerange, y=max.prop)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Max. proportion') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")


windows(90,60)
egg::ggarrange(ggb.sr, ggb.C, ggb.TL, ggb.maxTL, ggb.evenness, ggb.prop,labels = c("a", "b", "c","d","e","f"),ncol = 3, nrow = 2)


datAC <- subset(datAC, datAC$sr<23)
ggb.B <- ggplot(datAC, aes(x=homerange, y=biomass)) + geom_boxplot(aes(fill=homerange), outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Biomass') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.M <- ggplot(datAC, aes(x=homerange, y=metab)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Metabolism') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.P <- ggplot(datAC, aes(x=homerange, y=gainr)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Production') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.PB <- ggplot(datAC, aes(x=homerange, y=PB)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Productivity') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.sr <- ggplot(datAC, aes(x=homerange, y=sr)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Species Richness') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.C <- ggplot(datAC, aes(x=homerange, y=C)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Connectance') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.TL <- ggplot(datAC, aes(x=homerange, y=meanTL)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Mean Trophic Level') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.prop <- ggplot(datAC, aes(x=homerange, y=max.prop)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Max. proportion') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

ggb.evenness <- ggplot(datAC, aes(x=homerange, y=evenness)) + geom_boxplot(aes(fill=homerange),outlier.size = 1) + theme_bw() + 
  xlab('') + ylab('Evenness') + theme(axis.text=element_text(size=16),axis.title=element_text(size=16)) + 
  scale_x_discrete(breaks=c('A','C'), labels=c(expression(paste(alpha, '=0.5')),expression(paste(alpha, '=2')))) +scale_fill_manual(values=c('grey80','darkgoldenrod1'))+
  theme(panel.grid.minor = element_blank())+ theme(legend.position="none")

windows(100,90)
egg::ggarrange(ggb.B, ggb.M, ggb.P, ggb.PB, ggb.C, ggb.TL, ggb.prop, ggb.evenness, labels = c("a", "b", "c", "d",'e','f','g','h'),ncol = 3, nrow = 3)








