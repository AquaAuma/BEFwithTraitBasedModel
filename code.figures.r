rm(list=ls())
library(vegan)
library(stringr)
library(readr)
library(ggpubr)
library(egg)
library(ggplot2)
library(dplyr)
library(tidyr)
library(magrittr)

datall <- read.csv('SimulationsSmallAlpha.csv')

###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### PLOT FIG 1: Community Assembly
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

### Number of species for one run
dat.sr <- subset(datall, datall$run=='SP200_A11')
dat.sr$i <- seq(from=1,to=nrow(dat.sr), by=1)
plot(dat.sr$sr ~ dat.sr$i, type='l', xaxt='n', yaxt='n', xlab='', ylab='', ylim=c(0,40), lwd=4)
axis(2, at=c(10,20,30,40), labels=c(10,20,30,40), cex.axis=2.5)
abline(v=11, lty='dashed', lwd=2)
abline(v=51, lty='dashed', lwd=2)
abline(v=274, lty='dashed', lwd=2)


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### PLOT FIG 2: BEF relationships
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

dat2 <- subset(datall, datall$homerange=='A' & datall$poolsize==200)

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

ggB <- ggplot(data=dat3, aes(x=as.factor(sr),y=biomass)) + 
  geom_boxplot(outlier.size = 0, col='grey50') + geom_smooth(aes(x=sr, y=biomass), method='loess', span=1, lwd=1, color='grey50', fill=adjustcolor('grey', alpha=0.3)) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  #scale_y_discrete(breaks=c(0.1,0.2,0.3), labels=c(0,0.1,0.2,0.3,0.4)) +
  theme(panel.grid.minor = element_blank()) + xlab('')+ylab('Biomass')+
  theme(axis.text=element_text(size=16))

ggM <- ggplot(data=dat3, aes(x=as.factor(sr),y=metab)) + 
  geom_boxplot(outlier.size = 0, col='grey50') + geom_smooth(aes(x=sr, y=metab), method='loess', span=1, lwd=1, color='grey50', fill=adjustcolor('grey', alpha=0.3)) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  #scale_y_discrete(breaks=c(0,0.1,0.2,0.3,0.4), labels=c(0,0.1,0.2,0.3,0.4)) +
  theme(panel.grid.minor = element_blank()) + xlab('')+ylab('Metabolism')+
  theme(axis.text=element_text(size=16))

ggP <- ggplot(data=dat3, aes(x=as.factor(sr),y=gainr)) + 
  geom_boxplot(outlier.size = 0, col='grey50') + geom_smooth(aes(x=sr, y=gainr), method='loess', span=1, lwd=1, color='grey50', fill=adjustcolor('grey', alpha=0.3)) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  #scale_y_discrete(breaks=c(0,0.1,0.2,0.3,0.4), labels=c(0,0.1,0.2,0.3,0.4)) +
  theme(panel.grid.minor = element_blank()) + xlab('')+ylab('Production')+
  theme(axis.text=element_text(size=16))

ggPB <- ggplot(data=dat3, aes(x=as.factor(sr),y=PB1)) + 
  geom_boxplot(outlier.size = 0, col='grey50') + geom_smooth(aes(x=sr, y=PB1), method='loess', span=1, lwd=1, color='grey50', fill=adjustcolor('grey', alpha=0.3)) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  theme(panel.grid.minor = element_blank()) + xlab('')+ylab('Productivity')+
  theme(axis.text=element_text(size=16))

dat.sd <- aggregate(cbind(dat3$biomass, dat3$metab, dat3$gainr, dat3$PB1), by=list(dat3$sr), function(x) FUN=sd(x)/mean(x))
names(dat.sd) <- c('sr','biomass','metab','prod','PB')

ggB.cv <- ggplot(data=dat.sd, aes(x=sr, y=biomass)) + geom_line(col='grey50',size=1.5) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  theme(panel.grid.minor = element_blank())+xlab('Species Richness')+ylab('CV Biomass')+ ylim(0,0.8)+
  theme(axis.text=element_text(size=16))

ggM.cv <- ggplot(data=dat.sd, aes(x=sr, y=metab)) + geom_line(col='grey50',size=1.5) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  theme(panel.grid.minor = element_blank())+xlab('Species Richness')+ylab('CV Metabolism')+ ylim(0,0.8)+
  theme(axis.text=element_text(size=16))

ggP.cv <- ggplot(data=dat.sd, aes(x=sr, y=prod)) + geom_line(col='grey50',size=1.5) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  theme(panel.grid.minor = element_blank())+xlab('Species Richness')+ylab('CV Production')+ ylim(0,0.8)+
  theme(axis.text=element_text(size=16))

ggPB.cv <- ggplot(data=dat.sd, aes(x=sr, y=PB)) + geom_line(col='grey50',size=1.5) +
  theme_bw() + theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  theme(panel.grid.minor = element_blank())+xlab('Species Richness')+ylab('CV Productivity')+ ylim(0,0.8)+
  theme(axis.text=element_text(size=16))

windows(120,60)
egg::ggarrange(ggB, ggM, ggP, ggPB, ggB.cv, ggM.cv, ggP.cv, ggPB.cv, labels = c("a", "b", "c","d","e", "f", "g", "h"),ncol = 4, nrow = 2)



###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### PLOT FIG 3: Variance Partitioning
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
dat3$B.dom <- dat3$B.dom/max.resource


### Explain function 1: Biomass
###-----------------------------------------------------------------------------------

# Exploration
windows()
par(mfrow=c(3,3))
plot(biomass ~ sr, data=dat3)
plot(biomass ~ evesimpson, data=dat3)
plot(biomass ~ evesimpson.ab, data=dat3)
plot(biomass ~ maxTL, data=dat3)
plot(biomass ~ B.dom, data=dat3)
plot(biomass ~ top, data=dat3)
plot(biomass ~ NS00001, data=dat3)
plot(biomass ~ max.prop, data=dat3)
plot(biomass ~ connectance, data=dat3)

# Variance partitioning
dat4 <- subset(dat3, sr>1)
dat4$sr2 <- (dat4$sr-mean(dat4$sr))^2
dat4$NS2 <- (dat4$NS00001-mean(dat4$NS00001))^2

Dominance <- data.frame(dat4$evesimpson, dat4$B.dom, dat4$max.prop)
Horizon <- data.frame(dat4$sr, dat4$sr2)
Vertical <- data.frame(dat4$maxTL, dat4$meanTL, dat4$connectance, dat4$top)
var.modB <- varpart(dat4$biomass, Dominance, Horizon, Vertical, dat4$NS00001)


### Explain function 2: Metabolism
###-----------------------------------------------------------------------------------

# Exploration
windows()
par(mfrow=c(3,3))
plot(metab ~ sr, data=dat3)
plot(metab ~ evesimpson, data=dat3)
plot(metab ~ evesimpson.ab, data=dat3)
plot(metab ~ maxTL, data=dat3)
plot(metab ~ B.dom, data=dat3)
plot(metab ~ top, data=dat3)
plot(metab ~ NS00001, data=dat3)
plot(metab ~ max.prop, data=dat3)
plot(metab ~ connectance, data=dat3)

# Variance partitioning
dat4 <- subset(dat3, sr>1)
dat4$sr2 <- (dat4$sr-mean(dat4$sr))^2
dat4$NS2 <- (dat4$NS00001-mean(dat4$NS00001))^2
dat4$top2 <- (dat4$top - mean(dat4$top))^2

Dominance <- data.frame(dat4$evesimpson, dat4$B.dom, dat4$max.prop)
Horizon <- data.frame(dat4$sr, dat4$sr2)
Vertical <- data.frame(dat4$maxTL, dat4$meanTL, dat4$connectance, dat4$top, dat4$top2)
var.modM <- varpart(dat4$metab, Dominance, Horizon, Vertical, dat4$NS00001)


### Explain function 3: Production
###-----------------------------------------------------------------------------------

# Exploration
windows()
par(mfrow=c(3,3))
plot(gainr ~ sr, data=dat3)
plot(gainr ~ evesimpson, data=dat3)
plot(gainr ~ evesimpson.ab, data=dat3)
plot(gainr ~ maxTL, data=dat3)
plot(gainr ~ B.dom, data=dat3)
plot(gainr ~ top, data=dat3)
plot(gainr ~ NS00001, data=dat3)
plot(gainr ~ max.prop, data=dat3)
plot(gainr ~ connectance, data=dat3)

# Variance partitioning
dat4 <- subset(dat3, sr>1)
dat4$sr2 <- (dat4$sr-mean(dat4$sr))^2
dat4$NS2 <- (dat4$NS00001-mean(dat4$NS00001))^2
dat4$max.prop2 <- (dat4$max.prop-mean(dat4$max.prop))^2

Dominance <- data.frame(dat4$evesimpson, dat4$B.dom, dat4$max.prop,dat4$max.prop2)
Horizon <- data.frame(dat4$sr, dat4$sr2)
Vertical <- data.frame(dat4$maxTL, dat4$meanTL, dat4$connectance, dat4$top)
var.modP <- varpart(dat4$gainr, Dominance, Horizon, Vertical, dat4$NS00001)


### Explain function 4: Productivity
###-----------------------------------------------------------------------------------

# Exploration
windows()
par(mfrow=c(3,3))
plot(PB1 ~ sr, data=dat3)
plot(PB1 ~ evesimpson, data=dat3)
plot(PB1 ~ evesimpson.ab, data=dat3)
plot(PB1 ~ maxTL, data=dat3)
plot(PB1 ~ B.dom, data=dat3)
plot(PB1 ~ top, data=dat3)
plot(PB1 ~ NS00001, data=dat3)
plot(PB1 ~ max.prop, data=dat3)
plot(PB1 ~ connectance, data=dat3)

# Variance partitioning
dat4 <- subset(dat3, sr>1)
dat4$sr2 <- (dat4$sr-mean(dat4$sr))^2
dat4$NS2 <- (dat4$NS00001-mean(dat4$NS00001))^2
dat4$max.prop2 <- (dat4$max.prop-mean(dat4$max.prop))^2

Dominance <- data.frame(dat4$evesimpson, dat4$B.dom, dat4$max.prop, dat4$max.prop2)
Horizon <- data.frame(dat4$sr, dat4$sr2)
Vertical <- data.frame(dat4$maxTL, dat4$meanTL, dat4$connectance, dat4$top)
var.modPB <- varpart(dat4$PB1, Dominance, Horizon, Vertical, Trait)

# Plot all models
windows(80,50)
par(mfrow=c(2,2), mar=c(2,2,2,2))
plot(var.modB, Xnames=c('','','',''), bg=c('red','yellow','green','blue'), cutoff=0.01, digits=1)
plot(var.modM, Xnames=c('','','',''), bg=c('red','yellow','green','blue'), cutoff=0.01, digits=1)
plot(var.modP, Xnames=c('','','',''), bg=c('red','yellow','green','blue'), cutoff=0.01, digits=1)
plot(var.modPB, Xnames=c('','','',''), bg=c('red','yellow','green','blue'), cutoff=0.01, digits=1)


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### PLOT FIG 4: Dominance and structure
###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------

windows(120,120)
dat3$red <- 'black'
dat3[dat3$sr==10,]$red <- 'grey60'
dat3[dat3$sr %in% c(5,30),]$red <- 'grey30'
ggB <- ggplot(data=dat3, aes(x=as.factor(sr),y=biomass, color=red, fill=red)) + 
  geom_boxplot(outlier.size = 1) + 
  theme_bw() + theme(axis.text=element_text(size=30),axis.title=element_text(size=25))+
  scale_x_discrete(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40)) + 
  theme(panel.grid.minor = element_blank()) + xlab('')+ylab('')+
  scale_color_manual(values=c('black','black','black')) + 
  geom_smooth(aes(x=sr, y=biomass), method='loess', span=0.5, lwd=1, color='grey50', fill=adjustcolor('grey', alpha=0.3))+
  scale_fill_manual(values=c("white","grey80","grey50")) +
  theme(legend.position="none")+theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA)) # bg of the plot
ggB


###-----------------------------------------------------------------------------------
###-----------------------------------------------------------------------------------
### PLOT FIG 5: generalists and specialists communities
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
max.resource <- 101*510
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
datC$biomass <- datC$biomass/max.resource
datC$gainr <- datC$gainr/max.resource
datC$metab <- datC$metab/max.resource


datAC <- data.frame(rbind(cbind(as.character(datA$homerange), datA$sr, datA$biomass, datA$metab, datA$gainr, datA$PB1, datA$NS00001, datA$propR.used),
                          cbind(as.character(datC$homerange), datC$sr, datC$biomass, datC$metab, datC$gainr, datC$PB1, datC$NS00001, datC$propR.used)))
names(datAC) <- c('homerange','sr','biomass','metab','gainr','PB','NS','propR')
datAC$biomass <- as.numeric(as.vector(datAC$biomass))
datAC$metab <- as.numeric(as.vector(datAC$metab))
datAC$gainr <- as.numeric(as.vector(datAC$gainr))
datAC$PB <- as.numeric(as.vector(datAC$PB))
datAC$sr <- as.numeric(as.vector(datAC$sr))
datAC$NS <- as.numeric(as.vector(datAC$NS))
datAC$propR <- as.numeric(as.vector(datAC$propR))
datAC <- datAC[order(datAC$sr),]


datAC.stats <- datAC %>%
  group_by(sr, homerange) %>%
  summarise(n = n(),mean = mean(biomass),median = median(biomass),sd = sd(biomass)) %>%
  mutate(sem = sd / sqrt(n - 1),CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggAC.B <- ggplot(datAC.stats, aes(x=sr, y=mean, color = homerange)) +
  geom_point(data=datA, aes(x=sr, y=biomass), shape=1, col=adjustcolor('grey60', alpha=0.4)) + 
  geom_point(data=datC, aes(x=sr, y=biomass), shape=1, col=adjustcolor('darkgoldenrod1', alpha=0.4))+
  geom_line(aes(x=sr, y=mean, color=homerange)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=homerange),alpha=0.6)+ theme_bw()+
  scale_color_manual(values=c('grey60','darkgoldenrod1')) + 
  scale_fill_manual(values=c('grey60','darkgoldenrod1')) +
  theme(panel.grid.minor = element_blank())+xlab('')+ylab('Biomass')+
  theme(axis.text.x=element_blank(),axis.text=element_text(size=16),axis.title=element_text(size=16))+ 
  theme(legend.position="none")

datAC.stats <- datAC %>%
  group_by(sr, homerange) %>%
  summarise(n = n(),mean = mean(metab),median = median(metab),sd = sd(metab)) %>%
  mutate(sem = sd / sqrt(n - 1),CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggAC.M <- ggplot(datAC.stats, aes(x=sr, y=mean, color = homerange)) +
  geom_point(data=datA, aes(x=sr, y=metab), shape=1, col=adjustcolor('grey60', alpha=0.4)) + 
  geom_point(data=datC, aes(x=sr, y=metab), shape=1, col=adjustcolor('darkgoldenrod1', alpha=0.4))+
  geom_line(aes(x=sr, y=mean, color=homerange)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=homerange),alpha=0.6)+ theme_bw()+
  scale_color_manual(values=c('grey60','darkgoldenrod1')) + scale_fill_manual(values=c('grey60','darkgoldenrod1')) +
  theme(panel.grid.minor = element_blank())+xlab('')+ylab('Metabolism')+
  theme(axis.text.x=element_blank(),axis.text=element_text(size=16),axis.title=element_text(size=16))+ theme(legend.position="none")

datAC.stats <- datAC %>%
  group_by(sr, homerange) %>%
  summarise(n = n(),mean = mean(gainr),median = median(gainr),sd = sd(gainr)) %>%
  mutate(sem = sd / sqrt(n - 1),CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggAC.P <- ggplot(datAC.stats, aes(x=sr, y=mean, color = homerange)) +
  geom_point(data=datA, aes(x=sr, y=gainr), shape=1, col=adjustcolor('grey60', alpha=0.4)) + 
  geom_point(data=datC, aes(x=sr, y=gainr), shape=1, col=adjustcolor('darkgoldenrod1', alpha=0.4))+
  geom_line(aes(x=sr, y=mean, color=homerange)) +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=homerange),alpha=0.6)+ theme_bw()+
  scale_color_manual(values=c('grey60','darkgoldenrod1')) + scale_fill_manual(values=c('grey60','darkgoldenrod1')) +
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40))+
  theme(panel.grid.minor = element_blank())+xlab('')+ylab('Production')+
  theme(axis.text.x=element_blank(), axis.text=element_text(size=16),axis.title=element_text(size=16))+ theme(legend.position="none")

datAC.stats <- datAC %>%
  group_by(sr, homerange) %>%
  summarise(n = n(),mean = mean(PB),median = median(PB),sd = sd(PB)) %>%
  mutate(sem = sd / sqrt(n - 1),CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggAC.PB <- ggplot(datAC.stats, aes(x=sr, y=mean, color = homerange)) +
  geom_point(data=datA, aes(x=sr, y=PB1), shape=1, col=adjustcolor('grey60', alpha=0.4)) + 
  geom_point(data=datC, aes(x=sr, y=PB1), shape=1, col=adjustcolor('darkgoldenrod1', alpha=0.4)) +
  geom_line(aes(x=sr, y=mean, color=homerange))+
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=homerange),alpha=0.6)+ theme_bw()+
  scale_color_manual(values=c('grey60','darkgoldenrod1')) + 
  scale_fill_manual(values=c('grey60','darkgoldenrod1')) +
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40))+
  theme(panel.grid.minor = element_blank())+xlab('')+ylab('Productivity')+
  theme(axis.text.x=element_blank(), axis.text=element_text(size=16),axis.title=element_text(size=16))+ theme(legend.position="none")

datAC.stats <- datAC %>%
  group_by(sr, homerange) %>%
  summarise(n = n(),mean = mean(NS),median = median(NS),sd = sd(NS)) %>%
  mutate(sem = sd / sqrt(n - 1),CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggAC.NS <- ggplot(datAC.stats, aes(x=sr, y=mean, color = homerange)) +
  geom_point(data=datA, aes(x=sr, y=NS00001), shape=1, col=adjustcolor('grey60', alpha=0.4)) + 
  geom_point(data=datC, aes(x=sr, y=NS00001), shape=1, col=adjustcolor('darkgoldenrod1', alpha=0.4)) +
  geom_line(aes(x=sr, y=mean, color=homerange))+
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=homerange),alpha=0.6)+ theme_bw()+
  scale_color_manual(values=c('grey60','darkgoldenrod1')) + 
  scale_fill_manual(values=c('grey60','darkgoldenrod1')) +
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1), labels=c(0,0.5,1))+
  theme(panel.grid.minor = element_blank())+xlab('Species Richness')+ylab('% of Trait Space filled')+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+ theme(legend.position="none")

datAC.stats <- datAC %>%
  group_by(sr, homerange) %>%
  summarise(n = n(),mean = mean(propR),median = median(propR),sd = sd(propR)) %>%
  mutate(sem = sd / sqrt(n - 1),CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggAC.propR <- ggplot(datAC.stats, aes(x=sr, y=mean, color = homerange)) +
  geom_point(data=datA, aes(x=sr, y=propR.used), shape=1, col=adjustcolor('grey60', alpha=0.4)) + 
  geom_point(data=datC, aes(x=sr, y=propR.used), shape=1, col=adjustcolor('darkgoldenrod1', alpha=0.4)) +
  geom_line(aes(x=sr, y=mean, color=homerange))+
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=homerange),alpha=0.6)+ theme_bw()+
  scale_color_manual(values=c('grey60','darkgoldenrod1')) + 
  scale_fill_manual(values=c('grey60','darkgoldenrod1')) +
  scale_x_continuous(breaks=c(1,10,20,30,40), labels=c(1,10,20,30,40))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1), labels=c(0,0.5,1))+
  theme(panel.grid.minor = element_blank())+xlab('Species Richness')+ylab('% of Resource used')+
  theme(axis.text=element_text(size=16),axis.title=element_text(size=16))+ theme(legend.position="none")


windows(100,100)
egg::ggarrange(ggAC.B, ggAC.M, ggAC.P, ggAC.PB, ggAC.NS, ggAC.propR, labels = c('a','b','c','d','e','f'),nrow=3, ncol=2)

