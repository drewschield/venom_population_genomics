############################################################################
# Compare mean values from 'other' venom genes to genome background
############################################################################

### Goal: Look at mean estimates for 'other' venom genes stacked against the
### background distribution across autosomes.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Load dependencies-------------------------------------------------------

library(scales)
library(data.table)
library(dplyr)

options('stringsAsFactors'=FALSE)

### Read in venom gene means------------------------------------------------

cv.other.taj <- read.table('./otherVG.TajimaD_mean.cv1.txt',header=T)
co.other.taj <- read.table('./otherVG.TajimaD_mean.co1.txt',header=T)
cv.other.ihs <- read.table('./otherVG.ihs_mean.cv1.txt',header=T)
co.other.ihs <- read.table('./otherVG.ihs_mean.co1.txt',header=T)
cv.other.beta <- read.table('./otherVG/otherVG.beta_mean.cv1.txt',header=T)
co.other.beta <- read.table('./otherVG/otherVG.beta_mean.co1.txt',header=T)

### Read in genome-wide data------------------------------------------------

cv.genome.taj <- read.table('./tajima_d/cv.colorado.all.10kb.Tajima.D',header=T)
co.genome.taj <- read.table('./tajima_d/co.california.all.10kb.Tajima.D',header=T)
cv.genome.ihs <- read.table('./rehh/cv.all_ihs.10kb.txt',header=T)
co.genome.ihs <- read.table('./rehh/co.all_ihs.10kb.txt',header=T)
cv.genome.beta <- read.table('./beta/results/cv1.phased.all.betascores.10kb.txt',header=T)
co.genome.beta <- read.table('./beta/results/co1.phased.all.betascores.10kb.txt',header=T)

### Calculate genome-wide quantiles-----------------------------------------

quantile(cv.genome.taj$TajimaD,c(0.025,0.975),na.rm=T)
quantile(co.genome.taj$TajimaD,c(0.025,0.975),na.rm=T)
quantile(as.numeric(cv.genome.ihs$iHS),c(0.95),na.rm=T)
quantile(as.numeric(co.genome.ihs$iHS),c(0.95),na.rm=T)
quantile(as.numeric(cv.genome.beta$Beta1.),c(0.95),na.rm=T)
quantile(as.numeric(co.genome.beta$Beta1.),c(0.95),na.rm=T)

### Plot density distributions----------------------------------------------

par(mfrow=c(3,2))
plot(density(cv.genome.taj$TajimaD,na.rm=T),main=NA,xlab="Tajima's D",ylim=c(-0.05,0.5))
polygon(density(cv.genome.taj$TajimaD,na.rm=T),main=NA,col='lightgrey')
for (ven in cv.other.taj$TajimaD) {
  abline(v=ven,col='darkblue')
}
plot(density(co.genome.taj$TajimaD,na.rm=T),main=NA,xlab="Tajima's D",ylim=c(-0.05,0.5))
polygon(density(co.genome.taj$TajimaD,na.rm=T),main=NA,col='lightgrey')
for (ven in co.other.taj$TajimaD) {
  abline(v=ven,col='darkblue')
}
plot(density(as.numeric(cv.genome.ihs$iHS),na.rm=T),main=NA,xlab="iHS",ylim=c(-0.1,0.95))
polygon(density(as.numeric(cv.genome.ihs$iHS),na.rm=T),main=NA,col='lightgrey')
for (ven in cv.other.ihs$iHS) {
  abline(v=ven,col='darkblue')
}
plot(density(as.numeric(co.genome.ihs$iHS),na.rm=T),main=NA,xlab="iHS",ylim=c(-0.1,0.95))
polygon(density(as.numeric(co.genome.ihs$iHS),na.rm=T),main=NA,col='lightgrey')
for (ven in co.other.ihs$iHS) {
  abline(v=ven,col='darkblue')
}
plot(density(as.numeric(cv.genome.beta$Beta1.),na.rm=T),main=NA,xlab="beta",xlim=c(0,15),ylim=c(-0.2,2))
polygon(density(as.numeric(cv.genome.beta$Beta1.),na.rm=T),main=NA,col='lightgrey')
for (ven in cv.other.beta$beta) {
  abline(v=ven,col='darkblue')
}
plot(density(as.numeric(co.genome.beta$Beta1.),na.rm=T),main=NA,xlab="beta",xlim=c(0,15),ylim=c(-0.2,2))
polygon(density(as.numeric(co.genome.beta$Beta1.),na.rm=T),main=NA,col='lightgrey')
for (ven in co.other.beta$beta) {
  abline(v=ven,col='darkblue')
}
