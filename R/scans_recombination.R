############################################################################
# Recombination rate scans across venom gene regions
############################################################################

### Goal: evaluate recombination rates in venom gene regions

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & load dependencies-------------------------------

setwd('./recombination/sliding_windows/')
library(scales)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(zoo)

### Read in general data locations------------------------------------------

# Gene locations
svmp <- read.table('./resources/gene_SVMP.bed',header=F)
svsp <- read.table('./resources/gene_SVSP.bed',header=F)
pla2 <- read.table('./resources/gene_PLA2.bed',header=F)

# Region locations
svmp.reg <- read.table('./resources/region_SVMP_scaffold-mi1.bed',header=F)
svsp.reg <- read.table('./resources/region_SVSP_scaffold-mi2.bed',header=F)
pla2.reg <- read.table('./resources/region_PLA2_scaffold-mi7.bed',header=F)

### Read in recombination rates in sliding windows-------------------------

# Not pi-corrected for plotting scan purposes! Statistical tests are based
# on the pi-corrected values.

cv1.rho.100kb <- read.table('viridis.recomb.bpen10.windowed.100kb.centromereMask.txt',header=T)
cv1.rho.10kb <- read.table('viridis.recomb.bpen10.windowed.10kb.centromereMask.txt',header=T)

co1.rho.100kb <- read.table('oreganus.recomb.bpen10.windowed.100kb.centromereMask.txt',header=T)
co1.rho.10kb <- read.table('oreganus.recomb.bpen10.windowed.10kb.centromereMask.cnvMask.txt',header=T)

### Parse venom gene family chromosomes-------------------------------------

cv1.rho.100kb.mi1 <- cv1.rho.100kb[which(cv1.rho.100kb$chrom=='scaffold-mi1'),]
cv1.rho.10kb.mi1 <- cv1.rho.10kb[which(cv1.rho.10kb$chrom=='scaffold-mi1'),]
co1.rho.100kb.mi1 <- co1.rho.100kb[which(co1.rho.100kb$chrom=='scaffold-mi1'),]
co1.rho.10kb.mi1 <- co1.rho.10kb[which(co1.rho.10kb$chrom=='scaffold-mi1'),]

cv1.rho.100kb.mi2 <- cv1.rho.100kb[which(cv1.rho.100kb$chrom=='scaffold-mi2'),]
cv1.rho.10kb.mi2 <- cv1.rho.10kb[which(cv1.rho.10kb$chrom=='scaffold-mi2'),]
co1.rho.100kb.mi2 <- co1.rho.100kb[which(co1.rho.100kb$chrom=='scaffold-mi2'),]
co1.rho.10kb.mi2 <- co1.rho.10kb[which(co1.rho.10kb$chrom=='scaffold-mi2'),]

cv1.rho.100kb.mi7 <- cv1.rho.100kb[which(cv1.rho.100kb$chrom=='scaffold-mi7'),]
cv1.rho.10kb.mi7 <- cv1.rho.10kb[which(cv1.rho.10kb$chrom=='scaffold-mi7'),]
co1.rho.100kb.mi7 <- co1.rho.100kb[which(co1.rho.100kb$chrom=='scaffold-mi7'),]
co1.rho.10kb.mi7 <- co1.rho.10kb[which(co1.rho.10kb$chrom=='scaffold-mi7'),]

### Plot SVMP region--------------------------------------------------------

options(scipen=5)

par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv1.rho.10kb.mi1$start,cv1.rho.10kb.mi1$mean,pch=20,col=alpha('lightblue3',0.15),ylim=c(0,0.20),ylab='ρ',xlab='Chromosome Position (Mb)')
lines(cv1.rho.100kb.mi1$start,cv1.rho.100kb.mi1$mean)
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(co1.rho.10kb.mi1$start,co1.rho.10kb.mi1$mean,pch=20,col=alpha('lightblue3',0.15),ylim=c(0,0.25),ylab='ρ',xlab='Chromosome Position (Mb)')
lines(co1.rho.100kb.mi1$start,co1.rho.100kb.mi1$mean)
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')

lines(cv1.mi1.s,col='black')
plot(co1.rho.10kb.mi1$start,co1.rho.10kb.mi1$mean_corr,pch=20,col=alpha('lightblue3',0.15),ylim=c(0,60),ylab='Rho',xlab='Chromosome Position')
co1.mi1.s <- smooth.spline(co1.rho.100kb.mi1$start,co1.rho.100kb.mi1$mean_corr,spar=0.3)
lines(co1.mi1.s,col='black')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
# Zoomed, map points, 10 kb windows
plot(cv1.rho.10kb.mi1$start,cv1.rho.10kb.mi1$mean_corr,xlim=c(13600000,14700000),ylim=c(0,100),type='l')
plot(co1.rho.10kb.mi1$start,co1.rho.10kb.mi1$mean_corr,xlim=c(13600000,14700000),ylim=c(0,100),type='l')

plot(win.map.mi1.cv$V2,win.map.mi1.cv$V4,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.5),ylab='Rho',xlab='Chromosome Position')
lines(win.10kb.mi1.cv$start,win.10kb.mi1.cv$mean,type='l')
for (h in hot.mi1.cv$start) {
  abline(v=h,col='darkorange')
}
vert <- -0.02
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
plot(win.map.mi1.co$V2,win.map.mi1.co$V4,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.5),ylab='Rho',xlab='Chromosome Position')
lines(win.10kb.mi1.co$start,win.10kb.mi1.co$mean,type='l')
for (h in hot.mi1.co$start) {
  abline(v=h,col='darkorange')
}
vert <- -0.02
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

### Plot SVSP region--------------------------------------------------------

par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv1.rho.10kb.mi2$start,cv1.rho.10kb.mi2$mean,pch=20,col=alpha('aquamarine3',0.15),ylab='Rho',xlab='Chromosome Position',ylim=c(0,0.15))
lines(cv1.rho.100kb.mi2$start,cv1.rho.100kb.mi2$mean)
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(co1.rho.10kb.mi2$start,co1.rho.10kb.mi2$mean,pch=20,col=alpha('aquamarine3',0.15),ylab='Rho',xlab='Chromosome Position',ylim=c(0,0.3))
lines(co1.rho.100kb.mi2$start,co1.rho.100kb.mi2$mean)
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
# Zoomed, map points, 10 kb windows
plot(cv1.rho.10kb.mi2$start,cv1.rho.10kb.mi2$mean_corr,xlim=c(8350000,9200000),type='l',ylim=c(0,60))
plot(co1.rho.10kb.mi2$start,co1.rho.10kb.mi2$mean_corr,xlim=c(8350000,9200000),type='l',ylim=c(0,80))

plot(win.map.mi2.cv$V2,win.map.mi2.cv$V4,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.025,0.4),ylab='Rho',xlab='Chromosome Position')
lines(win.10kb.mi2.cv$start,win.10kb.mi2.cv$mean,type='l')
for (h in hot.mi2.cv$start) {
  abline(v=h,col='darkorange')
}
vert <- -0.02
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
plot(win.map.mi2.co$V2,win.map.mi2.co$V4,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.025,0.4),ylab='Rho',xlab='Chromosome Position')
lines(win.10kb.mi2.co$start,win.10kb.mi2.co$mean,type='l')
for (h in hot.mi2.co$start) {
  abline(v=h,col='darkorange')
}
vert <- -0.02
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}

### Plot PLA2 region--------------------------------------------------------

par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv1.rho.10kb.mi7$start,cv1.rho.10kb.mi7$mean,pch=20,col=alpha('maroon',0.15),ylab='Rho',xlab='Chromosome Position',ylim=c(0,0.15))
lines(cv1.rho.100kb.mi7$start,cv1.rho.100kb.mi7$mean)
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
plot(co1.rho.10kb.mi7$start,co1.rho.10kb.mi7$mean,pch=20,col=alpha('maroon',0.15),ylab='Rho',xlab='Chromosome Position',ylim=c(0,0.3))
lines(co1.rho.100kb.mi7$start,co1.rho.100kb.mi7$mean)
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
# Zoomed, map points, 10 kb windows

plot(win.map.mi7.cv$V2,win.map.mi7.cv$V4,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.025,0.15),ylab='Rho',xlab='Chromosome Position')
lines(win.1kb.mi7.cv$start,win.1kb.mi7.cv$mean,type='l')
for (h in hot.mi7.cv$start) {
  abline(v=h,col='darkorange')
}
vert <- -0.02
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
plot(win.map.mi7.co$V2,win.map.mi7.co$V4,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.025,0.15),ylab='Rho',xlab='Chromosome Position')
lines(win.1kb.mi7.co$start,win.1kb.mi7.co$mean,type='l')
for (h in hot.mi7.co$start) {
  abline(v=h,col='darkorange')
}
vert <- -0.02
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
