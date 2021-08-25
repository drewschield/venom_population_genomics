############################################################################
# Recombination rate distributions for venom gene regions
############################################################################

### Goal: evaluate recombination rates in venom gene regions and compare to
### chromosome/genomic background and non-venom paralogs.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & load dependencies-------------------------------

setwd('./recombination/')
library(scales)
library(data.table)
library(dplyr)

options('stringsAsFactors'=FALSE)

### Read in general coordinates---------------------------------------------

# Gene locations
svmp <- read.table('./resources/gene_SVMP.bed',header=F)
svsp <- read.table('./resources/gene_SVSP.bed',header=F)
pla2 <- read.table('./resources/gene_PLA2.bed',header=F)

# Region locations
svmp.reg <- read.table('./resources/region_SVMP_scaffold-mi1.bed',header=F)
svsp.reg <- read.table('./resources/region_SVSP_scaffold-mi2.bed',header=F)
pla2.reg <- read.table('./resources/region_PLA2_scaffold-mi7.bed',header=F)

### ------------------------------------------------------------------------
### Compare recombination rates in venom regions versus backgrounds
### ------------------------------------------------------------------------

# Update 06.08.21 - these are pi-corrected rates

### Read in data------------------------------------------------------------

r.cv <- read.table('./sliding_windows/viridis.recomb.bpen10.windowed.10kb.centromereMask.piCorrected.txt',header=T)
r.co <- read.table('./sliding_windows/oreganus.recomb.bpen10.windowed.10kb.centromereMask.cnvMask.piCorrected.txt',header=T)

r.mi1.cv <- r.cv[r.cv$chrom=='scaffold-mi1',]
r.mi2.cv <- r.cv[r.cv$chrom=='scaffold-mi2',]
r.mi7.cv <- r.cv[r.cv$chrom=='scaffold-mi7',]
r.mi1.co <- r.co[r.co$chrom=='scaffold-mi1',]
r.mi2.co <- r.co[r.co$chrom=='scaffold-mi2',]
r.mi7.co <- r.co[r.co$chrom=='scaffold-mi7',]

r.svmp.cv <- r.mi1.cv[which(r.mi1.cv$start>=svmp.reg$V2 & r.mi1.cv$start<=svmp.reg$V3),]
r.svsp.cv <- r.mi2.cv[which(r.mi2.cv$start>=svsp.reg$V2 & r.mi2.cv$start<=svsp.reg$V3),]
r.pla2.cv <- r.mi7.cv[which(r.mi7.cv$start>=pla2.reg$V2 & r.mi7.cv$start<=pla2.reg$V3),]
r.svmp.co <- r.mi1.co[which(r.mi1.co$start>=svmp.reg$V2 & r.mi1.co$start<=svmp.reg$V3),]
r.svsp.co <- r.mi2.co[which(r.mi2.co$start>=svsp.reg$V2 & r.mi2.co$start<=svsp.reg$V3),]
r.pla2.co <- r.mi7.co[which(r.mi7.co$start>=pla2.reg$V2 & r.mi7.co$start<=pla2.reg$V3),]

r.mi1.bk.cv <- r.mi1.cv[which(r.mi1.cv$start<svmp.reg$V2 | r.mi1.cv$start>svmp.reg$V3),]
r.mi2.bk.cv <- r.mi2.cv[which(r.mi2.cv$start<svsp.reg$V2 | r.mi2.cv$start>svsp.reg$V3),]
r.mi7.bk.cv <- r.mi7.cv[which(r.mi7.cv$start<pla2.reg$V2 | r.mi7.cv$start>pla2.reg$V3),]
r.mi1.bk.co <- r.mi1.co[which(r.mi1.co$start<svmp.reg$V2 | r.mi1.co$start>svmp.reg$V3),]
r.mi2.bk.co <- r.mi2.co[which(r.mi2.co$start<svsp.reg$V2 | r.mi2.co$start>svsp.reg$V3),]
r.mi7.bk.co <- r.mi7.co[which(r.mi7.co$start<pla2.reg$V2 | r.mi7.co$start>pla2.reg$V3),]

r.nv.svmp.cv <- read.table('./sliding_windows/non-venom/viridis.recomb.bpen10.windowed.10kb.non-venom_SVMP.piCorrected.txt',header=T)
r.nv.svsp.cv <- read.table('./sliding_windows/non-venom/viridis.recomb.bpen10.windowed.10kb.non-venom_SVSP.piCorrected.txt',header=T)
r.nv.pla2.cv <- read.table('./sliding_windows/non-venom/viridis.recomb.bpen10.windowed.10kb.non-venom_PLA2.piCorrected.txt',header=T)
r.nv.svmp.co <- read.table('./sliding_windows/non-venom/oreganus.recomb.bpen10.windowed.10kb.non-venom_SVMP.piCorrected.txt',header=T)
r.nv.svsp.co <- read.table('./sliding_windows/non-venom/oreganus.recomb.bpen10.windowed.10kb.non-venom_SVSP.piCorrected.txt',header=T)
r.nv.pla2.co <- read.table('./sliding_windows/non-venom/oreganus.recomb.bpen10.windowed.10kb.non-venom_PLA2.piCorrected.txt',header=T)

### Plot distributions------------------------------------------------------

par(mfrow=c(2,3))
boxplot(r.mi1.bk.cv$mean_corr,r.nv.svmp.cv$mean,r.svmp.cv$mean_corr,outline=F,ylab='ρ/π',main='SVMP',names=c('Chr. 9','NV','SVMP'),col='lightblue3')
boxplot(r.mi2.bk.cv$mean_corr,r.nv.svsp.cv$mean,r.svsp.cv$mean_corr,outline=F,ylab='ρ/π',main='SVSP',names=c('Chr. 10','NV','SVSP'),col='aquamarine3')
boxplot(r.mi7.bk.cv$mean_corr,r.nv.pla2.cv$mean,r.pla2.cv$mean_corr,outline=F,ylab='ρ/π',main='PLA2',names=c('Chr. 15','NV','PLA2'),col='maroon')

boxplot(r.mi1.bk.co$mean_corr,r.nv.svmp.co$mean,r.svmp.co$mean_corr,outline=F,ylab='ρ/π',main='SVMP',names=c('Chr. 9','NV','SVMP'),col='lightblue3')
boxplot(r.mi2.bk.co$mean_corr,r.nv.svsp.co$mean,r.svsp.co$mean_corr,outline=F,ylab='ρ/π',main='SVSP',names=c('Chr. 10','NV','SVSP'),col='aquamarine3')
boxplot(r.mi7.bk.co$mean_corr,r.nv.pla2.co$mean,r.pla2.co$mean_corr,outline=F,ylab='ρ/π',main='PLA2',names=c('Chr. 15','NV','PLA2'),col='maroon')

### Summary statistics (mean & SD)------------------------------------------

mean(r.svmp.cv$mean_corr,na.rm=T)
mean(r.svmp.co$mean_corr,na.rm=T)
mean(r.svsp.cv$mean_corr,na.rm=T)
mean(r.svsp.co$mean_corr,na.rm=T)
mean(r.pla2.cv$mean_corr,na.rm=T)
mean(r.pla2.co$mean_corr,na.rm=T)

sd(r.svmp.cv$mean_corr,na.rm=T)
sd(r.svmp.co$mean_corr,na.rm=T)
sd(r.svsp.cv$mean_corr,na.rm=T)
sd(r.svsp.co$mean_corr,na.rm=T)
sd(r.pla2.cv$mean_corr,na.rm=T)
sd(r.pla2.co$mean_corr,na.rm=T)

mean(r.mi1.bk.cv$mean_corr,na.rm=T)
mean(r.mi1.bk.co$mean_corr,na.rm=T)
mean(r.mi2.bk.cv$mean_corr,na.rm=T)
mean(r.mi2.bk.co$mean_corr,na.rm=T)
mean(r.mi7.bk.cv$mean_corr,na.rm=T)
mean(r.mi7.bk.co$mean_corr,na.rm=T)

sd(r.mi1.bk.cv$mean_corr,na.rm=T)
sd(r.mi1.bk.co$mean_corr,na.rm=T)
sd(r.mi2.bk.cv$mean_corr,na.rm=T)
sd(r.mi2.bk.co$mean_corr,na.rm=T)
sd(r.mi7.bk.cv$mean_corr,na.rm=T)
sd(r.mi7.bk.co$mean_corr,na.rm=T)

mean(r.nv.svmp.cv$mean,na.rm=T)
mean(r.nv.svmp.co$mean,na.rm=T)
mean(r.nv.svsp.cv$mean,na.rm=T)
mean(r.nv.svsp.co$mean,na.rm=T)
mean(r.nv.pla2.cv$mean,na.rm=T)
mean(r.nv.pla2.co$mean,na.rm=T)

sd(r.nv.svmp.cv$mean,na.rm=T)
sd(r.nv.svmp.co$mean,na.rm=T)
sd(r.nv.svsp.cv$mean,na.rm=T)
sd(r.nv.svsp.co$mean,na.rm=T)
sd(r.nv.pla2.cv$mean,na.rm=T)
sd(r.nv.pla2.co$mean,na.rm=T)

### Mann-Whitney U tests comparing venom regions to backgrounds-------------

# C. viridis
wilcox.test(r.mi1.bk.cv$mean_corr,r.svmp.cv$mean_corr)
wilcox.test(r.mi2.bk.cv$mean_corr,r.svsp.cv$mean_corr)
wilcox.test(r.mi7.bk.cv$mean_corr,r.pla2.cv$mean_corr)

wilcox.test(r.nv.svmp.cv$mean,r.svmp.cv$mean_corr)
wilcox.test(r.nv.svsp.cv$mean,r.svsp.cv$mean_corr)
wilcox.test(r.nv.pla2.cv$mean,r.pla2.cv$mean_corr)

# C. oreganus
wilcox.test(r.mi1.bk.co$mean_corr,r.svmp.co$mean_corr)
wilcox.test(r.mi2.bk.co$mean_corr,r.svsp.co$mean_corr)
wilcox.test(r.mi7.bk.co$mean_corr,r.pla2.co$mean_corr)

wilcox.test(r.nv.svmp.co$mean,r.svmp.co$mean_corr)
wilcox.test(r.nv.svsp.co$mean,r.svsp.co$mean_corr)
wilcox.test(r.nv.pla2.co$mean,r.pla2.co$mean_corr)

### Mann-Whitney U tests comparing venom regions to non-venom paralogs------

# C. viridis
wilcox.test(r.para.svmp.cv$V4,r.svmp.cv$V4)
wilcox.test(r.para.svsp.cv$V4,r.svsp.cv$V4)
wilcox.test(r.para.pla2.cv$V4,r.pla2.cv$V4)

# C. oreganus
wilcox.test(r.para.svmp.co$V4,r.svmp.co$V4)
wilcox.test(r.para.svsp.co$V4,r.svsp.co$V4)
wilcox.test(r.para.pla2.co$V4,r.pla2.co$V4)


### ------------------------------------------------------------------------
### Plot venom gene regions with recombination rates, genes, and hotspots
### ------------------------------------------------------------------------

### Read in data------------------------------------------------------------

# 100 kb windows
win.100kb.cv <- read.table('sliding_windows/viridis.recomb.bpen10.windowed.100kb.centromereMask.txt',header=T)
win.100kb.mi1.cv <- win.100kb.cv[which(win.100kb.cv$chrom=='scaffold-mi1'),]
win.100kb.mi2.cv <- win.100kb.cv[which(win.100kb.cv$chrom=='scaffold-mi2'),]
win.100kb.mi7.cv <- win.100kb.cv[which(win.100kb.cv$chrom=='scaffold-mi7'),]

win.100kb.co <- read.table('sliding_windows/oreganus.recomb.bpen10.windowed.100kb.centromereMask.txt',header=T)
win.100kb.mi1.co <- win.100kb.co[which(win.100kb.co$chrom=='scaffold-mi1'),]
win.100kb.mi2.co <- win.100kb.co[which(win.100kb.co$chrom=='scaffold-mi2'),]
win.100kb.mi7.co <- win.100kb.co[which(win.100kb.co$chrom=='scaffold-mi7'),]

# 10 kb windows
win.10kb.cv <- read.table('sliding_windows/viridis.recomb.bpen10.windowed.10kb.centromereMask.txt',header=T)
win.10kb.mi1.cv <- win.10kb.cv[which(win.10kb.cv$chrom=='scaffold-mi1'),]
win.10kb.mi2.cv <- win.10kb.cv[which(win.10kb.cv$chrom=='scaffold-mi2'),]
win.10kb.mi7.cv <- win.10kb.cv[which(win.10kb.cv$chrom=='scaffold-mi7'),]

win.10kb.co <- read.table('sliding_windows/oreganus.recomb.bpen10.windowed.10kb.centromereMask.cnvMask.txt',header=T)
win.10kb.mi1.co <- win.10kb.co[which(win.10kb.co$chrom=='scaffold-mi1'),]
win.10kb.mi2.co <- win.10kb.co[which(win.10kb.co$chrom=='scaffold-mi2'),]
win.10kb.mi7.co <- win.10kb.co[which(win.10kb.co$chrom=='scaffold-mi7'),]

# 1 kb windows
win.1kb.cv <- read.table('sliding_windows/viridis.recomb.bpen10.windowed.1kb.centromereMask.txt',header=T)
win.1kb.mi1.cv <- win.1kb.cv[which(win.1kb.cv$chrom=='scaffold-mi1'),]
win.1kb.mi2.cv <- win.1kb.cv[which(win.1kb.cv$chrom=='scaffold-mi2'),]
win.1kb.mi7.cv <- win.1kb.cv[which(win.1kb.cv$chrom=='scaffold-mi7'),]

win.1kb.co <- read.table('sliding_windows/oreganus.recomb.bpen10.windowed.1kb.centromereMask.cnvMask.txt',header=T)
win.1kb.mi1.co <- win.1kb.co[which(win.1kb.co$chrom=='scaffold-mi1'),]
win.1kb.mi2.co <- win.1kb.co[which(win.1kb.co$chrom=='scaffold-mi2'),]
win.1kb.mi7.co <- win.1kb.co[which(win.1kb.co$chrom=='scaffold-mi7'),]

# Recombination map intervals
win.map.mi1.cv <- read.table('recombination_maps/viridis.recomb.bpen10.scaffold-mi1.txt',header=F)
win.map.mi2.cv <- read.table('recombination_maps/viridis.recomb.bpen10.scaffold-mi2.txt',header=F)
win.map.mi7.cv <- read.table('recombination_maps/viridis.recomb.bpen10.scaffold-mi7.txt',header=F)

win.map.mi1.co <- read.table('recombination_maps/oreganus.recomb.bpen10.cnv_mask.scaffold-mi1.txt',header=F)
win.map.mi2.co <- read.table('recombination_maps/oreganus.recomb.bpen10.cnv_mask.scaffold-mi2.txt',header=F)
win.map.mi7.co <- read.table('recombination_maps/oreganus.recomb.bpen10.cnv_mask.scaffold-mi7.txt',header=F)

# Hotspot locations
hot.cv <- read.table('hotspots/hotspots.filtered.viridis.txt',header=T)
hot.mi1.cv <- hot.cv[which(hot.cv$chrom=='scaffold-mi1'),]
hot.mi2.cv <- hot.cv[which(hot.cv$chrom=='scaffold-mi2'),]
hot.mi7.cv <- hot.cv[which(hot.cv$chrom=='scaffold-mi7'),]

hot.co <- read.table('hotspots/hotspots.filtered.oreganus.txt',header=T)
hot.mi1.co <- hot.co[which(hot.co$chrom=='scaffold-mi1'),]
hot.mi2.co <- hot.co[which(hot.co$chrom=='scaffold-mi2'),]
hot.mi7.co <- hot.co[which(hot.co$chrom=='scaffold-mi7'),]

# For plotting to evaluate what's going on with LD
ld.svmp <- read.table('LD_SVMP_region.txt',header=F)

### Plot SVMP region--------------------------------------------------------

options(scipen=5)

par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(win.10kb.mi1.cv$start,win.10kb.mi1.cv$mean,pch=20,col=alpha('lightblue3',0.15),ylim=c(0,0.2),ylab='Rho',xlab='Chromosome Position')
lines(win.100kb.mi1.cv$start,win.100kb.mi1.cv$mean,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(win.10kb.mi1.co$start,win.10kb.mi1.co$mean,pch=20,col=alpha('lightblue3',0.15),ylim=c(0,0.3),ylab='Rho',xlab='Chromosome Position')
lines(win.100kb.mi1.co$start,win.100kb.mi1.co$mean,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
# Zoomed, map points, 10 kb windows
plot(win.1kb.mi1.cv$start,win.1kb.mi1.cv$mean,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.35),ylab='Rho',xlab='Chromosome Position')
lines(win.10kb.mi1.cv$start,win.10kb.mi1.cv$mean,type='l')
#lines(ld.svmp$V2,ld.svmp$V4,col='navy')
for (h in hot.mi1.cv$start) {
  abline(v=h,col='darkorange')
}
vert <- -0.02
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
plot(win.1kb.mi1.co$start,win.1kb.mi1.co$mean,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.5),ylab='Rho',xlab='Chromosome Position')
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
plot(win.10kb.mi2.cv$start,win.10kb.mi2.cv$mean,pch=20,col=alpha('aquamarine3',0.15),ylim=c(0,0.2),ylab='Rho',xlab='Chromosome Position')
lines(win.100kb.mi2.cv$start,win.100kb.mi2.cv$mean,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(win.10kb.mi2.co$start,win.10kb.mi2.co$mean,pch=20,col=alpha('aquamarine3',0.15),ylim=c(0,0.25),ylab='Rho',xlab='Chromosome Position')
lines(win.100kb.mi2.co$start,win.100kb.mi2.co$mean,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
# Zoomed, map points, 10 kb windows
plot(win.1kb.mi2.cv$start,win.1kb.mi2.cv$mean,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.025,0.2),ylab='Rho',xlab='Chromosome Position')
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
plot(win.10kb.mi7.cv$start,win.10kb.mi7.cv$mean,pch=20,col=alpha('maroon',0.15),ylim=c(0,0.1),ylab='Rho',xlab='Chromosome Position')
lines(win.100kb.mi7.cv$start,win.100kb.mi7.cv$mean,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
plot(win.10kb.mi7.co$start,win.10kb.mi7.co$mean,pch=20,col=alpha('maroon',0.15),ylim=c(0,0.1),ylab='Rho',xlab='Chromosome Position')
lines(win.100kb.mi7.co$start,win.100kb.mi7.co$mean,type='l')
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
