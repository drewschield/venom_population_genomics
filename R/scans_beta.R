############################################################################
# Genomic scans of ß statistics
############################################################################

### Goal: look for evidence of balancing selection using ß statistics, which
### measure clustering of intermediate-frequency polymorphisms surrounding a
### balanced allele. Estimates were obtaine using BetaScan (Siewert and Voight
### 2017).

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('./beta/results/')

library(scales)

options('stringsAsFactors'=FALSE)

### Read in general data locations------------------------------------------

# Gene locations
svmp <- read.table('./resources/gene_SVMP.bed',header=F)
svsp <- read.table('./resources/gene_SVSP.bed',header=F)
pla2 <- read.table('./resources/gene_PLA2.bed',header=F)

# Region locations
svmp.reg <- read.table('./resources/region_SVMP_scaffold-mi1.bed',header=F)
svsp.reg <- read.table('./resources/region_SVSP_scaffold-mi2.bed',header=F)
pla2.reg <- read.table('./resources/region_PLA2_scaffold-mi7.bed',header=F)

### Read in data------------------------------------------------------------

cv.beta.100kb <- read.table('cv1.phased.all.betascores.100kb.txt',header=T)
cv.beta.10kb <- read.table('cv1.phased.all.betascores.10kb.txt',header=T)
cv.beta.1kb <- read.table('cv1.phased.all.betascores.1kb.txt',header=T)

co.beta.100kb <- read.table('co1.phased.all.betascores.100kb.txt',header=T)
co.beta.10kb <- read.table('co1.phased.all.betascores.10kb.cnvMask.txt',header=T)
co.beta.1kb <- read.table('co1.phased.all.betascores.1kb.cnvMask.txt',header=T)

plot(cv.beta.10kb$Beta1.,type='l')

mean(as.numeric(cv.beta.10kb$Beta1),na.rm=T)
quantile(as.numeric(cv.beta.10kb$Beta1.),c(0.975,0.99),na.rm=T)
hist(as.numeric(cv.beta.10kb$Beta1.),breaks=100)

# Compare beta with iHS and rho
plot(cv.beta.10kb$Beta1.,cv1.rho.10kb$mean_corr,pch=20,ylim=c(0,60),col=alpha('black',0.05))
cor.test(as.numeric(cv.beta.10kb$Beta1.),cv1.rho.10kb$mean_corr,method='spearman')
plot(cv.beta.10kb$Beta1.,cv.ihs.10kb$iHS,pch=20,col=alpha('black',0.05))
cor.test(as.numeric(cv.beta.10kb$Beta1.),as.numeric(cv.ihs.10kb$iHS))
plot(cv.ihs.10kb$iHS,cv1.rho.10kb$mean_corr,pch=20,col=alpha('black',0.05),ylim=c(0,60))
cor.test(as.numeric(cv.ihs.10kb$iHS),as.numeric(cv1.rho.10kb$mean_corr),method='spearman')

### Parse venom chromosomes-------------------------------------------------

cv.beta.100kb.mi1 <- cv.beta.100kb[which(cv.beta.100kb$chrom=='scaffold-mi1'),]
cv.beta.100kb.mi2 <- cv.beta.100kb[which(cv.beta.100kb$chrom=='scaffold-mi2'),]
cv.beta.100kb.mi7 <- cv.beta.100kb[which(cv.beta.100kb$chrom=='scaffold-mi7'),]

cv.beta.10kb.mi1 <- cv.beta.10kb[which(cv.beta.10kb$chrom=='scaffold-mi1'),]
cv.beta.10kb.mi2 <- cv.beta.10kb[which(cv.beta.10kb$chrom=='scaffold-mi2'),]
cv.beta.10kb.mi7 <- cv.beta.10kb[which(cv.beta.10kb$chrom=='scaffold-mi7'),]

cv.beta.1kb.mi1 <- cv.beta.1kb[which(cv.beta.1kb$chrom=='scaffold-mi1'),]
cv.beta.1kb.mi2 <- cv.beta.1kb[which(cv.beta.1kb$chrom=='scaffold-mi2'),]
cv.beta.1kb.mi7 <- cv.beta.1kb[which(cv.beta.1kb$chrom=='scaffold-mi7'),]

co.beta.100kb.mi1 <- co.beta.100kb[which(co.beta.100kb$chrom=='scaffold-mi1'),]
co.beta.100kb.mi2 <- co.beta.100kb[which(co.beta.100kb$chrom=='scaffold-mi2'),]
co.beta.100kb.mi7 <- co.beta.100kb[which(co.beta.100kb$chrom=='scaffold-mi7'),]

co.beta.10kb.mi1 <- co.beta.10kb[which(co.beta.10kb$chrom=='scaffold-mi1'),]
co.beta.10kb.mi2 <- co.beta.10kb[which(co.beta.10kb$chrom=='scaffold-mi2'),]
co.beta.10kb.mi7 <- co.beta.10kb[which(co.beta.10kb$chrom=='scaffold-mi7'),]

co.beta.1kb.mi1 <- co.beta.1kb[which(co.beta.1kb$chrom=='scaffold-mi1'),]
co.beta.1kb.mi2 <- co.beta.1kb[which(co.beta.1kb$chrom=='scaffold-mi2'),]
co.beta.1kb.mi7 <- co.beta.1kb[which(co.beta.1kb$chrom=='scaffold-mi7'),]

cv.beta.250bp.mi7 <- read.table('cv1.phased.scaffold-mi7.betascores.250bp.txt',header=T)
co.beta.250bp.mi7 <- read.table('co1.phased.scaffold-mi7.betascores.250bp.cnvMask.txt',header=T)

### Plot scans--------------------------------------------------------------

# Set different panel arrangements for plotting together:
par(mfrow=c(2,2))
par(mfrow=c(3,1))

## SVMP; chromosome 9
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv.beta.10kb.mi1$start,cv.beta.10kb.mi1$Beta1.,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='ß',ylim=c(0,6))
lines(cv.beta.100kb.mi1$start,cv.beta.100kb.mi1$Beta1.,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(co.beta.10kb.mi1$start,co.beta.10kb.mi1$Beta1.,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='ß',ylim=c(0,6))
lines(co.beta.100kb.mi1$start,co.beta.100kb.mi1$Beta1.,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv.beta.1kb.mi1$start,cv.beta.1kb.mi1$Beta1.,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylab='ß',xlab='Chromosome Position',ylim=c(-2,15))
lines(cv.beta.10kb.mi1$start,cv.beta.10kb.mi1$Beta1.,type='l')
vert <- -1.2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
plot(co.beta.1kb.mi1$start,co.beta.1kb.mi1$Beta1.,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylab='ß',xlab='Chromosome Position',ylim=c(-2,15))
lines(co.beta.10kb.mi1$start,co.beta.10kb.mi1$Beta1.,type='l')
vert <- -1.2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

## SVSP; chromosome 10
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv.beta.10kb.mi2$start,cv.beta.10kb.mi2$Beta1.,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='ß',ylim=c(0,8))
lines(cv.beta.100kb.mi2$start,cv.beta.100kb.mi2$Beta1.,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(co.beta.10kb.mi2$start,co.beta.10kb.mi2$Beta1.,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='ß',ylim=c(0,8))
lines(co.beta.100kb.mi2$start,co.beta.100kb.mi2$Beta1.,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv.beta.1kb.mi2$start,cv.beta.1kb.mi2$Beta1.,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylab='ß',xlab='Chromosome Position',ylim=c(-2,15))
lines(cv.beta.10kb.mi2$start,cv.beta.10kb.mi2$Beta1.,type='l')
vert <- -1.2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
plot(co.beta.1kb.mi2$start,co.beta.1kb.mi2$Beta1.,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylab='ß',xlab='Chromosome Position',ylim=c(-2,15))
lines(co.beta.10kb.mi2$start,co.beta.10kb.mi2$Beta1.,type='l')
vert <- -1.2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}

## PLA2; chromosome 15
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv.beta.10kb.mi7$start,cv.beta.10kb.mi7$Beta1.,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='ß',ylim=c(-1,4))
lines(cv.beta.100kb.mi7$start,cv.beta.100kb.mi7$Beta1.,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
plot(co.beta.10kb.mi7$start,co.beta.10kb.mi7$Beta1.,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='ß',ylim=c(-1,4))
lines(co.beta.100kb.mi7$start,co.beta.100kb.mi7$Beta1.,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv.beta.250bp.mi7$start,cv.beta.250bp.mi7$Beta.,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylab='ß',xlab='Chromosome Position',ylim=c(-2,10))
lines(cv.beta.1kb.mi7$start,cv.beta.1kb.mi7$Beta1.,type='l')
vert <- -1.2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
plot(co.beta.250bp.mi7$start,co.beta.250bp.mi7$Beta1.,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylab='ß',xlab='Chromosome Position',ylim=c(-2,10))
lines(co.beta.1kb.mi7$start,co.beta.1kb.mi7$Beta1.,type='l')
vert <- -1.2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}

plot(cv.beta.100kb.mi1$Beta1.,cv.ihs.mi1.100kb$iHS,pch=20)
cor.test(cv.ihs.mi1.100kb$iHS,cv.beta.100kb.mi1$Beta1.)

plot(cv.beta.100kb.mi1$Beta1.,cv1.rho.100kb.mi1$mean_corr,pch=20)
cor.test(cv.beta.100kb.mi1$Beta1.,cv1.rho.100kb.mi1$mean_corr)

plot(cv.ihs.mi1.100kb$iHS,cv1.rho.100kb.mi1$mean_corr,pch=20)
cor.test(cv.ihs.mi1.100kb$iHS,cv1.rho.100kb.mi1$mean_corr)

plot(cv.beta.100kb$Beta1.,cv1.rho.100kb$mean_corr,pch=20)
cor.test(cv.beta.100kb$Beta1.,cv1.rho.100kb$mean_corr)
