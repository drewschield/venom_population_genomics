############################################################################
# Genomic scans haplotype diversity statistics
############################################################################

### Goal: look for evidence of selection based on extended haplotype lengths
### using the 'rehh' package and phased haplotypes for C. viridis and C.
### oreganus from the recombination study.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('./rehh/')

library(scales)

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

cv.ihs.100kb <- read.table('cv.all_ihs.100kb.txt',header=T)
cv.ihs.10kb <- read.table('cv.all_ihs.10kb.txt',header=T)
cv.ihs.1kb <- read.table('cv.all_ihs.1kb.txt',header=T)

co.ihs.100kb <- read.table('co.all_ihs.100kb.txt',header=T)
co.ihs.10kb <- read.table('co.all_ihs.10kb.cnvMask.txt',header=T)
co.ihs.1kb <- read.table('co.all_ihs.1kb.cnvMask.txt',header=T)

par(mfrow=c(2,1))
plot(cv.ihs.100kb$iHS,type='l')
plot(co.ihs.100kb$iHS,type='l')

cv.ihs.mi1.100kb <- cv.ihs.100kb[which(cv.ihs.100kb$chrom=='scaffold-mi1'),]
cv.ihs.mi2.100kb <- cv.ihs.100kb[which(cv.ihs.100kb$chrom=='scaffold-mi2'),]
cv.ihs.mi7.100kb <- cv.ihs.100kb[which(cv.ihs.100kb$chrom=='scaffold-mi7'),]

cv.ihs.mi1.10kb <- cv.ihs.10kb[which(cv.ihs.10kb$chrom=='scaffold-mi1'),]
cv.ihs.mi2.10kb <- cv.ihs.10kb[which(cv.ihs.10kb$chrom=='scaffold-mi2'),]
cv.ihs.mi7.10kb <- cv.ihs.10kb[which(cv.ihs.10kb$chrom=='scaffold-mi7'),]

cv.ihs.mi1.1kb <- cv.ihs.1kb[which(cv.ihs.1kb$chrom=='scaffold-mi1'),]
cv.ihs.mi2.1kb <- cv.ihs.1kb[which(cv.ihs.1kb$chrom=='scaffold-mi2'),]
cv.ihs.mi7.1kb <- cv.ihs.1kb[which(cv.ihs.1kb$chrom=='scaffold-mi7'),]

cv.ihs.mi7.250bp <- read.table('cv.mi7_ihs.250bp.txt',header=T)

co.ihs.mi1.100kb <- co.ihs.100kb[which(co.ihs.100kb$chrom=='scaffold-mi1'),]
co.ihs.mi2.100kb <- co.ihs.100kb[which(co.ihs.100kb$chrom=='scaffold-mi2'),]
co.ihs.mi7.100kb <- co.ihs.100kb[which(co.ihs.100kb$chrom=='scaffold-mi7'),]

co.ihs.mi1.10kb <- co.ihs.10kb[which(co.ihs.10kb$chrom=='scaffold-mi1'),]
co.ihs.mi2.10kb <- co.ihs.10kb[which(co.ihs.10kb$chrom=='scaffold-mi2'),]
co.ihs.mi7.10kb <- co.ihs.10kb[which(co.ihs.10kb$chrom=='scaffold-mi7'),]

co.ihs.mi1.1kb <- co.ihs.1kb[which(co.ihs.1kb$chrom=='scaffold-mi1'),]
co.ihs.mi2.1kb <- co.ihs.1kb[which(co.ihs.1kb$chrom=='scaffold-mi2'),]
co.ihs.mi7.1kb <- co.ihs.1kb[which(co.ihs.1kb$chrom=='scaffold-mi7'),]

co.ihs.mi7.250bp <- read.table('co.mi7_ihs.250bp.cnvMask.txt',header=T)


cv.ihs.mi1.100kb <- read.table('cv.mi1_ihs.100kb.txt',header=T)
cv.ihs.mi1.10kb <- read.table('cv.mi1_ihs.10kb.txt',header=T)
cv.ihs.mi1.1kb <- read.table('cv.mi1_ihs.1kb.txt',header=T)
cv.ihs.mi2.100kb <- read.table('cv.mi2_ihs.100kb.txt',header=T)
cv.ihs.mi2.10kb <- read.table('cv.mi2_ihs.10kb.txt',header=T)
cv.ihs.mi2.1kb <- read.table('cv.mi2_ihs.1kb.txt',header=T)
cv.ihs.mi7.100kb <- read.table('cv.mi7_ihs.100kb.txt',header=T)
cv.ihs.mi7.10kb <- read.table('cv.mi7_ihs.10kb.txt',header=T)
cv.ihs.mi7.1kb <- read.table('cv.mi7_ihs.1kb.txt',header=T)
cv.ihs.mi7.250bp <- read.table('cv.mi7_ihs.250bp.txt',header=T)

co.ihs.mi1.100kb <- read.table('co.mi1_ihs.100kb.txt',header=T)
co.ihs.mi1.10kb <- read.table('co.mi1_ihs.10kb.txt',header=T)
co.ihs.mi1.1kb <- read.table('co.mi1_ihs.1kb.txt',header=T)
co.ihs.mi2.100kb <- read.table('co.mi2_ihs.100kb.txt',header=T)
co.ihs.mi2.10kb <- read.table('co.mi2_ihs.10kb.txt',header=T)
co.ihs.mi2.1kb <- read.table('co.mi2_ihs.1kb.txt',header=T)
co.ihs.mi7.100kb <- read.table('co.mi7_ihs.100kb.txt',header=T)
co.ihs.mi7.10kb <- read.table('co.mi7_ihs.10kb.txt',header=T)
co.ihs.mi7.1kb <- read.table('co.mi7_ihs.1kb.txt',header=T)
co.ihs.mi7.250bp <- read.table('co.mi7_ihs.250bp.txt',header=T)

cvco.xpehh.mi1.100kb <- read.table('cvco.mi1_xpehh.100kb.txt',header=T)
cvco.xpehh.mi1.10kb <- read.table('cvco.mi1_xpehh.10kb.txt',header=T)
cvco.xpehh.mi1.1kb <- read.table('cvco.mi1_xpehh.1kb.txt',header=T)
cvco.xpehh.mi2.100kb <- read.table('cvco.mi2_xpehh.100kb.txt',header=T)
cvco.xpehh.mi2.10kb <- read.table('cvco.mi2_xpehh.10kb.txt',header=T)
cvco.xpehh.mi2.1kb <- read.table('cvco.mi2_xpehh.1kb.txt',header=T)
cvco.xpehh.mi7.100kb <- read.table('cvco.mi7_xpehh.100kb.txt',header=T)
cvco.xpehh.mi7.10kb <- read.table('cvco.mi7_xpehh.10kb.txt',header=T)
cvco.xpehh.mi7.1kb <- read.table('cvco.mi7_xpehh.1kb.txt',header=T)
cvco.xpehh.mi7.250bp <- read.table('cvco.mi7_xpehh.250bp.txt',header=T)

### Plot scans--------------------------------------------------------------

# Set different panel arrangements for plotting together:
par(mfrow=c(2,2))
par(mfrow=c(3,1))

###-------------------------------------------------------------------------
### iHS statistic
###-------------------------------------------------------------------------

## SVMP; chromosome 9
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv.ihs.mi1.10kb$start,cv.ihs.mi1.10kb$iHS,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='iHS')
lines(cv.ihs.mi1.100kb$start,cv.ihs.mi1.100kb$iHS,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(co.ihs.mi1.10kb$start,co.ihs.mi1.10kb$iHS,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='iHS')
lines(co.ihs.mi1.100kb$start,co.ihs.mi1.100kb$iHS,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv.ihs.mi1.1kb$start,cv.ihs.mi1.1kb$iHS,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-1.5,3),ylab='iHS',xlab='Chromosome Position')
lines(cv.ihs.mi1.10kb$start,cv.ihs.mi1.10kb$iHS,type='l')
vert <- -1.2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
plot(co.ihs.mi1.1kb$start,co.ihs.mi1.1kb$iHS,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-1.5,3),ylab='iHS',xlab='Chromosome Position')
lines(co.ihs.mi1.10kb$start,co.ihs.mi1.10kb$iHS,type='l')
vert <- -1.2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

## SVSP; chromosome 10
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv.ihs.mi2.10kb$start,cv.ihs.mi2.10kb$iHS,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='iHS')
lines(cv.ihs.mi2.100kb$start,cv.ihs.mi2.100kb$iHS,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(co.ihs.mi2.10kb$start,co.ihs.mi2.10kb$iHS,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='iHS')
lines(co.ihs.mi2.100kb$start,co.ihs.mi2.100kb$iHS,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv.ihs.mi2.1kb$start,cv.ihs.mi2.1kb$iHS,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-1.5,3),ylab='iHS',xlab='Chromosome Position')
lines(cv.ihs.mi2.10kb$start,cv.ihs.mi2.10kb$iHS,type='l')
vert <- -1.2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
plot(co.ihs.mi2.1kb$start,co.ihs.mi2.1kb$iHS,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-1.5,3),ylab='iHS',xlab='Chromosome Position')
lines(co.ihs.mi2.10kb$start,co.ihs.mi2.10kb$iHS,type='l')
vert <- -1.2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}

## PLA2; chromosome 15
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv.ihs.mi7.10kb$start,cv.ihs.mi7.10kb$iHS,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='iHS')
lines(cv.ihs.mi7.100kb$start,cv.ihs.mi7.100kb$iHS,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
plot(co.ihs.mi7.10kb$start,co.ihs.mi7.10kb$iHS,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='iHS')
lines(co.ihs.mi7.100kb$start,co.ihs.mi7.100kb$iHS,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
# Zoomed, 250 bp points, 1 kb windows
plot(cv.ihs.mi7.250bp$start,cv.ihs.mi7.250bp$iHS,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-1.5,3),ylab='iHS',xlab='Chromosome Position')
lines(cv.ihs.mi7.1kb$start,cv.ihs.mi7.1kb$iHS,type='l')
vert <- -1.2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
plot(co.ihs.mi7.250bp$start,co.ihs.mi7.250bp$iHS,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-1.5,3),ylab='iHS',xlab='Chromosome Position')
lines(co.ihs.mi7.1kb$start,co.ihs.mi7.1kb$iHS,type='l')
vert <- -1.2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
