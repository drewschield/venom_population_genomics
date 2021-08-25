############################################################################
# Genomic scans of fixed differences
############################################################################

### Goal: look for evidence of selection based on fixed differences
### between C. viridis and C.oreganus.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('./df/')

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

cv1co1.df.100kb <- read.table('window.100kb.df_prop.all.cv1.co1.txt',header=T)
cv1cv2.df.100kb <- read.table('window.100kb.df_prop.all.cv1.cv2.txt',header=T)
co1co2.df.100kb <- read.table('window.100kb.df_prop.all.co1.co2.txt',header=T)
cv1co1.df.10kb <- read.table('window.10kb.df_prop.all.cv1.co1.cnvMask.txt',header=T)
cv1cv2.df.10kb <- read.table('window.10kb.df_prop.all.cv1.cv2.txt',header=T)
co1co2.df.10kb <- read.table('window.10kb.df_prop.all.co1.co2.cnvMask.txt',header=T)
cv1co1.df.1kb <- read.table('window.1kb.df_prop.all.cv1.co1.cnvMask.txt',header=T)
cv1cv2.df.1kb <- read.table('window.1kb.df_prop.all.cv1.cv2.txt',header=T)
co1co2.df.1kb <- read.table('window.1kb.df_prop.all.co1.co2.cnvMask.txt',header=T)

cv1co1.df.mi1.10kb <- cv1co1.df.10kb[which(cv1co1.df.10kb$chrom=='scaffold-mi1'),]
cv1co1.df.mi2.10kb <- cv1co1.df.10kb[which(cv1co1.df.10kb$chrom=='scaffold-mi2'),]
cv1co1.df.mi7.10kb <- cv1co1.df.10kb[which(cv1co1.df.10kb$chrom=='scaffold-mi7'),]

cv1co1.df.mi1.1kb <- cv1co1.df.1kb[which(cv1co1.df.1kb$chrom=='scaffold-mi1'),]
cv1co1.df.mi2.1kb <- cv1co1.df.1kb[which(cv1co1.df.1kb$chrom=='scaffold-mi2'),]
cv1co1.df.mi7.1kb <- cv1co1.df.1kb[which(cv1co1.df.1kb$chrom=='scaffold-mi7'),]

cv1co1.df.mi7.250bp <- read.table('window.250bp.df_prop.scaffold-mi7.cv1.co1.cnvMask.txt',header=T)

## Calculate mean df from whole genome--------------------------------------

cv1co1.df.mean <- mean(cv1co1.df.10kb$df)
cv1cv2.df.mean <- mean(cv1cv2.df.10kb$df)
co1co2.df.mean <- mean(co1co2.df.10kb$df)

## SVMP; chromosome 9
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv1co1.df.mi1.10kb$start,cv1co1.df.mi1.10kb$df,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.2))
lines(cv1co1.df.mi1.100kb$start,cv1co1.df.mi1.100kb$df,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
abline(h=cv1co1.df.mean,lty=2,col='red')
plot(cv1cv2.df.mi1.10kb$start,cv1cv2.df.mi1.10kb$df,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.01))
lines(cv1cv2.df.mi1.100kb$start,cv1cv2.df.mi1.100kb$df,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
abline(h=cv1cv2.df.mean,lty=2,col='red')
plot(co1co2.df.mi1.10kb$start,co1co2.df.mi1.10kb$df,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.01))
lines(co1co2.df.mi1.100kb$start,co1co2.df.mi1.100kb$df,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
abline(h=co1co2.df.mean,lty=2,col='red')

plot(cv1co1.df.mi2.10kb$start,cv1co1.df.mi2.10kb$df,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.2))
lines(cv1co1.df.mi2.100kb$start,cv1co1.df.mi2.100kb$df,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
abline(h=cv1co1.df.mean,lty=2,col='red')
plot(cv1cv2.df.mi2.10kb$start,cv1cv2.df.mi2.10kb$df,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.01))
lines(cv1cv2.df.mi2.100kb$start,cv1cv2.df.mi2.100kb$df,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
abline(h=cv1cv2.df.mean,lty=2,col='red')
plot(co1co2.df.mi2.10kb$start,co1co2.df.mi2.10kb$df,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.2))
lines(co1co2.df.mi2.100kb$start,co1co2.df.mi2.100kb$df,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
abline(h=co1co2.df.mean,lty=2,col='red')


plot(cv1co1.df.mi7.10kb$start,cv1co1.df.mi7.10kb$df,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.2))
lines(cv1co1.df.mi7.100kb$start,cv1co1.df.mi7.100kb$df,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
abline(h=cv1co1.df.mean,lty=2,col='red')
plot(cv1cv2.df.mi7.10kb$start,cv1cv2.df.mi7.10kb$df,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.025))
lines(cv1cv2.df.mi7.100kb$start,cv1cv2.df.mi7.100kb$df,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
abline(h=cv1cv2.df.mean,lty=2,col='red')
plot(co1co2.df.mi7.10kb$start,co1co2.df.mi7.10kb$df,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='df',ylim=c(0,0.2))
lines(co1co2.df.mi7.100kb$start,co1co2.df.mi7.100kb$df,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
abline(h=co1co2.df.mean,lty=2,col='red')

# Zoomed, 1kb points, 10 kb windows
plot(cv1co1.df.mi1.1kb$start,cv1co1.df.mi1.1kb$df,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.2),ylab='df',xlab='Chromosome Position')
lines(cv1co1.df.mi1.10kb$start,cv1co1.df.mi1.10kb$df,type='l')
vert <- -0.015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.df.mean,lty=2,col='red')
plot(cv1cv2.df.mi1.1kb$start,cv1cv2.df.mi1.1kb$df,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.2),ylab='df',xlab='Chromosome Position')
lines(cv1cv2.df.mi1.10kb$start,cv1cv2.df.mi1.10kb$df,type='l')
vert <- -0.015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv1cv2.df.mean,lty=2,col='red')
plot(co1co2.df.mi1.1kb$start,co1co2.df.mi1.1kb$df,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.2),ylab='df',xlab='Chromosome Position')
lines(co1co2.df.mi1.10kb$start,co1co2.df.mi1.10kb$df,type='l')
vert <- -0.015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=co1co2.df.mean,lty=2,col='red')


plot(cv1co1.df.mi2.1kb$start,cv1co1.df.mi2.1kb$df,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.025,0.2),ylab='df',xlab='Chromosome Position')
lines(cv1co1.df.mi2.10kb$start,cv1co1.df.mi2.10kb$df,type='l')
vert <- -0.015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.df.mean,lty=2,col='red')

plot(cv1co1.df.mi7.250bp$start,cv1co1.df.mi7.250bp$df,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.025,0.35),ylab='df',xlab='Chromosome Position')
lines(cv1co1.df.mi7.1kb$start,cv1co1.df.mi7.1kb$df,type='l')
vert <- -0.015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.df.mean,lty=2,col='red')
