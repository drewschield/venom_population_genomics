############################################################################
# Genomic scans of neutrality tests (Tajima's D)
############################################################################

### Goal: examine patterns of Tajima's D across venom gene loci to look for
### evidence of departures from neutrality, as expected if positive and/or
### balancing selection have been acting on venom.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('./tajima_d/')

library(scales)
install.packages("cowplot")
library(cowplot)

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

# Genome-wide 10 kb windows
cv1.taj.10kb <- read.table('cv.colorado.all.10kb.Tajima.D',header=T)
cv2.taj.10kb <- read.table('cv.montana.all.10kb.Tajima.D',header=T)
co1.taj.10kb <- read.table('co.california.all.10kb.cnvMask.Tajima.D',header=T)
co2.taj.10kb <- read.table('co.idaho.all.10kb.cnvMask.Tajima.D',header=T)

# Genome-wide 1 kb windows
cv1.taj.1kb <- read.table('cv.colorado.all.1kb.Tajima.D',header=T)
cv2.taj.1kb <- read.table('cv.montana.all.1kb.Tajima.D',header=T)
co1.taj.1kb <- read.table('co.california.all.1kb.cnvMask.Tajima.D',header=T)
co2.taj.1kb <- read.table('co.idaho.all.1kb.cnvMask.Tajima.D',header=T)

# 100 kb windows
cv1.taj.mi1.100kb <- read.table('cv.colorado.scaffold-mi1.100kb.Tajima.D',header=T)
cv2.taj.mi1.100kb <- read.table('cv.montana.scaffold-mi1.100kb.Tajima.D',header=T)
co1.taj.mi1.100kb <- read.table('co.california.scaffold-mi1.100kb.Tajima.D',header=T)
co2.taj.mi1.100kb <- read.table('co.idaho.scaffold-mi1.100kb.Tajima.D',header=T)

cv1.taj.mi2.100kb <- read.table('cv.colorado.scaffold-mi2.100kb.Tajima.D',header=T)
cv2.taj.mi2.100kb <- read.table('cv.montana.scaffold-mi2.100kb.Tajima.D',header=T)
co1.taj.mi2.100kb <- read.table('co.california.scaffold-mi2.100kb.Tajima.D',header=T)
co2.taj.mi2.100kb <- read.table('co.idaho.scaffold-mi2.100kb.Tajima.D',header=T)

cv1.taj.mi7.100kb <- read.table('cv.colorado.scaffold-mi7.100kb.Tajima.D',header=T)
cv2.taj.mi7.100kb <- read.table('cv.montana.scaffold-mi7.100kb.Tajima.D',header=T)
co1.taj.mi7.100kb <- read.table('co.california.scaffold-mi7.100kb.Tajima.D',header=T)
co2.taj.mi7.100kb <- read.table('co.idaho.scaffold-mi7.100kb.Tajima.D',header=T)

# 10 kb windows
cv1.taj.mi1.10kb <- read.table('cv.colorado.scaffold-mi1.10kb.Tajima.D',header=T)
cv2.taj.mi1.10kb <- read.table('cv.montana.scaffold-mi1.10kb.Tajima.D',header=T)
co1.taj.mi1.10kb <- co1.taj.10kb[which(co1.taj.10kb$CHROM=='scaffold-mi1'),]
co2.taj.mi1.10kb <- co2.taj.10kb[which(co2.taj.10kb$CHROM=='scaffold-mi1'),]

cv1.taj.mi2.10kb <- read.table('cv.colorado.scaffold-mi2.10kb.Tajima.D',header=T)
cv2.taj.mi2.10kb <- read.table('cv.montana.scaffold-mi2.10kb.Tajima.D',header=T)
co1.taj.mi2.10kb <- co1.taj.10kb[which(co1.taj.10kb$CHROM=='scaffold-mi2'),]
co2.taj.mi2.10kb <- co2.taj.10kb[which(co2.taj.10kb$CHROM=='scaffold-mi2'),]

cv1.taj.mi7.10kb <- read.table('cv.colorado.scaffold-mi7.10kb.Tajima.D',header=T)
cv2.taj.mi7.10kb <- read.table('cv.montana.scaffold-mi7.10kb.Tajima.D',header=T)
co1.taj.mi7.10kb <- co1.taj.10kb[which(co1.taj.10kb$CHROM=='scaffold-mi7'),]
co2.taj.mi7.10kb <- co2.taj.10kb[which(co2.taj.10kb$CHROM=='scaffold-mi7'),]

# 1 kb windows
cv1.taj.mi1.1kb <- read.table('cv.colorado.scaffold-mi1.1kb.Tajima.D',header=T)
cv2.taj.mi1.1kb <- read.table('cv.montana.scaffold-mi1.1kb.Tajima.D',header=T)
co1.taj.mi1.1kb <- co1.taj.1kb[which(co1.taj.1kb$CHROM=='scaffold-mi1'),]
co2.taj.mi1.1kb <- co2.taj.1kb[which(co2.taj.1kb$CHROM=='scaffold-mi1'),]

cv1.taj.mi2.1kb <- read.table('cv.colorado.scaffold-mi2.1kb.Tajima.D',header=T)
cv2.taj.mi2.1kb <- read.table('cv.montana.scaffold-mi2.1kb.Tajima.D',header=T)
co1.taj.mi2.1kb <- co1.taj.1kb[which(co1.taj.1kb$CHROM=='scaffold-mi2'),]
co2.taj.mi2.1kb <- co2.taj.1kb[which(co2.taj.1kb$CHROM=='scaffold-mi2'),]

cv1.taj.mi7.1kb <- read.table('cv.colorado.scaffold-mi7.1kb.Tajima.D',header=T)
cv2.taj.mi7.1kb <- read.table('cv.montana.scaffold-mi7.1kb.Tajima.D',header=T)
co1.taj.mi7.1kb <- co1.taj.1kb[which(co1.taj.1kb$CHROM=='scaffold-mi7'),]
co2.taj.mi7.1kb <- co2.taj.1kb[which(co2.taj.1kb$CHROM=='scaffold-mi7'),]

# 250 bp windows
cv1.taj.mi7.250bp <- read.table('cv.colorado.scaffold-mi7.250bp.Tajima.D',header=T)
cv2.taj.mi7.250bp <- read.table('cv.montana.scaffold-mi7.250bp.Tajima.D',header=T)
co1.taj.mi7.250bp <- read.table('co.california.scaffold-mi7.250bp.cnvMask.Tajima.D',header=T)
co2.taj.mi7.250bp <- read.table('co.idaho.scaffold-mi7.250bp.cnvMask.Tajima.D',header=T)

### Calculated genome-wide summary statistics-------------------------------

mean(cv1.taj.10kb$TajimaD,na.rm=T)
mean(cv2.taj.10kb$TajimaD,na.rm=T)
mean(co1.taj.10kb$TajimaD,na.rm=T)
mean(co2.taj.10kb$TajimaD,na.rm=T)

sd(cv1.taj.10kb$TajimaD,na.rm=T)
sd(cv2.taj.10kb$TajimaD,na.rm=T)
sd(co1.taj.10kb$TajimaD,na.rm=T)
sd(co2.taj.10kb$TajimaD,na.rm=T)

### Statistical comparison of venom D vs genome-wide------------------------

cv1.taj.svmp <- cv1.taj.mi1.10kb[which(cv1.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & cv1.taj.mi1.10kb$BIN_START<svmp.reg$V3),]
cv2.taj.svmp <- cv2.taj.mi1.10kb[which(cv2.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & cv2.taj.mi1.10kb$BIN_START<svmp.reg$V3),]
co1.taj.svmp <- co1.taj.mi1.10kb[which(co1.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & co1.taj.mi1.10kb$BIN_START<svmp.reg$V3),]
co2.taj.svmp <- co2.taj.mi1.10kb[which(co2.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & co2.taj.mi1.10kb$BIN_START<svmp.reg$V3),]

mean(cv1.taj.svmp$TajimaD,na.rm=T)
mean(cv2.taj.svmp$TajimaD,na.rm=T)
mean(co1.taj.svmp$TajimaD,na.rm=T)
mean(co2.taj.svmp$TajimaD,na.rm=T)

sd(cv1.taj.svmp$TajimaD,na.rm=T)
sd(cv2.taj.svmp$TajimaD,na.rm=T)
sd(co1.taj.svmp$TajimaD,na.rm=T)
sd(co2.taj.svmp$TajimaD,na.rm=T)

t.test(cv1.taj.10kb$TajimaD,cv1.taj.svmp$TajimaD)
t.test(cv2.taj.10kb$TajimaD,cv2.taj.svmp$TajimaD)
t.test(co1.taj.10kb$TajimaD,co1.taj.svmp$TajimaD)
t.test(co2.taj.10kb$TajimaD,co2.taj.svmp$TajimaD)

### Statistical comparison of venom D vs non-venom paralogs-----------------

svmp.nv <- read.table('./resources/non-venom_paralogs_SVMP.bed',header=F)
svsp.nv <- read.table('./resources/non-venom_paralogs_SVSP.bed',header=F)
pla2.nv <- read.table('./resources/non-venom_paralogs_PLA2.bed',header=F)

cv1.taj.nv <- cv1.taj.1kb[which(cv1.taj.1kb$CHROM==svmp.nv$V1 & cv1.taj.1kb$BIN_START>=svmp.nv$V2 & cv1.taj.1kb$BIN_START<=svmp.nv$V3),]

boxplot(cv1.taj.svmp$TajimaD,cv1.taj.nv$TajimaD,cv1.taj.1kb$TajimaD)

### Plot scans--------------------------------------------------------------

###-------------------------------------------------------------------------
### SVMP; microchromosome 1
###-------------------------------------------------------------------------

## Tajima's D
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv1.taj.mi1.10kb$BIN_START,cv1.taj.mi1.10kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(cv1.taj.mi1.100kb$BIN_START,cv1.taj.mi1.100kb$TajimaD,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(cv2.taj.mi1.10kb$BIN_START,cv2.taj.mi1.10kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(cv2.taj.mi1.100kb$BIN_START,cv2.taj.mi1.100kb$TajimaD,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(co1.taj.mi1.10kb$BIN_START,co1.taj.mi1.10kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(co1.taj.mi1.100kb$BIN_START,co1.taj.mi1.100kb$TajimaD,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(co2.taj.mi1.10kb$BIN_START,co2.taj.mi1.10kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(co2.taj.mi1.100kb$BIN_START,co2.taj.mi1.100kb$TajimaD,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv1.taj.mi1.1kb$BIN_START,cv1.taj.mi1.1kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(cv1.taj.mi1.10kb$BIN_START,cv1.taj.mi1.10kb$TajimaD,type='l')
vert <- -2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
plot(cv2.taj.mi1.1kb$BIN_START,cv2.taj.mi1.1kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(cv2.taj.mi1.10kb$BIN_START,cv2.taj.mi1.10kb$TajimaD,type='l')
vert <- -2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
plot(co1.taj.mi1.1kb$BIN_START,co1.taj.mi1.1kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(co1.taj.mi1.10kb$BIN_START,co1.taj.mi1.10kb$TajimaD,type='l')
vert <- -2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
plot(co2.taj.mi1.1kb$BIN_START,co2.taj.mi1.1kb$TajimaD,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(co2.taj.mi1.10kb$BIN_START,co2.taj.mi1.10kb$TajimaD,type='l')
vert <- -2
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

###-------------------------------------------------------------------------
### SVSP; microchromosome 2
###-------------------------------------------------------------------------

## Tajima's D
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv1.taj.mi2.10kb$BIN_START,cv1.taj.mi2.10kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(cv1.taj.mi2.100kb$BIN_START,cv1.taj.mi2.100kb$TajimaD,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(cv2.taj.mi2.10kb$BIN_START,cv2.taj.mi2.10kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(cv2.taj.mi2.100kb$BIN_START,cv2.taj.mi2.100kb$TajimaD,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(co1.taj.mi2.10kb$BIN_START,co1.taj.mi2.10kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(co1.taj.mi2.100kb$BIN_START,co1.taj.mi2.100kb$TajimaD,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(co2.taj.mi2.10kb$BIN_START,co2.taj.mi2.10kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(co2.taj.mi2.100kb$BIN_START,co2.taj.mi2.100kb$TajimaD,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv1.taj.mi2.1kb$BIN_START,cv1.taj.mi2.1kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(cv1.taj.mi2.10kb$BIN_START,cv1.taj.mi2.10kb$TajimaD,type='l')
vert <- -2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
plot(cv2.taj.mi2.1kb$BIN_START,cv2.taj.mi2.1kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(cv2.taj.mi2.10kb$BIN_START,cv2.taj.mi2.10kb$TajimaD,type='l')
vert <- -2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
plot(co1.taj.mi2.1kb$BIN_START,co1.taj.mi2.1kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(co1.taj.mi2.10kb$BIN_START,co1.taj.mi2.10kb$TajimaD,type='l')
vert <- -2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
plot(co2.taj.mi2.1kb$BIN_START,co2.taj.mi2.1kb$TajimaD,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(co2.taj.mi2.10kb$BIN_START,co2.taj.mi2.10kb$TajimaD,type='l')
vert <- -2
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}


###-------------------------------------------------------------------------
### PLA2; microchromosome 7
###-------------------------------------------------------------------------

## Tajima's D
par(mfrow=c(2,2))
# 10 kb points, 100 kb windows
plot(cv1.taj.mi7.10kb$BIN_START,cv1.taj.mi7.10kb$TajimaD,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(cv1.taj.mi7.100kb$BIN_START,cv1.taj.mi7.100kb$TajimaD,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
plot(cv2.taj.mi7.10kb$BIN_START,cv2.taj.mi7.10kb$TajimaD,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(cv2.taj.mi7.100kb$BIN_START,cv2.taj.mi7.100kb$TajimaD,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
plot(co1.taj.mi7.10kb$BIN_START,co1.taj.mi7.10kb$TajimaD,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(co1.taj.mi7.100kb$BIN_START,co1.taj.mi7.100kb$TajimaD,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
plot(co2.taj.mi7.10kb$BIN_START,co2.taj.mi7.10kb$TajimaD,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab="Tajima's D")
lines(co2.taj.mi7.100kb$BIN_START,co2.taj.mi7.100kb$TajimaD,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
# Zoomed, 1kb points, 10 kb windows
plot(cv1.taj.mi7.250bp$BIN_START,cv1.taj.mi7.250bp$TajimaD,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(cv1.taj.mi7.1kb$BIN_START,cv1.taj.mi7.1kb$TajimaD,type='l')
vert <- -2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
plot(cv2.taj.mi7.250bp$BIN_START,cv2.taj.mi7.250bp$TajimaD,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-3.5,3.5),ylab="Tajima's D",xlab='Chromosome Position')
lines(cv2.taj.mi7.1kb$BIN_START,cv2.taj.mi7.1kb$TajimaD,type='l')
vert <- -2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
plot(co1.taj.mi7.250bp$BIN_START,co1.taj.mi7.250bp$TajimaD,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-3,3),ylab="Tajima's D",xlab='Chromosome Position')
lines(co1.taj.mi7.1kb$BIN_START,co1.taj.mi7.1kb$TajimaD,type='l')
vert <- -2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
plot(co2.taj.mi7.250bp$BIN_START,co2.taj.mi7.250bp$TajimaD,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-3.5,3.5),ylab="Tajima's D",xlab='Chromosome Position')
lines(co2.taj.mi7.1kb$BIN_START,co2.taj.mi7.1kb$TajimaD,type='l')
vert <- -2
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}


