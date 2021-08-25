############################################################################
# Comparison of selection statistics in venom regions to backgrounds
############################################################################

### Goal: examine evidence for selection in venom regions by comparing the
### distributions of various statistics in venom clusters to genome-wide,
### chromosome-specific, and non-venom paralog backgrounds.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

install.packages("cowplot")
library(scales)
library(ggplot2)
library(cowplot)
library(gridExtra)

# Gene locations
svmp <- read.table('./resources/gene_SVMP.bed',header=F)
svsp <- read.table('./resources/gene_SVSP.bed',header=F)
pla2 <- read.table('./resources/gene_PLA2.bed',header=F)

# Region locations
svmp.reg <- read.table('./resources/region_SVMP_scaffold-mi1.bed',header=F)
svsp.reg <- read.table('./resources/region_SVSP_scaffold-mi2.bed',header=F)
pla2.reg <- read.table('./resources/region_PLA2_scaffold-mi7.bed',header=F)

### 1. Tajima's D-----------------------------------------------------------

setwd('./tajima_d/')

### 1.1 Read in and parse Tajima's D results--------------------------------

# Genome-wide 10 kb windows
cv1.taj.10kb <- read.table('cv.colorado.all.10kb.Tajima.D',header=T)
cv2.taj.10kb <- read.table('cv.montana.all.10kb.Tajima.D',header=T)
co1.taj.10kb <- read.table('co.california.all.10kb.cnvMask.Tajima.D',header=T)
co2.taj.10kb <- read.table('co.idaho.all.10kb.cnvMask.Tajima.D',header=T)

# Venom chromosome 10/1 kb windows
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

cv1.taj.mi7.1kb <- read.table('cv.colorado.scaffold-mi7.1kb.Tajima.D',header=T)
cv2.taj.mi7.1kb <- read.table('cv.montana.scaffold-mi7.1kb.Tajima.D',header=T)
co1.taj.mi7.1kb <- co1.taj.1kb[which(co1.taj.1kb$CHROM=='scaffold-mi7'),]
co2.taj.mi7.1kb <- co2.taj.1kb[which(co2.taj.1kb$CHROM=='scaffold-mi7'),]

# Parse venom regions
cv1.taj.svmp <- cv1.taj.mi1.10kb[which(cv1.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & cv1.taj.mi1.10kb$BIN_START<=svmp.reg$V3),]
cv2.taj.svmp <- cv2.taj.mi1.10kb[which(cv2.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & cv2.taj.mi1.10kb$BIN_START<=svmp.reg$V3),]
co1.taj.svmp <- co1.taj.mi1.10kb[which(co1.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & co1.taj.mi1.10kb$BIN_START<=svmp.reg$V3),]
co2.taj.svmp <- co2.taj.mi1.10kb[which(co2.taj.mi1.10kb$BIN_START>=svmp.reg$V2 & co2.taj.mi1.10kb$BIN_START<=svmp.reg$V3),]

cv1.taj.svsp <- cv1.taj.mi2.10kb[which(cv1.taj.mi2.10kb$BIN_START>=svsp.reg$V2 & cv1.taj.mi2.10kb$BIN_START<=svsp.reg$V3),]
cv2.taj.svsp <- cv2.taj.mi2.10kb[which(cv2.taj.mi2.10kb$BIN_START>=svsp.reg$V2 & cv2.taj.mi2.10kb$BIN_START<=svsp.reg$V3),]
co1.taj.svsp <- co1.taj.mi2.10kb[which(co1.taj.mi2.10kb$BIN_START>=svsp.reg$V2 & co1.taj.mi2.10kb$BIN_START<=svsp.reg$V3),]
co2.taj.svsp <- co2.taj.mi2.10kb[which(co2.taj.mi2.10kb$BIN_START>=svsp.reg$V2 & co2.taj.mi2.10kb$BIN_START<=svsp.reg$V3),]

cv1.taj.pla2 <- cv1.taj.mi7.1kb[which(cv1.taj.mi7.1kb$BIN_START>=pla2.reg$V2 & cv1.taj.mi7.1kb$BIN_START<=pla2.reg$V3),]
cv2.taj.pla2 <- cv2.taj.mi7.1kb[which(cv2.taj.mi7.1kb$BIN_START>=pla2.reg$V2 & cv2.taj.mi7.1kb$BIN_START<=pla2.reg$V3),]
co1.taj.pla2 <- co1.taj.mi7.1kb[which(co1.taj.mi7.1kb$BIN_START>=pla2.reg$V2 & co1.taj.mi7.1kb$BIN_START<=pla2.reg$V3),]
co2.taj.pla2 <- co2.taj.mi7.1kb[which(co2.taj.mi7.1kb$BIN_START>=pla2.reg$V2 & co2.taj.mi7.1kb$BIN_START<=pla2.reg$V3),]

# Parse chromosome backgrounds
cv1.taj.mi1.bk <- cv1.taj.mi1.10kb[which(cv1.taj.mi1.10kb$BIN_START<svmp.reg$V2 | cv1.taj.mi1.10kb$BIN_START>svmp.reg$V3),]
cv2.taj.mi1.bk <- cv2.taj.mi1.10kb[which(cv2.taj.mi1.10kb$BIN_START<svmp.reg$V2 | cv2.taj.mi1.10kb$BIN_START>svmp.reg$V3),]
co1.taj.mi1.bk <- co1.taj.mi1.10kb[which(co1.taj.mi1.10kb$BIN_START<svmp.reg$V2 | co1.taj.mi1.10kb$BIN_START>svmp.reg$V3),]
co2.taj.mi1.bk <- co2.taj.mi1.10kb[which(co2.taj.mi1.10kb$BIN_START<svmp.reg$V2 | co2.taj.mi1.10kb$BIN_START>svmp.reg$V3),]

cv1.taj.mi2.bk <- cv1.taj.mi2.10kb[which(cv1.taj.mi2.10kb$BIN_START<svmp.reg$V2 | cv1.taj.mi2.10kb$BIN_START>svmp.reg$V3),]
cv2.taj.mi2.bk <- cv2.taj.mi2.10kb[which(cv2.taj.mi2.10kb$BIN_START<svmp.reg$V2 | cv2.taj.mi2.10kb$BIN_START>svmp.reg$V3),]
co1.taj.mi2.bk <- co1.taj.mi2.10kb[which(co1.taj.mi2.10kb$BIN_START<svmp.reg$V2 | co1.taj.mi2.10kb$BIN_START>svmp.reg$V3),]
co2.taj.mi2.bk <- co2.taj.mi2.10kb[which(co2.taj.mi2.10kb$BIN_START<svmp.reg$V2 | co2.taj.mi2.10kb$BIN_START>svmp.reg$V3),]

cv1.taj.mi7.bk <- cv1.taj.mi7.10kb[which(cv1.taj.mi7.10kb$BIN_START<svmp.reg$V2 | cv1.taj.mi7.10kb$BIN_START>svmp.reg$V3),]
cv2.taj.mi7.bk <- cv2.taj.mi7.10kb[which(cv2.taj.mi7.10kb$BIN_START<svmp.reg$V2 | cv2.taj.mi7.10kb$BIN_START>svmp.reg$V3),]
co1.taj.mi7.bk <- co1.taj.mi7.10kb[which(co1.taj.mi7.10kb$BIN_START<svmp.reg$V2 | co1.taj.mi7.10kb$BIN_START>svmp.reg$V3),]
co2.taj.mi7.bk <- co2.taj.mi7.10kb[which(co2.taj.mi7.10kb$BIN_START<svmp.reg$V2 | co2.taj.mi7.10kb$BIN_START>svmp.reg$V3),]

# Read in non-venom paralog backgrounds
cv1.taj.nv.svmp <- read.table('cv.colorado.non-venom_SVMP.10kb.Tajima.D',header=T)
cv2.taj.nv.svmp <- read.table('cv.montana.non-venom_SVMP.10kb.Tajima.D',header=T)
co1.taj.nv.svmp <- read.table('co.california.non-venom_SVMP.10kb.Tajima.D',header=T)
co2.taj.nv.svmp <- read.table('co.idaho.non-venom_SVMP.10kb.Tajima.D',header=T)

cv1.taj.nv.svsp <- read.table('cv.colorado.non-venom_SVSP.10kb.Tajima.D',header=T)
cv2.taj.nv.svsp <- read.table('cv.montana.non-venom_SVSP.10kb.Tajima.D',header=T)
co1.taj.nv.svsp <- read.table('co.california.non-venom_SVSP.10kb.Tajima.D',header=T)
co2.taj.nv.svsp <- read.table('co.idaho.non-venom_SVSP.10kb.Tajima.D',header=T)

cv1.taj.nv.pla2 <- read.table('cv.colorado.non-venom_PLA2.1kb.Tajima.D',header=T)
cv2.taj.nv.pla2 <- read.table('cv.montana.non-venom_PLA2.1kb.Tajima.D',header=T)
co1.taj.nv.pla2 <- read.table('co.california.non-venom_PLA2.1kb.Tajima.D',header=T)
co2.taj.nv.pla2 <- read.table('co.idaho.non-venom_PLA2.1kb.Tajima.D',header=T)

### 1.2 Calculate genome-wide summary statistics----------------------------

mean(cv1.taj.10kb$TajimaD,na.rm=T)
mean(cv2.taj.10kb$TajimaD,na.rm=T)
mean(co1.taj.10kb$TajimaD,na.rm=T)
mean(co2.taj.10kb$TajimaD,na.rm=T)

sd(cv1.taj.10kb$TajimaD,na.rm=T)
sd(cv2.taj.10kb$TajimaD,na.rm=T)
sd(co1.taj.10kb$TajimaD,na.rm=T)
sd(co2.taj.10kb$TajimaD,na.rm=T)

### 1.3 Statistical comparison of venom D vs backgrounds---------------------

## SVMP
mean(cv1.taj.svmp$TajimaD,na.rm=T)
mean(cv2.taj.svmp$TajimaD,na.rm=T)
mean(co1.taj.svmp$TajimaD,na.rm=T)
mean(co2.taj.svmp$TajimaD,na.rm=T)

sd(cv1.taj.svmp$TajimaD,na.rm=T)
sd(cv2.taj.svmp$TajimaD,na.rm=T)
sd(co1.taj.svmp$TajimaD,na.rm=T)
sd(co2.taj.svmp$TajimaD,na.rm=T)

## Chr 9
mean(cv1.taj.mi1.bk$TajimaD,na.rm=T)
sd(cv1.taj.mi1.bk$TajimaD,na.rm=T)
mean(cv2.taj.mi1.bk$TajimaD,na.rm=T)
sd(cv2.taj.mi1.bk$TajimaD,na.rm=T)
mean(co1.taj.mi1.bk$TajimaD,na.rm=T)
sd(co1.taj.mi1.bk$TajimaD,na.rm=T)
mean(co2.taj.mi1.bk$TajimaD,na.rm=T)
sd(co2.taj.mi1.bk$TajimaD,na.rm=T)

## Non-venom
mean(cv1.taj.nv.svmp$TajimaD,na.rm=T)
sd(cv1.taj.nv.svmp$TajimaD,na.rm=T)
mean(cv2.taj.nv.svmp$TajimaD,na.rm=T)
sd(cv2.taj.nv.svmp$TajimaD,na.rm=T)
mean(co1.taj.nv.svmp$TajimaD,na.rm=T)
sd(co1.taj.nv.svmp$TajimaD,na.rm=T)
mean(co2.taj.nv.svmp$TajimaD,na.rm=T)
sd(co2.taj.nv.svmp$TajimaD,na.rm=T)

## T-tests
# Genome-wide
t.test(cv1.taj.10kb$TajimaD,cv1.taj.svmp$TajimaD)
t.test(cv2.taj.10kb$TajimaD,cv2.taj.svmp$TajimaD)
t.test(co1.taj.10kb$TajimaD,co1.taj.svmp$TajimaD)
t.test(co2.taj.10kb$TajimaD,co2.taj.svmp$TajimaD)
# Chromosome 9
t.test(cv1.taj.mi1.bk$TajimaD,cv1.taj.svmp$TajimaD)
t.test(cv2.taj.mi1.bk$TajimaD,cv2.taj.svmp$TajimaD)
t.test(co1.taj.mi1.bk$TajimaD,co1.taj.svmp$TajimaD)
t.test(co2.taj.mi1.bk$TajimaD,co2.taj.svmp$TajimaD)
# Non-venom paralogs
t.test(cv1.taj.nv.svmp$TajimaD,cv1.taj.svmp$TajimaD)
t.test(cv2.taj.nv.svmp$TajimaD,cv2.taj.svmp$TajimaD)
t.test(co1.taj.nv.svmp$TajimaD,co1.taj.svmp$TajimaD)
t.test(co2.taj.nv.svmp$TajimaD,co2.taj.svmp$TajimaD)

## SVSP
mean(cv1.taj.10kb$TajimaD,na.rm=T)
sd(cv1.taj.10kb$TajimaD,na.rm=T)
mean(cv2.taj.10kb$TajimaD,na.rm=T)
sd(cv2.taj.10kb$TajimaD,na.rm=T)
mean(co1.taj.10kb$TajimaD,na.rm=T)
sd(co1.taj.10kb$TajimaD,na.rm=T)
mean(co2.taj.10kb$TajimaD,na.rm=T)
sd(co2.taj.10kb$TajimaD,na.rm=T)

mean(cv1.taj.mi2.10kb$TajimaD,na.rm=T)
sd(cv1.taj.mi2.10kb$TajimaD,na.rm=T)
mean(cv2.taj.mi2.10kb$TajimaD,na.rm=T)
sd(cv2.taj.mi2.10kb$TajimaD,na.rm=T)
mean(co1.taj.mi2.10kb$TajimaD,na.rm=T)
sd(co1.taj.mi2.10kb$TajimaD,na.rm=T)
mean(co2.taj.mi2.10kb$TajimaD,na.rm=T)
sd(co2.taj.mi2.10kb$TajimaD,na.rm=T)

mean(cv1.taj.nv.svsp$TajimaD,na.rm=T)
sd(cv1.taj.nv.svsp$TajimaD,na.rm=T)
mean(cv2.taj.nv.svsp$TajimaD,na.rm=T)
sd(cv2.taj.nv.svsp$TajimaD,na.rm=T)
mean(co1.taj.nv.svsp$TajimaD,na.rm=T)
sd(co1.taj.nv.svsp$TajimaD,na.rm=T)
mean(co2.taj.nv.svsp$TajimaD,na.rm=T)
sd(co2.taj.nv.svsp$TajimaD,na.rm=T)

mean(cv1.taj.svsp$TajimaD,na.rm=T)
sd(cv1.taj.svsp$TajimaD,na.rm=T)
mean(cv2.taj.svsp$TajimaD,na.rm=T)
sd(cv2.taj.svsp$TajimaD,na.rm=T)
mean(co1.taj.svsp$TajimaD,na.rm=T)
sd(co1.taj.svsp$TajimaD,na.rm=T)
mean(co2.taj.svsp$TajimaD,na.rm=T)
sd(co2.taj.svsp$TajimaD,na.rm=T)

## T tests
# Genome-wide
t.test(cv1.taj.10kb$TajimaD,cv1.taj.svsp$TajimaD)
t.test(cv2.taj.10kb$TajimaD,cv2.taj.svsp$TajimaD)
t.test(co1.taj.10kb$TajimaD,co1.taj.svsp$TajimaD)
t.test(co2.taj.10kb$TajimaD,co2.taj.svsp$TajimaD)
# Chromosome 10
t.test(cv1.taj.mi2.bk$TajimaD,cv1.taj.svsp$TajimaD)
t.test(cv2.taj.mi2.bk$TajimaD,cv2.taj.svsp$TajimaD)
t.test(co1.taj.mi2.bk$TajimaD,co1.taj.svsp$TajimaD)
t.test(co2.taj.mi2.bk$TajimaD,co2.taj.svsp$TajimaD)
# Non-venom paralogs
t.test(cv1.taj.nv.svsp$TajimaD,cv1.taj.svsp$TajimaD)
t.test(cv2.taj.nv.svsp$TajimaD,cv2.taj.svsp$TajimaD)
t.test(co1.taj.nv.svsp$TajimaD,co1.taj.svsp$TajimaD)
t.test(co2.taj.nv.svsp$TajimaD,co2.taj.svsp$TajimaD)

## PLA2
mean(cv1.taj.1kb$TajimaD,na.rm=T)
sd(cv1.taj.1kb$TajimaD,na.rm=T)
mean(cv2.taj.1kb$TajimaD,na.rm=T)
sd(cv2.taj.1kb$TajimaD,na.rm=T)
mean(co1.taj.1kb$TajimaD,na.rm=T)
sd(co1.taj.1kb$TajimaD,na.rm=T)
mean(co2.taj.1kb$TajimaD,na.rm=T)
sd(co2.taj.1kb$TajimaD,na.rm=T)

mean(cv1.taj.mi7.1kb$TajimaD,na.rm=T)
sd(cv1.taj.mi7.1kb$TajimaD,na.rm=T)
mean(cv2.taj.mi7.1kb$TajimaD,na.rm=T)
sd(cv2.taj.mi7.1kb$TajimaD,na.rm=T)
mean(co1.taj.mi7.1kb$TajimaD,na.rm=T)
sd(co1.taj.mi7.1kb$TajimaD,na.rm=T)
mean(co2.taj.mi7.1kb$TajimaD,na.rm=T)
sd(co2.taj.mi7.1kb$TajimaD,na.rm=T)

mean(cv1.taj.nv.pla2$TajimaD,na.rm=T)
sd(cv1.taj.nv.pla2$TajimaD,na.rm=T)
mean(cv2.taj.nv.pla2$TajimaD,na.rm=T)
sd(cv2.taj.nv.pla2$TajimaD,na.rm=T)
mean(co1.taj.nv.pla2$TajimaD,na.rm=T)
sd(co1.taj.nv.pla2$TajimaD,na.rm=T)
mean(co2.taj.nv.pla2$TajimaD,na.rm=T)
sd(co2.taj.nv.pla2$TajimaD,na.rm=T)

mean(cv1.taj.pla2$TajimaD,na.rm=T)
sd(cv1.taj.pla2$TajimaD,na.rm=T)
mean(cv2.taj.pla2$TajimaD,na.rm=T)
sd(cv2.taj.pla2$TajimaD,na.rm=T)
mean(co1.taj.pla2$TajimaD,na.rm=T)
sd(co1.taj.pla2$TajimaD,na.rm=T)
mean(co2.taj.pla2$TajimaD,na.rm=T)
sd(co2.taj.pla2$TajimaD,na.rm=T)

# T tests
# Genome-wide
t.test(cv1.taj.1kb$TajimaD,cv1.taj.pla2$TajimaD)
t.test(cv2.taj.1kb$TajimaD,cv2.taj.pla2$TajimaD)
t.test(co1.taj.1kb$TajimaD,co1.taj.pla2$TajimaD)
t.test(co2.taj.1kb$TajimaD,co2.taj.pla2$TajimaD)
# Chromosome 10
t.test(cv1.taj.mi7.bk$TajimaD,cv1.taj.pla2$TajimaD)
t.test(cv2.taj.mi7.bk$TajimaD,cv2.taj.pla2$TajimaD)
t.test(co1.taj.mi7.bk$TajimaD,co1.taj.pla2$TajimaD)
t.test(co2.taj.mi7.bk$TajimaD,co2.taj.pla2$TajimaD)
# Non-venom paralogs
t.test(cv1.taj.nv.pla2$TajimaD,cv1.taj.pla2$TajimaD)
t.test(cv2.taj.nv.pla2$TajimaD,cv2.taj.pla2$TajimaD)
t.test(co1.taj.nv.pla2$TajimaD,co1.taj.pla2$TajimaD)
t.test(co2.taj.nv.pla2$TajimaD,co2.taj.pla2$TajimaD)

### 1.4 Plot----------------------------------------------------------------

par(mfrow=c(3,2))
boxplot(cv1.taj.mi1.bk$TajimaD,cv1.taj.nv.svmp$TajimaD,cv1.taj.svmp$TajimaD,border=c('black','black','lightblue3'),names=c('Chr 9','NV','SVMP'),outline=F,ylab="Tajima's D")
boxplot(co1.taj.mi1.bk$TajimaD,co1.taj.nv.svmp$TajimaD,co1.taj.svmp$TajimaD,border=c('black','black','lightblue3'),names=c('Chr 9','NV','SVMP'),outline=F,ylab="Tajima's D")
boxplot(cv1.taj.mi2.bk$TajimaD,cv1.taj.nv.svsp$TajimaD,cv1.taj.svsp$TajimaD,border=c('black','black','aquamarine3'),names=c('Chr 10','NV','SVSP'),outline=F,ylab="Tajima's D")
boxplot(co1.taj.mi2.bk$TajimaD,co1.taj.nv.svsp$TajimaD,co1.taj.svsp$TajimaD,border=c('black','black','aquamarine3'),names=c('Chr 10','NV','SVSP'),outline=F,ylab="Tajima's D")
boxplot(cv1.taj.mi7.bk$TajimaD,cv1.taj.nv.pla2$TajimaD,cv1.taj.pla2$TajimaD,border=c('black','black','maroon'),names=c('Chr 15','NV','PLA2'),outline=F,ylab="Tajima's D")
boxplot(co1.taj.mi7.bk$TajimaD,co1.taj.nv.pla2$TajimaD,co1.taj.pla2$TajimaD,border=c('black','black','maroon'),names=c('Chr 15','NV','PLA2'),outline=F,ylab="Tajima's D")

par(mfrow=c(3,2))
boxplot(cv2.taj.mi1.bk$TajimaD,cv2.taj.nv.svmp$TajimaD,cv2.taj.svmp$TajimaD,border=c('black','black','lightblue3'),names=c('Chr 9','NV','SVMP'),outline=F,ylab="Tajima's D")
boxplot(co2.taj.mi1.bk$TajimaD,co2.taj.nv.svmp$TajimaD,co2.taj.svmp$TajimaD,border=c('black','black','lightblue3'),names=c('Chr 9','NV','SVMP'),outline=F,ylab="Tajima's D")
boxplot(cv2.taj.mi2.bk$TajimaD,cv2.taj.nv.svsp$TajimaD,cv2.taj.svsp$TajimaD,border=c('black','black','aquamarine3'),names=c('Chr 10','NV','SVSP'),outline=F,ylab="Tajima's D")
boxplot(co2.taj.mi2.bk$TajimaD,co2.taj.nv.svsp$TajimaD,co2.taj.svsp$TajimaD,border=c('black','black','aquamarine3'),names=c('Chr 10','NV','SVSP'),outline=F,ylab="Tajima's D")
boxplot(cv2.taj.mi7.bk$TajimaD,cv2.taj.nv.pla2$TajimaD,cv2.taj.pla2$TajimaD,border=c('black','black','maroon'),names=c('Chr 15','NV','PLA2'),outline=F,ylab="Tajima's D")
boxplot(co2.taj.mi7.bk$TajimaD,co2.taj.nv.pla2$TajimaD,co2.taj.pla2$TajimaD,border=c('black','black','maroon'),names=c('Chr 15','NV','PLA2'),outline=F,ylab="Tajima's D")

### 2. df--------------------------------------------------------------------

setwd('./df')

### 2.1 Read in and parse df results-----------------------------------------

# Genome-wide 10/1 kb windows
cv1co1.df.10kb <- read.table('window.10kb.df_prop.all.cv1.co1.cnvMask.txt',header=T)
cv1cv2.df.10kb <- read.table('window.10kb.df_prop.all.cv1.cv2.txt',header=T)
co1co2.df.10kb <- read.table('window.10kb.df_prop.all.co1.co2.cnvMask.txt',header=T)

cv1co1.df.1kb <- read.table('window.1kb.df_prop.all.cv1.co1.cnvMask.txt',header=T)
cv1cv2.df.1kb <- read.table('window.1kb.df_prop.all.cv1.cv2.txt',header=T)
co1co2.df.1kb <- read.table('window.1kb.df_prop.all.co1.co2.cnvMask.txt',header=T)

# Venom chromosome 10/1 kb windows
cv1co1.df.mi1.1kb <- cv1co1.df.1kb[which(cv1co1.df.1kb$chrom=='scaffold-mi1'),]
cv1cv2.df.mi1.1kb <- cv1cv2.df.1kb[which(cv1cv2.df.1kb$chrom=='scaffold-mi1'),]
co1co2.df.mi1.1kb <- co1co2.df.1kb[which(co1co2.df.1kb$chrom=='scaffold-mi1'),]

cv1co1.df.mi2.1kb <- cv1co1.df.1kb[which(cv1co1.df.1kb$chrom=='scaffold-mi2'),]
cv1cv2.df.mi2.1kb <- cv1cv2.df.1kb[which(cv1cv2.df.1kb$chrom=='scaffold-mi2'),]
co1co2.df.mi2.1kb <- co1co2.df.1kb[which(co1co2.df.1kb$chrom=='scaffold-mi2'),]

cv1co1.df.mi7.10kb <- cv1co1.df.10kb[which(cv1co1.df.10kb$chrom=='scaffold-mi7'),]
cv1cv2.df.mi7.10kb <- cv1cv2.df.10kb[which(cv1cv2.df.10kb$chrom=='scaffold-mi7'),]
co1co2.df.mi7.10kb <- co1co2.df.10kb[which(co1co2.df.10kb$chrom=='scaffold-mi7'),]

cv1co1.df.mi7.1kb <- cv1co1.df.1kb[which(cv1co1.df.1kb$chrom=='scaffold-mi7'),]
cv1cv2.df.mi7.1kb <- cv1cv2.df.1kb[which(cv1cv2.df.1kb$chrom=='scaffold-mi7'),]
co1co2.df.mi7.1kb <- co1co2.df.1kb[which(co1co2.df.1kb$chrom=='scaffold-mi7'),]

# Parse venom regions
cv1co1.df.svmp <- cv1co1.df.mi1.1kb[which(cv1co1.df.mi1.1kb$start>=svmp.reg$V2 & cv1co1.df.mi1.1kb$end<=svmp.reg$V3),]
cv1cv2.df.svmp <- cv1cv2.df.mi1.1kb[which(cv1cv2.df.mi1.1kb$start>=svmp.reg$V2 & cv1cv2.df.mi1.1kb$end<=svmp.reg$V3),]
co1co2.df.svmp <- co1co2.df.mi1.1kb[which(co1co2.df.mi1.1kb$start>=svmp.reg$V2 & co1co2.df.mi1.1kb$end<=svmp.reg$V3),]

cv1co1.df.svsp <- cv1co1.df.mi2.1kb[which(cv1co1.df.mi2.1kb$start>=svsp.reg$V2 & cv1co1.df.mi2.1kb$end<=svsp.reg$V3),]
cv1cv2.df.svsp <- cv1cv2.df.mi2.1kb[which(cv1cv2.df.mi2.1kb$start>=svsp.reg$V2 & cv1cv2.df.mi2.1kb$end<=svsp.reg$V3),]
co1co2.df.svsp <- co1co2.df.mi2.1kb[which(co1co2.df.mi2.1kb$start>=svsp.reg$V2 & co1co2.df.mi2.1kb$end<=svsp.reg$V3),]

cv1co1.df.pla2 <- cv1co1.df.mi7.1kb[which(cv1co1.df.mi7.1kb$start>=pla2.reg$V2 & cv1co1.df.mi7.1kb$end<=pla2.reg$V3),]
cv1cv2.df.pla2 <- cv1cv2.df.mi7.1kb[which(cv1cv2.df.mi7.1kb$start>=pla2.reg$V2 & cv1cv2.df.mi7.1kb$end<=pla2.reg$V3),]
co1co2.df.pla2 <- co1co2.df.mi7.1kb[which(co1co2.df.mi7.1kb$start>=pla2.reg$V2 & co1co2.df.mi7.1kb$end<=pla2.reg$V3),]

# Parse chromosome backgrounds
cv1co1.df.mi1.bk <- cv1co1.df.mi1.1kb[which(cv1co1.df.mi1.1kb$start<svmp.reg$V2 & cv1co1.df.mi1.1kb$end<svmp.reg$V2 | cv1co1.df.mi1.1kb$start>svmp.reg$V3 & cv1co1.df.mi1.1kb$end>svmp.reg$V3),]
cv1cv2.df.mi1.bk <- cv1cv2.df.mi1.1kb[which(cv1cv2.df.mi1.1kb$start<svmp.reg$V2 & cv1cv2.df.mi1.1kb$end<svmp.reg$V2 | cv1cv2.df.mi1.1kb$start>svmp.reg$V3 & cv1cv2.df.mi1.1kb$end>svmp.reg$V3),]
co1co2.df.mi1.bk <- co1co2.df.mi1.1kb[which(co1co2.df.mi1.1kb$start<svmp.reg$V2 & co1co2.df.mi1.1kb$end<svmp.reg$V2 | co1co2.df.mi1.1kb$start>svmp.reg$V3 & co1co2.df.mi1.1kb$end>svmp.reg$V3),]

cv1co1.df.mi2.bk <- cv1co1.df.mi2.1kb[which(cv1co1.df.mi2.1kb$start<svsp.reg$V2 & cv1co1.df.mi2.1kb$end<svsp.reg$V2 | cv1co1.df.mi2.1kb$start>svsp.reg$V3 & cv1co1.df.mi2.1kb$end>svsp.reg$V3),]
cv1cv2.df.mi2.bk <- cv1cv2.df.mi2.1kb[which(cv1cv2.df.mi2.1kb$start<svsp.reg$V2 & cv1cv2.df.mi2.1kb$end<svsp.reg$V2 | cv1cv2.df.mi2.1kb$start>svsp.reg$V3 & cv1cv2.df.mi2.1kb$end>svsp.reg$V3),]
co1co2.df.mi2.bk <- co1co2.df.mi2.1kb[which(co1co2.df.mi2.1kb$start<svsp.reg$V2 & co1co2.df.mi2.1kb$end<svsp.reg$V2 | co1co2.df.mi2.1kb$start>svsp.reg$V3 & co1co2.df.mi2.1kb$end>svsp.reg$V3),]

cv1co1.df.mi7.bk <- cv1co1.df.mi7.1kb[which(cv1co1.df.mi7.1kb$start<pla2.reg$V2 & cv1co1.df.mi7.1kb$end<pla2.reg$V2 | cv1co1.df.mi7.1kb$start>pla2.reg$V3 & cv1co1.df.mi7.1kb$end>pla2.reg$V3),]
cv1cv2.df.mi7.bk <- cv1cv2.df.mi7.1kb[which(cv1cv2.df.mi7.1kb$start<pla2.reg$V2 & cv1cv2.df.mi7.1kb$end<pla2.reg$V2 | cv1cv2.df.mi7.1kb$start>pla2.reg$V3 & cv1cv2.df.mi7.1kb$end>pla2.reg$V3),]
co1co2.df.mi7.bk <- co1co2.df.mi7.1kb[which(co1co2.df.mi7.1kb$start<pla2.reg$V2 & co1co2.df.mi7.1kb$end<pla2.reg$V2 | co1co2.df.mi7.1kb$start>pla2.reg$V3 & co1co2.df.mi7.1kb$end>pla2.reg$V3),]

# Read in non-venom paralog backgrounds

cv1co1.df.nv.svmp <- read.table('window.1kb.df_prop.non-venom_SVMP.cv1.co1.txt',header=T)
cv1cv2.df.nv.svmp <- read.table('window.1kb.df_prop.non-venom_SVMP.cv1.cv2.txt',header=T)
co1co2.df.nv.svmp <- read.table('window.1kb.df_prop.non-venom_SVMP.co1.co2.txt',header=T)

cv1co1.df.nv.svsp <- read.table('window.1kb.df_prop.non-venom_SVSP.cv1.co1.txt',header=T)
cv1cv2.df.nv.svsp <- read.table('window.1kb.df_prop.non-venom_SVSP.cv1.cv2.txt',header=T)
co1co2.df.nv.svsp <- read.table('window.1kb.df_prop.non-venom_SVSP.co1.co2.txt',header=T)

cv1co1.df.nv.pla2 <- read.table('window.1kb.df_prop.non-venom_PLA2.cv1.co1.txt',header=T)
cv1cv2.df.nv.pla2 <- read.table('window.1kb.df_prop.non-venom_PLA2.cv1.cv2.txt',header=T)
co1co2.df.nv.pla2 <- read.table('window.1kb.df_prop.non-venom_PLA2.co1.co2.txt',header=T)

### 2.2 Calculate genome-wide summary statistics----------------------------

mean(cv1co1.df.10kb$df,na.rm=T)
mean(cv1co1.df.mi1.10kb$df,na.rm=T)
mean(cv1co1.df.nv.svmp$df)
mean(cv1co1.df.svmp$df,na.rm=T)

sd(cv1co1.df.10kb$df,na.rm=T)
sd(cv1co1.df.mi1.10kb$df,na.rm=T)
sd(cv1co1.df.nv.svmp$df)
sd(cv1co1.df.svmp$df,na.rm=T)

mean(cv1co1.df.10kb$df,na.rm=T)
mean(cv1co1.df.mi2.10kb$df,na.rm=T)
mean(cv1co1.df.nv.svsp$df,na.rm=T)
mean(cv1co1.df.svsp$df,na.rm=T)

sd(cv1co1.df.10kb$df,na.rm=T)
sd(cv1co1.df.mi2.10kb$df,na.rm=T)
sd(cv1co1.df.nv.svsp$df,na.rm=T)
sd(cv1co1.df.svsp$df,na.rm=T)

mean(cv1co1.df.1kb$df,na.rm=T)
mean(cv1co1.df.mi7.1kb$df,na.rm=T)
mean(cv1co1.df.nv.pla2$df,na.rm=T)
mean(cv1co1.df.pla2$df,na.rm=T)

sd(cv1co1.df.1kb$df,na.rm=T)
sd(cv1co1.df.mi7.1kb$df,na.rm=T)
sd(cv1co1.df.nv.pla2$df,na.rm=T)
sd(cv1co1.df.pla2$df,na.rm=T)

### 2.3 Statistical comparison of venom df vs backgrounds-------------------

wilcox.test(cv1co1.df.10kb$df,cv1co1.df.svmp$df)
wilcox.test(cv1co1.df.mi1.bk$df,cv1co1.df.svmp$df)
wilcox.test(cv1co1.df.nv.svmp$df,cv1co1.df.svmp$df)

wilcox.test(cv1co1.df.10kb$df,cv1co1.df.svsp$df)
wilcox.test(cv1co1.df.mi2.bk$df,cv1co1.df.svsp$df)
wilcox.test(cv1co1.df.nv.svsp$df,cv1co1.df.svsp$df)

wilcox.test(cv1co1.df.1kb$df,cv1co1.df.pla2$df)
wilcox.test(cv1co1.df.mi7.bk$df,cv1co1.df.pla2$df)
wilcox.test(cv1co1.df.nv.pla2$df,cv1co1.df.pla2$df)

### 2.4 Plot--------------------------------------------------------------------

par(mfrow=c(3,2))
boxplot(cv1co1.df.mi1.bk$df,cv1co1.df.nv.svmp$df,cv1co1.df.svmp$df,border=c('black','black','lightblue3'),names=c('Chr. 9','NV','SVMP'),outline=F,ylab="df")
boxplot(cv1co1.df.mi2.bk$df,cv1co1.df.nv.svsp$df,cv1co1.df.svsp$df,border=c('black','black','aquamarine3'),names=c('Chr. 10','NV','SVSP'),outline=F,ylab="df")
boxplot(cv1co1.df.mi7.bk$df,cv1co1.df.nv.pla2$df,cv1co1.df.pla2$df,border=c('black','black','maroon'),names=c('Chr. 15','NV','PLA2'),outline=F,ylab="df")


### 3. iHS------------------------------------------------------------------

setwd('./rehh/')

### 3.1 Read in and parse data----------------------------------------------

# Genome-wide windows
cv1.ihs.100kb <- read.table('cv.all_ihs.100kb.txt',header=T)
co1.ihs.100kb <- read.table('co.all_ihs.100kb.txt',header=T)
cv1.ihs.10kb <- read.table('cv.all_ihs.10kb.txt',header=T)
co1.ihs.10kb <- read.table('co.all_ihs.10kb.cnvMask.txt',header=T)
cv1.ihs.1kb <- read.table('cv.all_ihs.1kb.txt',header=T)
co1.ihs.1kb <- read.table('co.all_ihs.1kb.cnvMask.txt',header=T)

# Venom chromosome 10/1 kb windows
cv1.ihs.mi1.10kb <- cv1.ihs.10kb[which(cv1.ihs.10kb$chrom=='scaffold-mi1'),]
co1.ihs.mi1.10kb <- co1.ihs.10kb[which(co1.ihs.10kb$chrom=='scaffold-mi1'),]
cv1.ihs.mi2.10kb <- cv1.ihs.10kb[which(cv1.ihs.10kb$chrom=='scaffold-mi2'),]
co1.ihs.mi2.10kb <- co1.ihs.10kb[which(co1.ihs.10kb$chrom=='scaffold-mi2'),]
cv1.ihs.mi7.10kb <- cv1.ihs.10kb[which(cv1.ihs.10kb$chrom=='scaffold-mi7'),]
co1.ihs.mi7.10kb <- co1.ihs.10kb[which(co1.ihs.10kb$chrom=='scaffold-mi7'),]
cv1.ihs.mi7.1kb <- cv1.ihs.1kb[which(cv1.ihs.1kb$chrom=='scaffold-mi7'),]
co1.ihs.mi7.1kb <- co1.ihs.1kb[which(co1.ihs.1kb$chrom=='scaffold-mi7'),]

# Parse venom regions
cv1.ihs.svmp <- cv1.ihs.mi1.10kb[which(cv1.ihs.mi1.10kb$start>=svmp.reg$V2 & cv1.ihs.mi1.10kb$end<=svmp.reg$V3),]
co1.ihs.svmp <- co1.ihs.mi1.10kb[which(co1.ihs.mi1.10kb$start>=svmp.reg$V2 & co1.ihs.mi1.10kb$end<=svmp.reg$V3),]
cv1.ihs.svsp <- cv1.ihs.mi2.10kb[which(cv1.ihs.mi2.10kb$start>=svsp.reg$V2 & cv1.ihs.mi2.10kb$end<=svsp.reg$V3),]
co1.ihs.svsp <- co1.ihs.mi2.10kb[which(co1.ihs.mi2.10kb$start>=svsp.reg$V2 & co1.ihs.mi2.10kb$end<=svsp.reg$V3),]
cv1.ihs.pla2 <- cv1.ihs.mi7.1kb[which(cv1.ihs.mi7.1kb$start>=pla2.reg$V2 & cv1.ihs.mi7.1kb$end<=pla2.reg$V3),]
co1.ihs.pla2 <- co1.ihs.mi7.1kb[which(co1.ihs.mi7.1kb$start>=pla2.reg$V2 & co1.ihs.mi7.1kb$end<=pla2.reg$V3),]

# Parse chromosome backgrounds
cv1.ihs.mi1.bk <- cv1.ihs.mi1.10kb[which(cv1.ihs.mi1.10kb$start<svmp.reg$V2 & cv1.ihs.mi1.10kb$end<svmp.reg$V2 | cv1.ihs.mi1.10kb$start>svmp.reg$V3 & cv1.ihs.mi1.10kb$end>svmp.reg$V3),]
co1.ihs.mi1.bk <- co1.ihs.mi1.10kb[which(co1.ihs.mi1.10kb$start<svmp.reg$V2 & co1.ihs.mi1.10kb$end<svmp.reg$V2 | co1.ihs.mi1.10kb$start>svmp.reg$V3 & co1.ihs.mi1.10kb$end>svmp.reg$V3),]
cv1.ihs.mi2.bk <- cv1.ihs.mi2.10kb[which(cv1.ihs.mi2.10kb$start<svmp.reg$V2 & cv1.ihs.mi2.10kb$end<svmp.reg$V2 | cv1.ihs.mi2.10kb$start>svmp.reg$V3 & cv1.ihs.mi2.10kb$end>svmp.reg$V3),]
co1.ihs.mi2.bk <- co1.ihs.mi2.10kb[which(co1.ihs.mi2.10kb$start<svmp.reg$V2 & co1.ihs.mi2.10kb$end<svmp.reg$V2 | co1.ihs.mi2.10kb$start>svmp.reg$V3 & co1.ihs.mi2.10kb$end>svmp.reg$V3),]
cv1.ihs.mi7.bk <- cv1.ihs.mi7.10kb[which(cv1.ihs.mi7.10kb$start<svmp.reg$V2 & cv1.ihs.mi7.10kb$end<svmp.reg$V2 | cv1.ihs.mi7.10kb$start>svmp.reg$V3 & cv1.ihs.mi7.10kb$end>svmp.reg$V3),]
co1.ihs.mi7.bk <- co1.ihs.mi7.10kb[which(co1.ihs.mi7.10kb$start<svmp.reg$V2 & co1.ihs.mi7.10kb$end<svmp.reg$V2 | co1.ihs.mi7.10kb$start>svmp.reg$V3 & co1.ihs.mi7.10kb$end>svmp.reg$V3),]

# Read in non-venom paralog backgrounds
cv1.ihs.nv.svmp <- read.table('cv.non-venom_SVMP_ihs.10kb.txt',header=T)
co1.ihs.nv.svmp <- read.table('co.non-venom_SVMP_ihs.10kb.txt',header=T)
cv1.ihs.nv.svsp <- read.table('cv.non-venom_SVSP_ihs.10kb.txt',header=T)
co1.ihs.nv.svsp <- read.table('co.non-venom_SVSP_ihs.10kb.txt',header=T)
cv1.ihs.nv.pla2 <- read.table('cv.non-venom_PLA2_ihs.1kb.txt',header=T)
co1.ihs.nv.pla2 <- read.table('co.non-venom_PLA2_ihs.1kb.txt',header=T)

### 3.2 Calculate summary statistics----------------------------------------

# SVMP
mean(as.numeric(cv1.ihs.10kb$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.mi1.10kb$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.nv.svmp$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.svmp$iHS),na.rm=T)

sd(as.numeric(cv1.ihs.10kb$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.mi1.10kb$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.nv.svmp$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.svmp$iHS),na.rm=T)

mean(as.numeric(co1.ihs.10kb$iHS),na.rm=T)
mean(as.numeric(co1.ihs.mi1.10kb$iHS),na.rm=T)
mean(as.numeric(co1.ihs.nv.svmp$iHS),na.rm=T)
mean(as.numeric(co1.ihs.svmp$iHS),na.rm=T)

sd(as.numeric(co1.ihs.10kb$iHS),na.rm=T)
sd(as.numeric(co1.ihs.mi1.10kb$iHS),na.rm=T)
sd(as.numeric(co1.ihs.nv.svmp$iHS),na.rm=T)
sd(as.numeric(co1.ihs.svmp$iHS),na.rm=T)

# SVSP
mean(as.numeric(cv1.ihs.10kb$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.mi2.10kb$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.nv.svsp$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.svsp$iHS),na.rm=T)

sd(as.numeric(cv1.ihs.10kb$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.mi2.10kb$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.nv.svsp$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.svsp$iHS),na.rm=T)

mean(as.numeric(co1.ihs.10kb$iHS),na.rm=T)
mean(as.numeric(co1.ihs.mi2.10kb$iHS),na.rm=T)
mean(as.numeric(co1.ihs.nv.svsp$iHS),na.rm=T)
mean(as.numeric(co1.ihs.svsp$iHS),na.rm=T)

sd(as.numeric(co1.ihs.10kb$iHS),na.rm=T)
sd(as.numeric(co1.ihs.mi2.10kb$iHS),na.rm=T)
sd(as.numeric(co1.ihs.nv.svsp$iHS),na.rm=T)
sd(as.numeric(co1.ihs.svsp$iHS),na.rm=T)

# PLA2
mean(as.numeric(cv1.ihs.10kb$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.mi7.10kb$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.nv.pla2$iHS),na.rm=T)
mean(as.numeric(cv1.ihs.pla2$iHS),na.rm=T)

sd(as.numeric(cv1.ihs.10kb$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.mi7.10kb$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.nv.pla2$iHS),na.rm=T)
sd(as.numeric(cv1.ihs.pla2$iHS),na.rm=T)

mean(as.numeric(co1.ihs.10kb$iHS),na.rm=T)
mean(as.numeric(co1.ihs.mi7.10kb$iHS),na.rm=T)
mean(as.numeric(co1.ihs.nv.pla2$iHS),na.rm=T)
mean(as.numeric(co1.ihs.pla2$iHS),na.rm=T)

sd(as.numeric(co1.ihs.10kb$iHS),na.rm=T)
sd(as.numeric(co1.ihs.mi7.10kb$iHS),na.rm=T)
sd(as.numeric(co1.ihs.nv.pla2$iHS),na.rm=T)
sd(as.numeric(co1.ihs.pla2$iHS),na.rm=T)

# Quantiles
quantile(as.numeric(cv1.ihs.100kb$iHS),c(0.95,0.975,0.99),na.rm=T)
quantile(as.numeric(co1.ihs.100kb$iHS),c(0.95,0.975,0.99),na.rm=T)

### 3.3 Statistical comparisons of venom iHS vs backgrounds-----------------

t.test(as.numeric(cv1.ihs.10kb$iHS),as.numeric(cv1.ihs.svmp$iHS))
t.test(as.numeric(co1.ihs.10kb$iHS),as.numeric(co1.ihs.svmp$iHS))
t.test(as.numeric(cv1.ihs.mi1.bk$iHS),as.numeric(cv1.ihs.svmp$iHS))
t.test(as.numeric(co1.ihs.mi1.bk$iHS),as.numeric(co1.ihs.svmp$iHS))
t.test(as.numeric(cv1.ihs.nv.svmp$iHS),as.numeric(cv1.ihs.svmp$iHS))
t.test(as.numeric(co1.ihs.nv.svmp$iHS),as.numeric(co1.ihs.svmp$iHS))

t.test(as.numeric(cv1.ihs.10kb$iHS),as.numeric(cv1.ihs.svsp$iHS))
t.test(as.numeric(co1.ihs.10kb$iHS),as.numeric(co1.ihs.svsp$iHS))
t.test(as.numeric(cv1.ihs.mi2.bk$iHS),as.numeric(cv1.ihs.svsp$iHS))
t.test(as.numeric(co1.ihs.mi2.bk$iHS),as.numeric(co1.ihs.svsp$iHS))
t.test(as.numeric(cv1.ihs.nv.svsp$iHS),as.numeric(cv1.ihs.svsp$iHS))
t.test(as.numeric(co1.ihs.nv.svsp$iHS),as.numeric(co1.ihs.svsp$iHS))

t.test(as.numeric(cv1.ihs.10kb$iHS),as.numeric(cv1.ihs.pla2$iHS))
t.test(as.numeric(co1.ihs.10kb$iHS),as.numeric(co1.ihs.pla2$iHS))
t.test(as.numeric(cv1.ihs.mi7.bk$iHS),as.numeric(cv1.ihs.pla2$iHS))
t.test(as.numeric(co1.ihs.mi7.bk$iHS),as.numeric(co1.ihs.pla2$iHS))
t.test(as.numeric(cv1.ihs.nv.pla2$iHS),as.numeric(cv1.ihs.pla2$iHS))
t.test(as.numeric(co1.ihs.nv.pla2$iHS),as.numeric(co1.ihs.pla2$iHS))

### 3.4 Plot----------------------------------------------------------------

par(mfrow=c(3,2))
boxplot(as.numeric(cv1.ihs.mi1.bk$iHS),as.numeric(cv1.ihs.nv.svmp$iHS),as.numeric(cv1.ihs.svmp$iHS),border=c('black','black','lightblue3'),names=c('Chr. 9','NV','SVMP'),outline=F,ylab="ihs",ylim=c(-1,2))
boxplot(as.numeric(co1.ihs.mi1.bk$iHS),as.numeric(co1.ihs.nv.svmp$iHS),as.numeric(co1.ihs.svmp$iHS),border=c('black','black','lightblue3'),names=c('Chr. 9','NV','SVMP'),outline=F,ylab="ihs",ylim=c(-1,2))
boxplot(as.numeric(cv1.ihs.mi2.bk$iHS),as.numeric(cv1.ihs.nv.svsp$iHS),as.numeric(cv1.ihs.svsp$iHS),border=c('black','black','aquamarine3'),names=c('Chr. 10','NV','SVSP'),outline=F,ylab="ihs",ylim=c(-1,2))
boxplot(as.numeric(co1.ihs.mi2.bk$iHS),as.numeric(co1.ihs.nv.svsp$iHS),as.numeric(co1.ihs.svsp$iHS),border=c('black','black','aquamarine3'),names=c('Chr. 10','NV','SVSP'),outline=F,ylab="ihs",ylim=c(-1,2))
boxplot(as.numeric(cv1.ihs.mi7.bk$iHS),as.numeric(cv1.ihs.nv.pla2$iHS),as.numeric(cv1.ihs.pla2$iHS),border=c('black','black','maroon'),names=c('Chr. 15','NV','PLA2'),outline=F,ylab="ihs",ylim=c(-1,2))
boxplot(as.numeric(co1.ihs.mi7.bk$iHS),as.numeric(co1.ihs.nv.pla2$iHS),as.numeric(co1.ihs.pla2$iHS),border=c('black','black','maroon'),names=c('Chr. 15','NV','PLA2'),outline=F,ylab="ihs",ylim=c(-1,2))


### 4. ß--------------------------------------------------------------------

setwd('./beta/results/')

### 4.1 Read in and parse data----------------------------------------------

# Genome-wide 10 kb windows
cv1.beta.10kb <- read.table('cv1.phased.all.betascores.10kb.txt',header=T)
co1.beta.10kb <- read.table('co1.phased.all.betascores.10kb.cnvMask.txt',header=T)

cv1.beta.1kb <- read.table('cv1.phased.all.betascores.1kb.txt',header=T)
co1.beta.1kb <- read.table('co1.phased.all.betascores.1kb.cnvMask.txt',header=T)

# Venom chromosome 10/1 kb windows
cv1.beta.mi1.10kb <- cv1.beta.10kb[which(cv1.beta.10kb$chrom=='scaffold-mi1'),]
co1.beta.mi1.10kb <- co1.beta.10kb[which(co1.beta.10kb$chrom=='scaffold-mi1'),]
cv1.beta.mi2.10kb <- cv1.beta.10kb[which(cv1.beta.10kb$chrom=='scaffold-mi2'),]
co1.beta.mi2.10kb <- co1.beta.10kb[which(co1.beta.10kb$chrom=='scaffold-mi2'),]
cv1.beta.mi7.10kb <- cv1.beta.10kb[which(cv1.beta.10kb$chrom=='scaffold-mi7'),]
co1.beta.mi7.10kb <- co1.beta.10kb[which(co1.beta.10kb$chrom=='scaffold-mi7'),]
cv1.beta.mi7.1kb <- cv1.beta.1kb[which(cv1.beta.1kb$chrom=='scaffold-mi7'),]
co1.beta.mi7.1kb <- co1.beta.1kb[which(co1.beta.1kb$chrom=='scaffold-mi7'),]

# Parse venom regions
cv1.beta.svmp <- cv1.beta.mi1.10kb[which(cv1.beta.mi1.10kb$start>=svmp.reg$V2 & cv1.beta.mi1.10kb$end<=svmp.reg$V3),]
co1.beta.svmp <- co1.beta.mi1.10kb[which(co1.beta.mi1.10kb$start>=svmp.reg$V2 & co1.beta.mi1.10kb$end<=svmp.reg$V3),]
cv1.beta.svsp <- cv1.beta.mi2.10kb[which(cv1.beta.mi2.10kb$start>=svsp.reg$V2 & cv1.beta.mi2.10kb$end<=svsp.reg$V3),]
co1.beta.svsp <- co1.beta.mi2.10kb[which(co1.beta.mi2.10kb$start>=svsp.reg$V2 & co1.beta.mi2.10kb$end<=svsp.reg$V3),]
cv1.beta.pla2 <- cv1.beta.mi7.1kb[which(cv1.beta.mi7.1kb$start>=pla2.reg$V2 & cv1.beta.mi7.1kb$end<=pla2.reg$V3),]
co1.beta.pla2 <- co1.beta.mi7.1kb[which(co1.beta.mi7.1kb$start>=pla2.reg$V2 & co1.beta.mi7.1kb$end<=pla2.reg$V3),]

# Parse chromosome backgrounds
cv1.beta.mi1.bk <- cv1.beta.mi1.10kb[which(cv1.beta.mi1.10kb$start<svmp.reg$V2 & cv1.beta.mi1.10kb$end<svmp.reg$V2 | cv1.beta.mi1.10kb$start>svmp.reg$V3 & cv1.beta.mi1.10kb$end>svmp.reg$V3),]
co1.beta.mi1.bk <- co1.beta.mi1.10kb[which(co1.beta.mi1.10kb$start<svmp.reg$V2 & co1.beta.mi1.10kb$end<svmp.reg$V2 | co1.beta.mi1.10kb$start>svmp.reg$V3 & co1.beta.mi1.10kb$end>svmp.reg$V3),]
cv1.beta.mi2.bk <- cv1.beta.mi2.10kb[which(cv1.beta.mi2.10kb$start<svmp.reg$V2 & cv1.beta.mi2.10kb$end<svmp.reg$V2 | cv1.beta.mi2.10kb$start>svmp.reg$V3 & cv1.beta.mi2.10kb$end>svmp.reg$V3),]
co1.beta.mi2.bk <- co1.beta.mi2.10kb[which(co1.beta.mi2.10kb$start<svmp.reg$V2 & co1.beta.mi2.10kb$end<svmp.reg$V2 | co1.beta.mi2.10kb$start>svmp.reg$V3 & co1.beta.mi2.10kb$end>svmp.reg$V3),]
cv1.beta.mi7.bk <- cv1.beta.mi7.10kb[which(cv1.beta.mi7.10kb$start<svmp.reg$V2 & cv1.beta.mi7.10kb$end<svmp.reg$V2 | cv1.beta.mi7.10kb$start>svmp.reg$V3 & cv1.beta.mi7.10kb$end>svmp.reg$V3),]
co1.beta.mi7.bk <- co1.beta.mi7.10kb[which(co1.beta.mi7.10kb$start<svmp.reg$V2 & co1.beta.mi7.10kb$end<svmp.reg$V2 | co1.beta.mi7.10kb$start>svmp.reg$V3 & co1.beta.mi7.10kb$end>svmp.reg$V3),]

# Read in non-venom paralog backgrounds
cv1.beta.nv.svmp <- read.table('cv1.phased.non-venom_SVMP.betascores.10kb.txt',header=T)
co1.beta.nv.svmp <- read.table('co1.phased.non-venom_SVMP.betascores.10kb.txt',header=T)
cv1.beta.nv.svsp <- read.table('cv1.phased.non-venom_SVSP.betascores.10kb.txt',header=T)
co1.beta.nv.svsp <- read.table('co1.phased.non-venom_SVSP.betascores.10kb.txt',header=T)
cv1.beta.nv.pla2 <- read.table('cv1.phased.non-venom_PLA2.betascores.1kb.txt',header=T)
co1.beta.nv.pla2 <- read.table('co1.phased.non-venom_PLA2.betascores.1kb.txt',header=T)

### 4.2 Calculated summary statistics---------------------------------------

## SVMP
mean(as.numeric(cv1.beta.svmp$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.svmp$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.svmp$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.svmp$Beta1.),na.rm=T)

mean(as.numeric(cv1.beta.mi1.bk$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.mi1.bk$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.mi1.bk$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.mi1.bk$Beta1.),na.rm=T)

mean(as.numeric(cv1.beta.nv.svmp$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.nv.svmp$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.nv.svmp$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.nv.svmp$Beta1.),na.rm=T)

## SVSP
mean(as.numeric(cv1.beta.svsp$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.svsp$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.svsp$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.svsp$Beta1.),na.rm=T)

mean(as.numeric(cv1.beta.mi2.bk$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.mi2.bk$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.mi2.bk$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.mi2.bk$Beta1.),na.rm=T)

mean(as.numeric(cv1.beta.nv.svsp$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.nv.svsp$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.nv.svsp$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.nv.svsp$Beta1.),na.rm=T)

## PLA2 
mean(as.numeric(cv1.beta.pla2$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.pla2$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.pla2$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.pla2$Beta1.),na.rm=T)

mean(as.numeric(cv1.beta.mi7.bk$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.mi7.bk$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.mi7.bk$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.mi7.bk$Beta1.),na.rm=T)

mean(as.numeric(cv1.beta.nv.pla2$Beta1.),na.rm=T)
mean(as.numeric(co1.beta.nv.pla2$Beta1.),na.rm=T)
sd(as.numeric(cv1.beta.nv.pla2$Beta1.),na.rm=T)
sd(as.numeric(co1.beta.nv.pla2$Beta1.),na.rm=T)

### 4.3 Statistical comparisons of ß in venom regions vs backgrounds--------

wilcox.test(as.numeric(cv1.beta.10kb$Beta1.),as.numeric(cv1.beta.svmp$Beta1.))
wilcox.test(as.numeric(co1.beta.10kb$Beta1.),as.numeric(co1.beta.svmp$Beta1.))
wilcox.test(as.numeric(cv1.beta.mi1.bk$Beta1.),as.numeric(cv1.beta.svmp$Beta1.))
wilcox.test(as.numeric(co1.beta.mi1.bk$Beta1.),as.numeric(co1.beta.svmp$Beta1.))
wilcox.test(as.numeric(cv1.beta.nv.svmp$Beta1.),as.numeric(cv1.beta.svmp$Beta1.))
wilcox.test(as.numeric(co1.beta.nv.svmp$Beta1.),as.numeric(co1.beta.svmp$Beta1.))

wilcox.test(as.numeric(cv1.beta.10kb$Beta1.),as.numeric(cv1.beta.svsp$Beta1.))
wilcox.test(as.numeric(co1.beta.10kb$Beta1.),as.numeric(co1.beta.svsp$Beta1.))
wilcox.test(as.numeric(cv1.beta.mi2.bk$Beta1.),as.numeric(cv1.beta.svsp$Beta1.))
wilcox.test(as.numeric(co1.beta.mi2.bk$Beta1.),as.numeric(co1.beta.svsp$Beta1.))
wilcox.test(as.numeric(cv1.beta.nv.svsp$Beta1.),as.numeric(cv1.beta.svsp$Beta1.))
wilcox.test(as.numeric(co1.beta.nv.svsp$Beta1.),as.numeric(co1.beta.svsp$Beta1.))

wilcox.test(as.numeric(cv1.beta.10kb$Beta1.),as.numeric(cv1.beta.pla2$Beta1.))
wilcox.test(as.numeric(co1.beta.10kb$Beta1.),as.numeric(co1.beta.pla2$Beta1.))
wilcox.test(as.numeric(cv1.beta.mi7.bk$Beta1.),as.numeric(cv1.beta.pla2$Beta1.))
wilcox.test(as.numeric(co1.beta.mi7.bk$Beta1.),as.numeric(co1.beta.pla2$Beta1.))
wilcox.test(as.numeric(cv1.beta.nv.pla2$Beta1.),as.numeric(cv1.beta.pla2$Beta1.))
wilcox.test(as.numeric(co1.beta.nv.pla2$Beta1.),as.numeric(co1.beta.pla2$Beta1.))

### 4.4 Plot----------------------------------------------------------------

par(mfrow=c(3,2))
boxplot(as.numeric(cv1.beta.mi1.bk$Beta1.),as.numeric(cv1.beta.nv.svmp$Beta1.),as.numeric(cv1.beta.svmp$Beta1.),border=c('black','black','lightblue3'),names=c('Chr. 9','NV','SVMP'),outline=F,ylab="ß")
boxplot(as.numeric(co1.beta.mi1.bk$Beta1.),as.numeric(co1.beta.nv.svmp$Beta1.),as.numeric(co1.beta.svmp$Beta1.),border=c('black','black','lightblue3'),names=c('Chr. 9','NV','SVMP'),outline=F,ylab="ß")
boxplot(as.numeric(cv1.beta.mi2.bk$Beta1.),as.numeric(cv1.beta.nv.svsp$Beta1.),as.numeric(cv1.beta.svsp$Beta1.),border=c('black','black','aquamarine3'),names=c('Chr. 10','NV','SVSP'),outline=F,ylab="ß")
boxplot(as.numeric(co1.beta.mi2.bk$Beta1.),as.numeric(co1.beta.nv.svsp$Beta1.),as.numeric(co1.beta.svsp$Beta1.),border=c('black','black','aquamarine3'),names=c('Chr. 10','NV','SVSP'),outline=F,ylab="ß")
boxplot(as.numeric(cv1.beta.mi7.bk$Beta1.),as.numeric(cv1.beta.nv.pla2$Beta1.),as.numeric(cv1.beta.pla2$Beta1.),border=c('black','black','maroon'),names=c('Chr. 15','NV','PLA2'),outline=F,ylab="ß")
boxplot(as.numeric(co1.beta.mi7.bk$Beta1.),as.numeric(co1.beta.nv.pla2$Beta1.),as.numeric(co1.beta.pla2$Beta1.),border=c('black','black','maroon'),names=c('Chr. 15','NV','PLA2'),outline=F,ylab="ß")

### 5 Boxplots comparing gene means to chromosome backgrounds---------------

# Make vectors of gene means
cv1.taj.svmp.m <- c(0.760,0.260,0.793,0.456,-0.002,0.329,0.432,-0.549,0.062,-0.072,-0.644)
cv2.taj.svmp.m <- c(0.784,0.381,1.288,0.859,0.194,0.358,1.110,0.515,0.870,0.682,0.583)
co1.taj.svmp.m <- c(0.623,1.094,2.074,0.964)
co2.taj.svmp.m <- c(0.730,1.191,2.256,1.522)

cv1.taj.svsp.m <- c(-0.744,-0.024,-0.186,0.380,1.038,0.626,0.181,0.112,0.219,0.413,0.686)
cv2.taj.svsp.m <- c(0.515,0.529,0.235,-0.169,1.550,1.691,0.916,-0.203,0.228,0.428,0.712)
co1.taj.svsp.m <- c(0.021,0.106,0.386,0.167,-0.021,0.815,0.722)
co2.taj.svsp.m <- c(0.798,0.122,0.810,1.352,0.153,1.064,1.381)

cv1.taj.pla2.m <- c(-0.650,-0.997,0.392,0.127)
cv2.taj.pla2.m <- c(0.220,0.751,2.136,2.505)
co1.taj.pla2.m <- c(-0.643)
co2.taj.pla2.m <- c(0.231)

cv1co1.df.svmp.m <- c(0.001,0.029,0.035,0.006)
cv1co1.df.svsp.m <- c(0.000,0.000,0.001,0.000,0.001,0.010,0.001)
cv1co1.df.pla2.m <- c(0.000)

cv1.ihs.svmp.m <- c(0.778,1.042,0.948,1.041,1.241,1.390,1.292,1.573,1.374,0.978,0.676)
co1.ihs.svmp.m <- c(0.854,0.793,0.610,-0.190)

cv1.ihs.svsp.m <- c(1.366,1.087,0.710,0.741,0.420,0.932,0.743,0.771,0.933,0.721,0.969)
co1.ihs.svsp.m <- c(1.239,1.013,0.725,0.905,0.814,0.573,0.905)

cv1.ihs.pla2.m <- c(1.184,0.867,1.059,1.153)
co1.ihs.pla2.m <- c(1.035)

cv1.beta.svmp.m <- c(2.376,1.580,1.325,1.067,1.909,1.195,1.326,1.748,1.259,1.248,0.541)
co1.beta.svmp.m <- c(0.819,1.707,1.680,0.597)

cv1.beta.svsp.m <- c(0.850,1.464,0.177,0.988,3.272,1.794,0.759,0.305,2.401,1.430,1.239)
co1.beta.svsp.m <- c(1.791,0.706,1.322,0.410,0.355,0.544,0.777)

cv1.beta.pla2.m <- c(1.401,0.848,-0.060,2.958)
co1.beta.pla2.m <- c(1.900)

# Plot

par(mfrow=c(4,3))
boxplot(cv1.taj.mi1.bk$TajimaD,cv2.taj.mi1.bk$TajimaD,co1.taj.mi1.bk$TajimaD,co2.taj.mi1.bk$TajimaD,col='lightblue3',outline=F,ylab="Tajima's D",names=c('CV1','CV2','CO1','CO2'))
for (p in cv1.taj.svmp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv2.taj.svmp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.taj.svmp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co2.taj.svmp.m){
  points(4,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1.taj.mi2.bk$TajimaD,cv2.taj.mi2.bk$TajimaD,co1.taj.mi2.bk$TajimaD,co2.taj.mi2.bk$TajimaD,col='aquamarine3',outline=F,ylab="Tajima's D",names=c('CV1','CV2','CO1','CO2'))
for (p in cv1.taj.svsp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv2.taj.svsp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.taj.svsp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co2.taj.svsp.m){
  points(4,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1.taj.mi7.bk$TajimaD,cv2.taj.mi7.bk$TajimaD,co1.taj.mi7.bk$TajimaD,co2.taj.mi7.bk$TajimaD,col='maroon',outline=F,ylab="Tajima's D",names=c('CV1','CV2','CO1','CO2'))
for (p in cv1.taj.pla2.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv2.taj.pla2.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.taj.pla2.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co2.taj.pla2.m){
  points(4,p,pch=24,lwd=1.5,col='navy')
}

boxplot(cv1co1.df.mi1.bk$df,col='lightblue3',outline=F,ylab="df",names=c('CV1 vs CO1'))
for (p in cv1co1.df.svmp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.df.mi2.bk$df,col='aquamarine3',outline=F,ylab="df",names=c('CV1 vs CO1'))
for (p in cv1co1.df.svsp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.df.mi7.bk$df,col='maroon',outline=F,ylab="df",names=c('CV1 vs CO1'))
for (p in cv1co1.df.pla2.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(as.numeric(cv1.ihs.mi1.bk$iHS),as.numeric(co1.ihs.mi1.bk$iHS),col='lightblue3',outline=F,ylab="iHS",names=c('CV1','CO1'),ylim=c(-1,2))
for (p in cv1.ihs.svmp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.ihs.svmp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(as.numeric(cv1.ihs.mi2.bk$iHS),as.numeric(co1.ihs.mi2.bk$iHS),col='aquamarine3',outline=F,ylab="iHS",names=c('CV1','CO1'),ylim=c(-1,2))
for (p in cv1.ihs.svsp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.ihs.svsp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(as.numeric(cv1.ihs.mi7.bk$iHS),as.numeric(co1.ihs.mi7.bk$iHS),col='maroon',outline=F,ylab="iHS",names=c('CV1','CO1'),ylim=c(-1,2))
for (p in cv1.ihs.pla2.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.ihs.pla2.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(as.numeric(cv1.beta.mi1.bk$Beta1.),as.numeric(co1.beta.mi1.bk$Beta1.),col='lightblue3',outline=F,ylab="ß",names=c('CV1','CO1'),ylim=c(-1,3.5))
for (p in cv1.beta.svmp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.beta.svmp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(as.numeric(cv1.beta.mi2.bk$Beta1.),as.numeric(co1.beta.mi2.bk$Beta1.),col='aquamarine3',outline=F,ylab="ß",names=c('CV1','CO1'),ylim=c(-1,3.5))
for (p in cv1.beta.svsp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.beta.svsp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(as.numeric(cv1.beta.mi7.bk$Beta1.),as.numeric(co1.beta.mi7.bk$Beta1.),col='maroon',outline=F,ylab="ß",names=c('CV1','CO1'),ylim=c(-1,3.5))
for (p in cv1.beta.pla2.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.beta.pla2.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
