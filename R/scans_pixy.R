############################################################################
# Genomic scans of genetic diversity and differentitation (pixy)
############################################################################

### Goal: examine genomic scans of pi, dxy, and fst between rattlesnake
### species/populations across venom gene regions.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('./pixy/')

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

### Calculate genome-wide means---------------------------------------------

pi.100kb <- read.table('pixy.all.100kb_pi.txt',header=T)
dxy.100kb <- read.table('pixy.all.100kb_dxy.txt',header=T)
fst.100kb <- read.table('pixy.all.100kb_fst.txt',header=T)

pi.10kb <- read.table('pixy.all.10kb_pi.txt',header=T)
dxy.10kb <- read.table('pixy.all.10kb_dxy.txt',header=T)
fst.10kb <- read.table('pixy.all.10kb_fst.txt',header=T)

pi.100kb.cv1 <- pi.100kb %>% filter(str_detect(pop, "CV1"))
pi.100kb.cv2 <- pi.100kb %>% filter(str_detect(pop, "CV2"))
pi.100kb.co1 <- pi.100kb %>% filter(str_detect(pop, "CO1"))
pi.100kb.co2 <- pi.100kb %>% filter(str_detect(pop, "CO2"))

pi.10kb.cv1 <- pi.10kb %>% filter(str_detect(pop,'CV1'))
pi.10kb.cv2 <- pi.10kb %>% filter(str_detect(pop,'CV2'))
pi.10kb.co1 <- pi.10kb %>% filter(str_detect(pop,'CO1'))
pi.10kb.co2 <- pi.10kb %>% filter(str_detect(pop,'CO2'))

dxy.100kb.cv1co1 <- dxy.100kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CO1'))
dxy.100kb.cv1cv2 <- dxy.100kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CV2'))
dxy.100kb.co1co2 <- dxy.100kb %>% filter(str_detect(pop1, 'CO1') & str_detect(pop2, 'CO2'))

dxy.10kb.cv1co1 <- dxy.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CO1'))
dxy.10kb.cv1cv2 <- dxy.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CV2'))
dxy.10kb.co1co2 <- dxy.10kb %>% filter(str_detect(pop1, 'CO1') & str_detect(pop2, 'CO2'))

fst.100kb.cv1co1 <- fst.100kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CO1'))
fst.100kb.cv1cv2 <- fst.100kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CV2'))
fst.100kb.co1co2 <- fst.100kb %>% filter(str_detect(pop1, 'CO1') & str_detect(pop2, 'CO2'))

fst.10kb.cv1co1 <- fst.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CO1'))
fst.10kb.cv1cv2 <- fst.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CV2'))
fst.10kb.co1co2 <- fst.10kb %>% filter(str_detect(pop1, 'CO1') & str_detect(pop2, 'CO2'))

cv1pi.avg <- mean(pi.10kb.cv1$avg_pi)
cv2pi.avg <- mean(pi.10kb.cv2$avg_pi)
co1pi.avg <- mean(pi.10kb.co1$avg_pi)
co2pi.avg <- mean(pi.10kb.co2$avg_pi)

cv1cv2.dxy.avg <- mean(dxy.10kb.cv1cv2$avg_dxy,na.rm=T)
co1co2.dxy.avg <- mean(dxy.10kb.co1co2$avg_dxy,na.rm=T)
cv1co1.dxy.avg <- mean(dxy.10kb.cv1co1$avg_dxy,na.rm=T)

cv1cv2.fst.avg <- mean(fst.10kb.cv1cv2$avg_wc_fst,na.rm=T)
co1co2.fst.avg <- mean(fst.10kb.co1co2$avg_wc_fst,na.rm=T)
cv1co1.fst.avg <- mean(fst.10kb.cv1co1$avg_wc_fst,na.rm=T)

###-------------------------------------------------------------------------
### SVMP; microchromosome 1
###-------------------------------------------------------------------------

### Read in data------------------------------------------------------------

# 100 kb windows
pi.mi1.100kb <- read.table('pixy.scaffold-mi1.100kb_pi.txt',header=T)
dxy.mi1.100kb <- read.table('pixy.scaffold-mi1.100kb_dxy.txt',header=T)
fst.mi1.100kb <- read.table('pixy.scaffold-mi1.100kb_fst.txt',header=T)

# 10 kb windows
pi.mi1.10kb <- read.table('pixy.scaffold-mi1.10kb_pi.txt',header=T)
dxy.mi1.10kb <- read.table('pixy.scaffold-mi1.10kb_dxy.txt',header=T)
fst.mi1.10kb <- read.table('pixy.scaffold-mi1.10kb_fst.txt',header=T)

# 1 kb windows
pi.mi1.1kb <- read.table('pixy.scaffold-mi1.1kb_pi.txt',header=T)
dxy.mi1.1kb <- read.table('pixy.scaffold-mi1.1kb_dxy.txt',header=T)
fst.mi1.1kb <- read.table('pixy.scaffold-mi1.1kb_fst.txt',header=T)

### Parse data--------------------------------------------------------------

# 100 kb windows
# Pi
cv1pi.mi1.100kb <- pi.mi1.100kb[which(pi.mi1.100kb$pop=='CV1'),]
cv2pi.mi1.100kb <- pi.mi1.100kb[which(pi.mi1.100kb$pop=='CV2'),]
co1pi.mi1.100kb <- pi.mi1.100kb[which(pi.mi1.100kb$pop=='CO1'),]
co1pi.mi1.100kb <- pi.mi1.100kb[which(pi.mi1.100kb$pop=='CO1'),]
co2pi.mi1.100kb <- pi.mi1.100kb[which(pi.mi1.100kb$pop=='CO2'),]
capi.mi1.100kb <- pi.mi1.100kb[which(pi.mi1.100kb$pop=='CA'),]
crpi.mi1.100kb <- pi.mi1.100kb[which(pi.mi1.100kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CV1' & dxy.mi1.100kb$pop2=='CV2'),]
cv1co1.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CV1' & dxy.mi1.100kb$pop2=='CO1'),]
co1co2.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CO1' & dxy.mi1.100kb$pop2=='CO2'),]
cv2co2.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CV2' & dxy.mi1.100kb$pop2=='CO2'),]
cv2co1.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CV2' & dxy.mi1.100kb$pop2=='CO1'),]
cv1ca.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CA' & dxy.mi1.100kb$pop2=='CV1'),]
cv2ca.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CA' & dxy.mi1.100kb$pop2=='CV2'),]
co1ca.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CA' & dxy.mi1.100kb$pop2=='CO1'),]
co2ca.dxy.mi1.100kb <- dxy.mi1.100kb[which(dxy.mi1.100kb$pop1=='CA' & dxy.mi1.100kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi1.100kb <- fst.mi1.100kb[which(fst.mi1.100kb$pop1=='CV1' & fst.mi1.100kb$pop2=='CV2'),]
cv1co1.fst.mi1.100kb <- fst.mi1.100kb[which(fst.mi1.100kb$pop1=='CV1' & fst.mi1.100kb$pop2=='CO1'),]
co1co2.fst.mi1.100kb <- fst.mi1.100kb[which(fst.mi1.100kb$pop1=='CO1' & fst.mi1.100kb$pop2=='CO2'),]
cv2co2.fst.mi1.100kb <- fst.mi1.100kb[which(fst.mi1.100kb$pop1=='CV2' & fst.mi1.100kb$pop2=='CO2'),]

# 10 kb windows
# Pi
cv1pi.mi1.10kb <- pi.mi1.10kb[which(pi.mi1.10kb$pop=='CV1'),]
cv2pi.mi1.10kb <- pi.mi1.10kb[which(pi.mi1.10kb$pop=='CV2'),]
co1pi.mi1.10kb <- read.table('pi.scaffold-mi1.10kb.co1.cnvMask.txt',header=T)
co2pi.mi1.10kb <- read.table('pi.scaffold-mi1.10kb.co2.cnvMask.txt',header=T)
capi.mi1.10kb <- pi.mi1.10kb[which(pi.mi1.10kb$pop=='CA'),]
crpi.mi1.10kb <- pi.mi1.10kb[which(pi.mi1.10kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi1.10kb <- dxy.mi1.10kb[which(dxy.mi1.10kb$pop1=='CV1' & dxy.mi1.10kb$pop2=='CV2'),]
cv1co1.dxy.mi1.10kb <- read.table('dxy.scaffold-mi1.10kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi1.10kb <- read.table('dxy.scaffold-mi1.10kb.co1co2.cnvMask.txt',header=T)
cv2co2.dxy.mi1.10kb <- dxy.mi1.10kb[which(dxy.mi1.10kb$pop1=='CV2' & dxy.mi1.10kb$pop2=='CO2'),]
cv2co1.dxy.mi1.10kb <- dxy.mi1.10kb[which(dxy.mi1.10kb$pop1=='CV2' & dxy.mi1.10kb$pop2=='CO1'),]
cv1ca.dxy.mi1.10kb <- dxy.mi1.10kb[which(dxy.mi1.10kb$pop1=='CA' & dxy.mi1.10kb$pop2=='CV1'),]
cv2ca.dxy.mi1.10kb <- dxy.mi1.10kb[which(dxy.mi1.10kb$pop1=='CA' & dxy.mi1.10kb$pop2=='CV2'),]
co1ca.dxy.mi1.10kb <- dxy.mi1.10kb[which(dxy.mi1.10kb$pop1=='CA' & dxy.mi1.10kb$pop2=='CO1'),]
co2ca.dxy.mi1.10kb <- dxy.mi1.10kb[which(dxy.mi1.10kb$pop1=='CA' & dxy.mi1.10kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi1.10kb <- fst.mi1.10kb[which(fst.mi1.10kb$pop1=='CV1' & fst.mi1.10kb$pop2=='CV2'),]
cv1co1.fst.mi1.10kb <- read.table('fst.scaffold-mi1.10kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi1.10kb <- read.table('fst.scaffold-mi1.10kb.co1co2.cnvMask.txt',header=T)

# 1 kb windows
# Pi
cv1pi.mi1.1kb <- pi.mi1.1kb[which(pi.mi1.1kb$pop=='CV1'),]
cv2pi.mi1.1kb <- pi.mi1.1kb[which(pi.mi1.1kb$pop=='CV2'),]
co1pi.mi1.1kb <- read.table('pi.scaffold-mi1.1kb.co1.cnvMask.txt',header=T)
co2pi.mi1.1kb <- read.table('pi.scaffold-mi1.1kb.co2.cnvMask.txt',header=T)
capi.mi1.1kb <- pi.mi1.1kb[which(pi.mi1.1kb$pop=='CA'),]
crpi.mi1.1kb <- pi.mi1.1kb[which(pi.mi1.1kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi1.1kb <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CV1' & dxy.mi1.1kb$pop2=='CV2'),]
cv1co1.dxy.mi1.1kb <- read.table('dxy.scaffold-mi1.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi1.1kb <- read.table('dxy.scaffold-mi1.1kb.co1co2.cnvMask.txt',header=T)
cv2co2.dxy.mi1.1kb <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CV2' & dxy.mi1.1kb$pop2=='CO2'),]
cv2co1.dxy.mi1.1kb <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CV2' & dxy.mi1.1kb$pop2=='CO1'),]
cv1ca.dxy.mi1.1kb <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CA' & dxy.mi1.1kb$pop2=='CV1'),]
cv2ca.dxy.mi1.1kb <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CA' & dxy.mi1.1kb$pop2=='CV2'),]
co1ca.dxy.mi1.1kb <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CA' & dxy.mi1.1kb$pop2=='CO1'),]
co2ca.dxy.mi1.1kb <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CA' & dxy.mi1.1kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi1.1kb <- fst.mi1.1kb[which(fst.mi1.1kb$pop1=='CV1' & fst.mi1.1kb$pop2=='CV2'),]
cv1co1.fst.mi1.1kb <- read.table('fst.scaffold-mi1.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi1.1kb <- read.table('fst.scaffold-mi1.1kb.co1co2.cnvMask.txt',header=T)

### Plot scans--------------------------------------------------------------

# Export 4 row, 1 column scans at 9 x 3.5
par(mfrow=c(6,1))
## pi
# 10 kb points, 100 kb windows
plot(cv1pi.mi1.10kb$window_pos_1,cv1pi.mi1.10kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(cv1pi.mi1.100kb$window_pos_1,cv1pi.mi1.100kb$avg_pi,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
cv1pi.mi1.avg <- mean(cv1pi.mi1.10kb$avg_pi)
abline(h=cv1pi.mi1.avg,lty=3,col='red')
abline(h=cv1pi.avg,lty=3,col='blue')

plot(cv2pi.mi1.10kb$window_pos_1,cv2pi.mi1.10kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(cv2pi.mi1.100kb$window_pos_1,cv2pi.mi1.100kb$avg_pi,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
cv2pi.mi1.avg <- mean(cv2pi.mi1.10kb$avg_pi)
abline(h=cv2pi.mi1.avg,lty=3,col='red')
abline(h=cv2pi.avg,lty=3,col='blue')

plot(co1pi.mi1.10kb$window_pos_1,co1pi.mi1.10kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(co1pi.mi1.100kb$window_pos_1,co1pi.mi1.100kb$avg_pi,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
co1pi.mi1.avg <- mean(co1pi.mi1.10kb$avg_pi,na.rm=T)
abline(h=co1pi.mi1.avg,lty=3,col='red')
abline(h=co1pi.avg,lty=3,col='blue')

plot(co2pi.mi1.10kb$window_pos_1,co2pi.mi1.10kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(co2pi.mi1.100kb$window_pos_1,co2pi.mi1.100kb$avg_pi,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
co2pi.mi1.avg <- mean(co2pi.mi1.10kb$avg_pi,na.rm=T)
abline(h=co2pi.mi1.avg,lty=3,col='red')
abline(h=co2pi.avg,lty=3,col='blue')

# Zoomed, 1kb points, 10 kb windows
par(mfrow=c(6,1))
plot(cv1pi.mi1.1kb$window_pos_1,cv1pi.mi1.1kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.002,0.015),ylab='pi',xlab='Chromosome Position')
lines(cv1pi.mi1.10kb$window_pos_1,cv1pi.mi1.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv1pi.mi1.avg,lty=3,col='red')
abline(h=cv1pi.avg,lty=3,col='blue')

plot(cv2pi.mi1.1kb$window_pos_1,cv2pi.mi1.1kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.002,0.015),ylab='pi',xlab='Chromosome Position')
lines(cv2pi.mi1.10kb$window_pos_1,cv2pi.mi1.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv2pi.mi1.avg,lty=3,col='red')
abline(h=cv2pi.avg,lty=3,col='blue')

plot(co1pi.mi1.1kb$window_pos_1,co1pi.mi1.1kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.002,0.03),ylab='pi',xlab='Chromosome Position')
lines(co1pi.mi1.10kb$window_pos_1,co1pi.mi1.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=co1pi.mi1.avg,lty=3,col='red')
abline(h=co1pi.avg,lty=3,col='blue')

plot(co2pi.mi1.1kb$window_pos_1,co2pi.mi1.1kb$avg_pi,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.002,0.03),ylab='pi',xlab='Chromosome Position')
lines(co2pi.mi1.10kb$window_pos_1,co2pi.mi1.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=co2pi.mi1.avg,lty=3,col='red')
abline(h=co2pi.avg,lty=3,col='blue')

## dxy
par(mfrow=c(3,1))
# 10 kb points, 100 kb windows
plot(cv1cv2.dxy.mi1.10kb$window_pos_1,cv1cv2.dxy.mi1.10kb$avg_dxy,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(cv1cv2.dxy.mi1.100kb$window_pos_1,cv1cv2.dxy.mi1.100kb$avg_dxy,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
cv1cv2.dxy.mi1.avg <- mean(cv1cv2.dxy.mi1.10kb$avg_dxy,na.rm=T)
abline(h=cv1cv2.dxy.mi1.avg,lty=3,col='red')
abline(h=cv1cv2.dxy.avg,lty=3,col='blue')

plot(co1co2.dxy.mi1.10kb$window_pos_1,co1co2.dxy.mi1.10kb$avg_dxy,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(co1co2.dxy.mi1.100kb$window_pos_1,co1co2.dxy.mi1.100kb$avg_dxy,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
co1co2.dxy.mi1.avg <- mean(co1co2.dxy.mi1.10kb$avg_dxy,na.rm=T)
abline(h=co1co2.dxy.mi1.avg,lty=3,col='red')
abline(h=co1co2.dxy.avg,lty=3,col='blue')

plot(cv1co1.dxy.mi1.10kb$window_pos_1,cv1co1.dxy.mi1.10kb$avg_dxy,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(cv1co1.dxy.mi1.100kb$window_pos_1,cv1co1.dxy.mi1.100kb$avg_dxy,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
cv1co1.dxy.mi1.avg <- mean(cv1co1.dxy.mi1.10kb$avg_dxy,na.rm=T)
abline(h=cv1co1.dxy.mi1.avg,lty=3,col='red')
abline(h=cv1co1.dxy.avg,lty=3,col='blue')

# Zoomed, 1kb points, 10 kb windows
plot(cv1cv2.dxy.mi1.1kb$window_pos_1,cv1cv2.dxy.mi1.1kb$avg_dxy,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.002,0.04),ylab='dxy',xlab='Chromosome Position')
lines(cv1cv2.dxy.mi1.10kb$window_pos_1,cv1cv2.dxy.mi1.10kb$avg_dxy,type='l')
vert <- -0.0015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv1cv2.dxy.mi1.avg,lty=3,col='red')
abline(h=cv1cv2.dxy.avg,lty=3,col='blue')

plot(co1co2.dxy.mi1.1kb$window_pos_1,co1co2.dxy.mi1.1kb$avg_dxy,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.002,0.04),ylab='dxy',xlab='Chromosome Position')
lines(co1co2.dxy.mi1.10kb$window_pos_1,co1co2.dxy.mi1.10kb$avg_dxy,type='l')
vert <- -0.0015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=co1co2.dxy.mi1.avg,lty=3,col='red')
abline(h=co1co2.dxy.avg,lty=3,col='blue')

plot(cv1co1.dxy.mi1.1kb$window_pos_1,cv1co1.dxy.mi1.1kb$avg_dxy,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.002,0.04),ylab='dxy',xlab='Chromosome Position')
lines(cv1co1.dxy.mi1.10kb$window_pos_1,cv1co1.dxy.mi1.10kb$avg_dxy,type='l')
vert <- -0.0015
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.dxy.mi1.avg,lty=3,col='red')
abline(h=cv1co1.dxy.avg,lty=3,col='blue')

## fst
par(mfrow=c(3,1))
# 10 kb points, 100 kb windows
plot(cv1cv2.fst.mi1.10kb$window_pos_1,cv1cv2.fst.mi1.10kb$avg_wc_fst,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(cv1cv2.fst.mi1.100kb$window_pos_1,cv1cv2.fst.mi1.100kb$avg_wc_fst,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
cv1cv2.fst.mi1.avg <- mean(cv1cv2.fst.mi1.10kb$avg_wc_fst,na.rm=T)
abline(h=cv1cv2.fst.mi1.avg,lty=3,col='red')
abline(h=cv1cv2.fst.avg,lty=3,col='blue')
plot(co1co2.fst.mi1.10kb$window_pos_1,co1co2.fst.mi1.10kb$avg_wc_fst,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(co1co2.fst.mi1.100kb$window_pos_1,co1co2.fst.mi1.100kb$avg_wc_fst,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
co1co2.fst.mi1.avg <- mean(co1co2.fst.mi1.10kb$avg_wc_fst,na.rm=T)
abline(h=co1co2.fst.mi1.avg,lty=3,col='red')
abline(h=co1co2.fst.avg,lty=3,col='blue')

plot(cv1co1.fst.mi1.10kb$window_pos_1,cv1co1.fst.mi1.10kb$avg_wc_fst,pch=20,col=alpha('lightblue3',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(cv1co1.fst.mi1.100kb$window_pos_1,cv1co1.fst.mi1.100kb$avg_wc_fst,type='l')
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
cv1co1.fst.mi1.avg <- mean(cv1co1.fst.mi1.10kb$avg_wc_fst,na.rm=T)
abline(h=cv1co1.fst.mi1.avg,lty=3,col='red')
abline(h=cv1co1.fst.avg,lty=3,col='blue')

# Zoomed, 1kb points, 10 kb windows
plot(cv1cv2.fst.mi1.1kb$window_pos_1,cv1cv2.fst.mi1.1kb$avg_wc_fst,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(cv1cv2.fst.mi1.10kb$window_pos_1,cv1cv2.fst.mi1.10kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv1cv2.fst.mi1.avg,lty=3,col='red')
abline(h=cv1cv2.fst.avg,lty=3,col='blue')

plot(co1co2.fst.mi1.1kb$window_pos_1,co1co2.fst.mi1.1kb$avg_wc_fst,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(co1co2.fst.mi1.10kb$window_pos_1,co1co2.fst.mi1.10kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=co1co2.fst.mi1.avg,lty=3,col='red')
abline(h=co1co2.fst.avg,lty=3,col='blue')

plot(cv1co1.fst.mi1.1kb$window_pos_1,cv1co1.fst.mi1.1kb$avg_wc_fst,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(cv1co1.fst.mi1.10kb$window_pos_1,cv1co1.fst.mi1.10kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.fst.mi1.avg,lty=3,col='red')
abline(h=cv1co1.fst.avg,lty=3,col='blue')


###-------------------------------------------------------------------------
### SVSP; microchromosome 2
###-------------------------------------------------------------------------

### Read in data------------------------------------------------------------

# 100 kb windows
pi.mi2.100kb <- read.table('pixy.scaffold-mi2.100kb_pi.txt',header=T)
dxy.mi2.100kb <- read.table('pixy.scaffold-mi2.100kb_dxy.txt',header=T)
fst.mi2.100kb <- read.table('pixy.scaffold-mi2.100kb_fst.txt',header=T)

# 10 kb windows
pi.mi2.10kb <- read.table('pixy.scaffold-mi2.10kb_pi.txt',header=T)
dxy.mi2.10kb <- read.table('pixy.scaffold-mi2.10kb_dxy.txt',header=T)
fst.mi2.10kb <- read.table('pixy.scaffold-mi2.10kb_fst.txt',header=T)

# 1 kb windows
pi.mi2.1kb <- read.table('pixy.scaffold-mi2.1kb_pi.txt',header=T)
dxy.mi2.1kb <- read.table('pixy.scaffold-mi2.1kb_dxy.txt',header=T)
fst.mi2.1kb <- read.table('pixy.scaffold-mi2.1kb_fst.txt',header=T)

### Parse data--------------------------------------------------------------

# 100 kb windows
# Pi
cv1pi.mi2.100kb <- pi.mi2.100kb[which(pi.mi2.100kb$pop=='CV1'),]
cv2pi.mi2.100kb <- pi.mi2.100kb[which(pi.mi2.100kb$pop=='CV2'),]
co1pi.mi2.100kb <- pi.mi2.100kb[which(pi.mi2.100kb$pop=='CO1'),]
co2pi.mi2.100kb <- pi.mi2.100kb[which(pi.mi2.100kb$pop=='CO2'),]
capi.mi2.100kb <- pi.mi2.100kb[which(pi.mi2.100kb$pop=='CA'),]
crpi.mi2.100kb <- pi.mi2.100kb[which(pi.mi2.100kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CV1' & dxy.mi2.100kb$pop2=='CV2'),]
cv1co1.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CV1' & dxy.mi2.100kb$pop2=='CO1'),]
co1co2.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CO1' & dxy.mi2.100kb$pop2=='CO2'),]
cv2co2.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CV2' & dxy.mi2.100kb$pop2=='CO2'),]
cv2co1.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CV2' & dxy.mi2.100kb$pop2=='CO1'),]
cv1ca.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CA' & dxy.mi2.100kb$pop2=='CV1'),]
cv2ca.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CA' & dxy.mi2.100kb$pop2=='CV2'),]
co1ca.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CA' & dxy.mi2.100kb$pop2=='CO1'),]
co2ca.dxy.mi2.100kb <- dxy.mi2.100kb[which(dxy.mi2.100kb$pop1=='CA' & dxy.mi2.100kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi2.100kb <- fst.mi2.100kb[which(fst.mi2.100kb$pop1=='CV1' & fst.mi2.100kb$pop2=='CV2'),]
cv1co1.fst.mi2.100kb <- fst.mi2.100kb[which(fst.mi2.100kb$pop1=='CV1' & fst.mi2.100kb$pop2=='CO1'),]
co1co2.fst.mi2.100kb <- fst.mi2.100kb[which(fst.mi2.100kb$pop1=='CO1' & fst.mi2.100kb$pop2=='CO2'),]
cv2co2.fst.mi2.100kb <- fst.mi2.100kb[which(fst.mi2.100kb$pop1=='CV2' & fst.mi2.100kb$pop2=='CO2'),]

# 10 kb windows
# Pi
cv1pi.mi2.10kb <- pi.mi2.10kb[which(pi.mi2.10kb$pop=='CV1'),]
cv2pi.mi2.10kb <- pi.mi2.10kb[which(pi.mi2.10kb$pop=='CV2'),]
co1pi.mi2.10kb <- read.table('pi.scaffold-mi2.10kb.co1.cnvMask.txt',header=T)
co2pi.mi2.10kb <- read.table('pi.scaffold-mi2.10kb.co2.cnvMask.txt',header=T)
capi.mi2.10kb <- pi.mi2.10kb[which(pi.mi2.10kb$pop=='CA'),]
crpi.mi2.10kb <- pi.mi2.10kb[which(pi.mi2.10kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi2.10kb <- dxy.mi2.10kb[which(dxy.mi2.10kb$pop1=='CV1' & dxy.mi2.10kb$pop2=='CV2'),]
cv1co1.dxy.mi2.10kb <- read.table('dxy.scaffold-mi2.10kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi2.10kb <- read.table('dxy.scaffold-mi2.10kb.co1co2.cnvMask.txt',header=T)
cv2co2.dxy.mi2.10kb <- dxy.mi2.10kb[which(dxy.mi2.10kb$pop1=='CV2' & dxy.mi2.10kb$pop2=='CO2'),]
cv2co1.dxy.mi2.10kb <- dxy.mi2.10kb[which(dxy.mi2.10kb$pop1=='CV2' & dxy.mi2.10kb$pop2=='CO1'),]
cv1ca.dxy.mi2.10kb <- dxy.mi2.10kb[which(dxy.mi2.10kb$pop1=='CA' & dxy.mi2.10kb$pop2=='CV1'),]
cv2ca.dxy.mi2.10kb <- dxy.mi2.10kb[which(dxy.mi2.10kb$pop1=='CA' & dxy.mi2.10kb$pop2=='CV2'),]
co1ca.dxy.mi2.10kb <- dxy.mi2.10kb[which(dxy.mi2.10kb$pop1=='CA' & dxy.mi2.10kb$pop2=='CO1'),]
co2ca.dxy.mi2.10kb <- dxy.mi2.10kb[which(dxy.mi2.10kb$pop1=='CA' & dxy.mi2.10kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi2.10kb <- fst.mi2.10kb[which(fst.mi2.10kb$pop1=='CV1' & fst.mi2.10kb$pop2=='CV2'),]
cv1co1.fst.mi2.10kb <- read.table('fst.scaffold-mi2.10kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi2.10kb <- read.table('fst.scaffold-mi2.10kb.co1co2.cnvMask.txt',header=T)
cv2co2.fst.mi2.10kb <- fst.mi2.10kb[which(fst.mi2.10kb$pop1=='CV2' & fst.mi2.10kb$pop2=='CO2'),]

# 1 kb windows
# Pi
cv1pi.mi2.1kb <- pi.mi2.1kb[which(pi.mi2.1kb$pop=='CV1'),]
cv2pi.mi2.1kb <- pi.mi2.1kb[which(pi.mi2.1kb$pop=='CV2'),]
co1pi.mi2.1kb <- read.table('pi.scaffold-mi2.1kb.co1.cnvMask.txt',header=T)
co2pi.mi2.1kb <- read.table('pi.scaffold-mi2.1kb.co2.cnvMask.txt',header=T)
capi.mi2.1kb <- pi.mi2.1kb[which(pi.mi2.1kb$pop=='CA'),]
crpi.mi2.1kb <- pi.mi2.1kb[which(pi.mi2.1kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi2.1kb <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CV1' & dxy.mi2.1kb$pop2=='CV2'),]
cv1co1.dxy.mi2.1kb <- read.table('dxy.scaffold-mi2.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi2.1kb <- read.table('dxy.scaffold-mi2.1kb.co1co2.cnvMask.txt',header=T)
cv2co2.dxy.mi2.1kb <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CV2' & dxy.mi2.1kb$pop2=='CO2'),]
cv2co1.dxy.mi2.1kb <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CV2' & dxy.mi2.1kb$pop2=='CO1'),]
cv1ca.dxy.mi2.1kb <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CA' & dxy.mi2.1kb$pop2=='CV1'),]
cv2ca.dxy.mi2.1kb <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CA' & dxy.mi2.1kb$pop2=='CV2'),]
co1ca.dxy.mi2.1kb <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CA' & dxy.mi2.1kb$pop2=='CO1'),]
co2ca.dxy.mi2.1kb <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CA' & dxy.mi2.1kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi2.1kb <- fst.mi2.1kb[which(fst.mi2.1kb$pop1=='CV1' & fst.mi2.1kb$pop2=='CV2'),]
cv1co1.fst.mi2.1kb <- read.table('fst.scaffold-mi2.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi2.1kb <- read.table('fst.scaffold-mi2.1kb.co1co2.cnvMask.txt',header=T)
cv2co2.fst.mi2.1kb <- fst.mi2.1kb[which(fst.mi2.1kb$pop1=='CV2' & fst.mi2.1kb$pop2=='CO2'),]

### Plot scans--------------------------------------------------------------

## pi
par(mfrow=c(6,1))
# 10 kb points, 100 kb windows
plot(cv1pi.mi2.10kb$window_pos_1,cv1pi.mi2.10kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='pi',ylim=c(0,0.015))
lines(cv1pi.mi2.100kb$window_pos_1,cv1pi.mi2.100kb$avg_pi,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
cv1pi.mi2.avg <- mean(cv1pi.mi2.10kb$avg_pi,na.rm=T)
abline(h=cv1pi.mi2.avg,lty=3,col='red')
abline(h=cv1pi.avg,lty=3,col='blue')

plot(cv2pi.mi2.10kb$window_pos_1,cv2pi.mi2.10kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(cv2pi.mi2.100kb$window_pos_1,cv2pi.mi2.100kb$avg_pi,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
cv2pi.mi2.avg <- mean(cv2pi.mi2.10kb$avg_pi,na.rm=T)
abline(h=cv2pi.mi2.avg,lty=3,col='red')
abline(h=cv2pi.avg,lty=3,col='blue')

plot(co1pi.mi2.10kb$window_pos_1,co1pi.mi2.10kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='pi',ylim=c(0,0.015))
lines(co1pi.mi2.100kb$window_pos_1,co1pi.mi2.100kb$avg_pi,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
co1pi.mi2.avg <- mean(co1pi.mi2.10kb$avg_pi,na.rm=T)
abline(h=co1pi.mi2.avg,lty=3,col='red')
abline(h=co1pi.avg,lty=3,col='blue')

plot(co2pi.mi2.10kb$window_pos_1,co2pi.mi2.10kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(co2pi.mi2.100kb$window_pos_1,co2pi.mi2.100kb$avg_pi,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
co2pi.mi2.avg <- mean(co2pi.mi2.10kb$avg_pi,na.rm=T)
abline(h=co2pi.mi2.avg,lty=3,col='red')
abline(h=co2pi.avg,lty=3,col='blue')

# Zoomed, 1kb points, 10 kb windows
plot(cv1pi.mi2.1kb$window_pos_1,cv1pi.mi2.1kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.002,0.03),ylab='pi',xlab='Chromosome Position')
lines(cv1pi.mi2.10kb$window_pos_1,cv1pi.mi2.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=cv1pi.mi2.avg,lty=3,col='red')
abline(h=cv1pi.avg,lty=3,col='blue')

plot(cv2pi.mi2.1kb$window_pos_1,cv2pi.mi2.1kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.002,0.03),ylab='pi',xlab='Chromosome Position')
lines(cv2pi.mi2.10kb$window_pos_1,cv2pi.mi2.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=cv2pi.mi2.avg,lty=3,col='red')
abline(h=cv2pi.avg,lty=3,col='blue')

plot(co1pi.mi2.1kb$window_pos_1,co1pi.mi2.1kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.002,0.03),ylab='pi',xlab='Chromosome Position')
lines(co1pi.mi2.10kb$window_pos_1,co1pi.mi2.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=co1pi.mi2.avg,lty=3,col='red')
abline(h=co1pi.avg,lty=3,col='blue')

plot(co2pi.mi2.1kb$window_pos_1,co2pi.mi2.1kb$avg_pi,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.002,0.03),ylab='pi',xlab='Chromosome Position')
lines(co2pi.mi2.10kb$window_pos_1,co2pi.mi2.10kb$avg_pi,type='l')
vert <- -0.0015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=co2pi.mi2.avg,lty=3,col='red')
abline(h=co2pi.avg,lty=3,col='blue')

## dxy
par(mfrow=c(3,1))
# 10 kb points, 100 kb windows
plot(cv1cv2.dxy.mi2.10kb$window_pos_1,cv1cv2.dxy.mi2.10kb$avg_dxy,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(cv1cv2.dxy.mi2.100kb$window_pos_1,cv1cv2.dxy.mi2.100kb$avg_dxy,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
cv1cv2.dxy.mi2.avg <- mean(cv1cv2.dxy.mi2.10kb$avg_dxy,na.rm=T)
abline(h=cv1cv2.dxy.mi2.avg,lty=3,col='red')
abline(h=cv1cv2.dxy.avg,lty=3,col='blue')

plot(co1co2.dxy.mi2.10kb$window_pos_1,co1co2.dxy.mi2.10kb$avg_dxy,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(co1co2.dxy.mi2.100kb$window_pos_1,co1co2.dxy.mi2.100kb$avg_dxy,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
co1co2.dxy.mi2.avg <- mean(co1co2.dxy.mi2.10kb$avg_dxy,na.rm=T)
abline(h=co1co2.dxy.mi2.avg,lty=3,col='red')
abline(h=co1co2.dxy.avg,lty=3,col='blue')

plot(cv1co1.dxy.mi2.10kb$window_pos_1,cv1co1.dxy.mi2.10kb$avg_dxy,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='dxy',ylim=c(0,0.025))
lines(cv1co1.dxy.mi2.100kb$window_pos_1,cv1co1.dxy.mi2.100kb$avg_dxy,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
cv1co1.dxy.mi2.avg <- mean(cv1co1.dxy.mi2.10kb$avg_dxy,na.rm=T)
abline(h=cv1co1.dxy.mi2.avg,lty=3,col='red')
abline(h=cv1co1.dxy.avg,lty=3,col='blue')

# Zoomed, 1kb points, 10 kb windows
plot(cv1cv2.dxy.mi2.1kb$window_pos_1,cv1cv2.dxy.mi2.1kb$avg_dxy,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.002,0.045),ylab='dxy',xlab='Chromosome Position')
lines(cv1cv2.dxy.mi2.10kb$window_pos_1,cv1cv2.dxy.mi2.10kb$avg_dxy,type='l')
vert <- -0.0015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=cv1cv2.dxy.mi2.avg,lty=3,col='red')
abline(h=cv1cv2.dxy.avg,lty=3,col='blue')

plot(co1co2.dxy.mi2.1kb$window_pos_1,co1co2.dxy.mi2.1kb$avg_dxy,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.002,0.045),ylab='dxy',xlab='Chromosome Position')
lines(co1co2.dxy.mi2.10kb$window_pos_1,co1co2.dxy.mi2.10kb$avg_dxy,type='l')
vert <- -0.0015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=co1co2.dxy.mi2.avg,lty=3,col='red')
abline(h=co1co2.dxy.avg,lty=3,col='blue')

plot(cv1co1.dxy.mi2.1kb$window_pos_1,cv1co1.dxy.mi2.1kb$avg_dxy,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.002,0.045),ylab='dxy',xlab='Chromosome Position')
lines(cv1co1.dxy.mi2.10kb$window_pos_1,cv1co1.dxy.mi2.10kb$avg_dxy,type='l')
vert <- -0.0015
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.dxy.mi2.avg,lty=3,col='red')
abline(h=cv1co1.dxy.avg,lty=3,col='blue')

## fst
par(mfrow=c(3,1))
# 10 kb points, 100 kb windows
plot(cv1cv2.fst.mi2.10kb$window_pos_1,cv1cv2.fst.mi2.10kb$avg_wc_fst,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(cv1cv2.fst.mi2.100kb$window_pos_1,cv1cv2.fst.mi2.100kb$avg_wc_fst,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
cv1cv2.fst.mi2.avg <- mean(cv1cv2.fst.mi2.10kb$avg_wc_fst,na.rm=T)
abline(h=cv1cv2.fst.mi2.avg,lty=3,col='red')
abline(h=cv1cv2.fst.avg,lty=3,col='blue')

plot(co1co2.fst.mi2.10kb$window_pos_1,co1co2.fst.mi2.10kb$avg_wc_fst,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(co1co2.fst.mi2.100kb$window_pos_1,co1co2.fst.mi2.100kb$avg_wc_fst,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
co1co2.fst.mi2.avg <- mean(co1co2.fst.mi2.10kb$avg_wc_fst,na.rm=T)
abline(h=co1co2.fst.mi2.avg,lty=3,col='red')
abline(h=co1co2.fst.avg,lty=3,col='blue')

plot(cv1co1.fst.mi2.10kb$window_pos_1,cv1co1.fst.mi2.10kb$avg_wc_fst,pch=20,col=alpha('aquamarine3',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(cv1co1.fst.mi2.100kb$window_pos_1,cv1co1.fst.mi2.100kb$avg_wc_fst,type='l')
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
cv1co1.fst.mi2.avg <- mean(cv1co1.fst.mi2.10kb$avg_wc_fst,na.rm=T)
abline(h=cv1co1.fst.mi2.avg,lty=3,col='red')
abline(h=cv1co1.fst.avg,lty=3,col='blue')

# Zoomed, 1kb points, 10 kb windows
plot(cv1cv2.fst.mi2.1kb$window_pos_1,cv1cv2.fst.mi2.1kb$avg_wc_fst,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(cv1cv2.fst.mi2.10kb$window_pos_1,cv1cv2.fst.mi2.10kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=cv1cv2.fst.mi2.avg,lty=3,col='red')
abline(h=cv1cv2.fst.avg,lty=3,col='blue')

plot(co1co2.fst.mi2.1kb$window_pos_1,co1co2.fst.mi2.1kb$avg_wc_fst,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(co1co2.fst.mi2.10kb$window_pos_1,co1co2.fst.mi2.10kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=co1co2.fst.mi2.avg,lty=3,col='red')
abline(h=co1co2.fst.avg,lty=3,col='blue')

plot(cv1co1.fst.mi2.1kb$window_pos_1,cv1co1.fst.mi2.1kb$avg_wc_fst,pch=20,col=alpha('aquamarine3',0.15),xlim=c(8350000,9200000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(cv1co1.fst.mi2.10kb$window_pos_1,cv1co1.fst.mi2.10kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in svsp) {
  segments(x0=svsp$V2,y0=vert,x1=svsp$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.fst.mi2.avg,lty=3,col='red')
abline(h=cv1co1.fst.avg,lty=3,col='blue')

###-------------------------------------------------------------------------
### PLA2; microchromosome 7
###-------------------------------------------------------------------------

### Read in data------------------------------------------------------------

# 100 kb windows
pi.mi7.100kb <- read.table('pixy.scaffold-mi7.100kb_pi.txt',header=T)
dxy.mi7.100kb <- read.table('pixy.scaffold-mi7.100kb_dxy.txt',header=T)
fst.mi7.100kb <- read.table('pixy.scaffold-mi7.100kb_fst.txt',header=T)

# 10 kb windows
pi.mi7.10kb <- read.table('pixy.scaffold-mi7.10kb_pi.txt',header=T)
dxy.mi7.10kb <- read.table('pixy.scaffold-mi7.10kb_dxy.txt',header=T)
fst.mi7.10kb <- read.table('pixy.scaffold-mi7.10kb_fst.txt',header=T)

# 1 kb windows
pi.mi7.1kb <- read.table('pixy.scaffold-mi7.1kb_pi.txt',header=T)
dxy.mi7.1kb <- read.table('pixy.scaffold-mi7.1kb_dxy.txt',header=T)
fst.mi7.1kb <- read.table('pixy.scaffold-mi7.1kb_fst.txt',header=T)

# 250 bp windows
pi.mi7.250bp <- read.table('pixy.scaffold-mi7.250bp_pi.txt',header=T)
dxy.mi7.250bp <- read.table('pixy.scaffold-mi7.250bp_dxy.txt',header=T)
fst.mi7.250bp <- read.table('pixy.scaffold-mi7.250bp_fst.txt',header=T)

### Parse data--------------------------------------------------------------

# 100 kb windows
# Pi
cv1pi.mi7.100kb <- pi.mi7.100kb[which(pi.mi7.100kb$pop=='CV1'),]
cv2pi.mi7.100kb <- pi.mi7.100kb[which(pi.mi7.100kb$pop=='CV2'),]
co1pi.mi7.100kb <- pi.mi7.100kb[which(pi.mi7.100kb$pop=='CO1'),]
co2pi.mi7.100kb <- pi.mi7.100kb[which(pi.mi7.100kb$pop=='CO2'),]
capi.mi7.100kb <- pi.mi7.100kb[which(pi.mi7.100kb$pop=='CA'),]
crpi.mi7.100kb <- pi.mi7.100kb[which(pi.mi7.100kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CV1' & dxy.mi7.100kb$pop2=='CV2'),]
cv1co1.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CV1' & dxy.mi7.100kb$pop2=='CO1'),]
co1co2.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CO1' & dxy.mi7.100kb$pop2=='CO2'),]
cv2co2.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CV2' & dxy.mi7.100kb$pop2=='CO2'),]
cv2co1.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CV2' & dxy.mi7.100kb$pop2=='CO1'),]
cv1ca.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CA' & dxy.mi7.100kb$pop2=='CV1'),]
cv2ca.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CA' & dxy.mi7.100kb$pop2=='CV2'),]
co1ca.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CA' & dxy.mi7.100kb$pop2=='CO1'),]
co2ca.dxy.mi7.100kb <- dxy.mi7.100kb[which(dxy.mi7.100kb$pop1=='CA' & dxy.mi7.100kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi7.100kb <- fst.mi7.100kb[which(fst.mi7.100kb$pop1=='CV1' & fst.mi7.100kb$pop2=='CV2'),]
cv1co1.fst.mi7.100kb <- fst.mi7.100kb[which(fst.mi7.100kb$pop1=='CV1' & fst.mi7.100kb$pop2=='CO1'),]
co1co2.fst.mi7.100kb <- fst.mi7.100kb[which(fst.mi7.100kb$pop1=='CO1' & fst.mi7.100kb$pop2=='CO2'),]
cv2co2.fst.mi7.100kb <- fst.mi7.100kb[which(fst.mi7.100kb$pop1=='CV2' & fst.mi7.100kb$pop2=='CO2'),]

# 10 kb windows
# Pi
cv1pi.mi7.10kb <- pi.mi7.10kb[which(pi.mi7.10kb$pop=='CV1'),]
cv2pi.mi7.10kb <- pi.mi7.10kb[which(pi.mi7.10kb$pop=='CV2'),]
co1pi.mi7.10kb <- read.table('pi.scaffold-mi7.10kb.co1.cnvMask.txt',header=T)
co2pi.mi7.10kb <- read.table('pi.scaffold-mi7.10kb.co2.cnvMask.txt',header=T)
capi.mi7.10kb <- pi.mi7.10kb[which(pi.mi7.10kb$pop=='CA'),]
crpi.mi7.10kb <- pi.mi7.10kb[which(pi.mi7.10kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi7.10kb <- dxy.mi7.10kb[which(dxy.mi7.10kb$pop1=='CV1' & dxy.mi7.10kb$pop2=='CV2'),]
cv1co1.dxy.mi7.10kb <- read.table('dxy.scaffold-mi7.10kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi7.10kb <- read.table('dxy.scaffold-mi7.10kb.co1co2.cnvMask.txt',header=T)
cv2co2.dxy.mi7.10kb <- dxy.mi7.10kb[which(dxy.mi7.10kb$pop1=='CV2' & dxy.mi7.10kb$pop2=='CO2'),]
cv2co1.dxy.mi7.10kb <- dxy.mi7.10kb[which(dxy.mi7.10kb$pop1=='CV2' & dxy.mi7.10kb$pop2=='CO1'),]
cv1ca.dxy.mi7.10kb <- dxy.mi7.10kb[which(dxy.mi7.10kb$pop1=='CA' & dxy.mi7.10kb$pop2=='CV1'),]
cv2ca.dxy.mi7.10kb <- dxy.mi7.10kb[which(dxy.mi7.10kb$pop1=='CA' & dxy.mi7.10kb$pop2=='CV2'),]
co1ca.dxy.mi7.10kb <- dxy.mi7.10kb[which(dxy.mi7.10kb$pop1=='CA' & dxy.mi7.10kb$pop2=='CO1'),]
co2ca.dxy.mi7.10kb <- dxy.mi7.10kb[which(dxy.mi7.10kb$pop1=='CA' & dxy.mi7.10kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi7.10kb <- fst.mi7.10kb[which(fst.mi7.10kb$pop1=='CV1' & fst.mi7.10kb$pop2=='CV2'),]
cv1co1.fst.mi7.10kb <- read.table('fst.scaffold-mi7.10kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi7.10kb <- read.table('fst.scaffold-mi7.10kb.co1co2.cnvMask.txt',header=T)
cv2co2.fst.mi7.10kb <- fst.mi7.10kb[which(fst.mi7.10kb$pop1=='CV2' & fst.mi7.10kb$pop2=='CO2'),]

# 1 kb windows
cv1pi.mi7.1kb <- pi.mi7.1kb[which(pi.mi7.1kb$pop=='CV1'),]
cv2pi.mi7.1kb <- pi.mi7.1kb[which(pi.mi7.1kb$pop=='CV2'),]
co1pi.mi7.1kb <- read.table('pi.scaffold-mi7.1kb.co1.cnvMask.txt',header=T)
co2pi.mi7.1kb <- read.table('pi.scaffold-mi7.1kb.co2.cnvMask.txt',header=T)
capi.mi7.1kb <- pi.mi7.1kb[which(pi.mi7.1kb$pop=='CA'),]
crpi.mi7.1kb <- pi.mi7.1kb[which(pi.mi7.1kb$pop=='CR'),]
# dxy
cv1cv2.dxy.mi7.1kb <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CV1' & dxy.mi7.1kb$pop2=='CV2'),]
cv1co1.dxy.mi7.1kb <- read.table('dxy.scaffold-mi7.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi7.1kb <- read.table('dxy.scaffold-mi7.1kb.co1co2.cnvMask.txt',header=T)
cv2co2.dxy.mi7.1kb <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CV2' & dxy.mi7.1kb$pop2=='CO2'),]
cv2co1.dxy.mi7.1kb <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CV2' & dxy.mi7.1kb$pop2=='CO1'),]
cv1ca.dxy.mi7.1kb <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CA' & dxy.mi7.1kb$pop2=='CV1'),]
cv2ca.dxy.mi7.1kb <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CA' & dxy.mi7.1kb$pop2=='CV2'),]
co1ca.dxy.mi7.1kb <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CA' & dxy.mi7.1kb$pop2=='CO1'),]
co2ca.dxy.mi7.1kb <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CA' & dxy.mi7.1kb$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi7.1kb <- fst.mi7.1kb[which(fst.mi7.1kb$pop1=='CV1' & fst.mi7.1kb$pop2=='CV2'),]
cv1co1.fst.mi7.1kb <- read.table('fst.scaffold-mi7.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi7.1kb <- read.table('fst.scaffold-mi7.1kb.co1co2.cnvMask.txt',header=T)
cv2co2.fst.mi7.1kb <- fst.mi7.1kb[which(fst.mi7.1kb$pop1=='CV2' & fst.mi7.1kb$pop2=='CO2'),]

# 250 bp windows
# Pi
cv1pi.mi7.250bp <- pi.mi7.250bp[which(pi.mi7.250bp$pop=='CV1'),]
cv2pi.mi7.250bp <- pi.mi7.250bp[which(pi.mi7.250bp$pop=='CV2'),]
co1pi.mi7.250bp <- read.table('pi.scaffold-mi7.250bp.co1.cnvMask.txt',header=T)
co2pi.mi7.250bp <- read.table('pi.scaffold-mi7.250bp.co2.cnvMask.txt',header=T)
capi.mi7.250bp <- pi.mi7.250bp[which(pi.mi7.250bp$pop=='CA'),]
crpi.mi7.250bp <- pi.mi7.250bp[which(pi.mi7.250bp$pop=='CR'),]
# dxy
cv1cv2.dxy.mi7.250bp <- dxy.mi7.250bp[which(dxy.mi7.250bp$pop1=='CV1' & dxy.mi7.250bp$pop2=='CV2'),]
cv1co1.dxy.mi7.250bp <- read.table('dxy.scaffold-mi7.250bp.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi7.250bp <- read.table('dxy.scaffold-mi7.250bp.co1co2.cnvMask.txt',header=T)
cv2co2.dxy.mi7.250bp <- dxy.mi7.250bp[which(dxy.mi7.250bp$pop1=='CV2' & dxy.mi7.250bp$pop2=='CO2'),]
cv2co1.dxy.mi7.250bp <- dxy.mi7.250bp[which(dxy.mi7.250bp$pop1=='CV2' & dxy.mi7.250bp$pop2=='CO1'),]
cv1ca.dxy.mi7.250bp <- dxy.mi7.250bp[which(dxy.mi7.250bp$pop1=='CA' & dxy.mi7.250bp$pop2=='CV1'),]
cv2ca.dxy.mi7.250bp <- dxy.mi7.250bp[which(dxy.mi7.250bp$pop1=='CA' & dxy.mi7.250bp$pop2=='CV2'),]
co1ca.dxy.mi7.250bp <- dxy.mi7.250bp[which(dxy.mi7.250bp$pop1=='CA' & dxy.mi7.250bp$pop2=='CO1'),]
co2ca.dxy.mi7.250bp <- dxy.mi7.250bp[which(dxy.mi7.250bp$pop1=='CA' & dxy.mi7.250bp$pop2=='CO2'),]
# Fst
cv1cv2.fst.mi7.250bp <- fst.mi7.250bp[which(fst.mi7.250bp$pop1=='CV1' & fst.mi7.250bp$pop2=='CV2'),]
cv1co1.fst.mi7.250bp <- read.table('fst.scaffold-mi7.250bp.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi7.250bp <- read.table('fst.scaffold-mi7.250bp.co1co2.cnvMask.txt',header=T)
cv2co2.fst.mi7.250bp <- fst.mi7.250bp[which(fst.mi7.250bp$pop1=='CV2' & fst.mi7.250bp$pop2=='CO2'),]

### Plot scans--------------------------------------------------------------

## pi
par(mfrow=c(6,1))
# 10 kb points, 100 kb windows
plot(cv1pi.mi7.10kb$window_pos_1,cv1pi.mi7.10kb$avg_pi,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(cv1pi.mi7.100kb$window_pos_1,cv1pi.mi7.100kb$avg_pi,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
cv1pi.mi7.avg <- mean(cv1pi.mi7.10kb$avg_pi,na.rm=T)
abline(h=cv1pi.mi7.avg,lty=3,col='blue')
abline(h=cv1pi.avg,lty=3,col='red')

plot(cv2pi.mi7.10kb$window_pos_1,cv2pi.mi7.10kb$avg_pi,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(cv2pi.mi7.100kb$window_pos_1,cv2pi.mi7.100kb$avg_pi,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
cv2pi.mi7.avg <- mean(cv2pi.mi7.10kb$avg_pi,na.rm=T)
abline(h=cv2pi.mi7.avg,lty=3,col='blue')
abline(h=cv2pi.avg,lty=3,col='red')

plot(co1pi.mi7.10kb$window_pos_1,co1pi.mi7.10kb$avg_pi,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(co1pi.mi7.100kb$window_pos_1,co1pi.mi7.100kb$avg_pi,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
co1pi.mi7.avg <- mean(co1pi.mi7.10kb$avg_pi,na.rm=T)
abline(h=co1pi.mi7.avg,lty=3,col='blue')
abline(h=co1pi.avg,lty=3,col='red')

plot(co2pi.mi7.10kb$window_pos_1,co2pi.mi7.10kb$avg_pi,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='pi')
lines(co2pi.mi7.100kb$window_pos_1,co2pi.mi7.100kb$avg_pi,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
co2pi.mi7.avg <- mean(co2pi.mi7.10kb$avg_pi,na.rm=T)
abline(h=co2pi.mi7.avg,lty=3,col='blue')
abline(h=co2pi.avg,lty=3,col='red')

# Zoomed, 250 bp points, 1 kb windows
plot(cv1pi.mi7.250bp$window_pos_1,cv1pi.mi7.250bp$avg_pi,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.002,0.025),ylab='pi',xlab='Chromosome Position')
lines(cv1pi.mi7.1kb$window_pos_1,cv1pi.mi7.1kb$avg_pi,type='l')
vert <- -0.0015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=cv1pi.mi7.avg,lty=3,col='blue')
abline(h=cv1pi.avg,lty=3,col='red')

plot(cv2pi.mi7.250bp$window_pos_1,cv2pi.mi7.250bp$avg_pi,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.002,0.025),ylab='pi',xlab='Chromosome Position')
lines(cv2pi.mi7.1kb$window_pos_1,cv2pi.mi7.1kb$avg_pi,type='l')
vert <- -0.0015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=cv2pi.mi7.avg,lty=3,col='blue')
abline(h=cv2pi.avg,lty=3,col='red')

plot(co1pi.mi7.250bp$window_pos_1,co1pi.mi7.250bp$avg_pi,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.002,0.025),ylab='pi',xlab='Chromosome Position')
lines(co1pi.mi7.1kb$window_pos_1,co1pi.mi7.1kb$avg_pi,type='l')
vert <- -0.0015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=co1pi.mi7.avg,lty=3,col='blue')
abline(h=co1pi.avg,lty=3,col='red')

plot(co2pi.mi7.250bp$window_pos_1,co2pi.mi7.250bp$avg_pi,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.002,0.025),ylab='pi',xlab='Chromosome Position')
lines(co2pi.mi7.1kb$window_pos_1,co2pi.mi7.1kb$avg_pi,type='l')
vert <- -0.0015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=co2pi.mi7.avg,lty=3,col='blue')
abline(h=co2pi.avg,lty=3,col='red')

## dxy
par(mfrow=c(3,1))
# 10 kb points, 100 kb windows
plot(cv1cv2.dxy.mi7.10kb$window_pos_1,cv1cv2.dxy.mi7.10kb$avg_dxy,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(cv1cv2.dxy.mi7.100kb$window_pos_1,cv1cv2.dxy.mi7.100kb$avg_dxy,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
cv1cv2.dxy.mi7.avg <- mean(cv1cv2.dxy.mi7.10kb$avg_dxy,na.rm=T)
abline(h=cv1cv2.dxy.mi7.avg,lty=3,col='blue')
abline(h=cv1cv2.dxy.avg,lty=3,col='red')

plot(co1co2.dxy.mi7.10kb$window_pos_1,co1co2.dxy.mi7.10kb$avg_dxy,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(co1co2.dxy.mi7.100kb$window_pos_1,co1co2.dxy.mi7.100kb$avg_dxy,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
co1co2.dxy.mi7.avg <- mean(co1co2.dxy.mi7.10kb$avg_dxy,na.rm=T)
abline(h=co1co2.dxy.mi7.avg,lty=3,col='blue')
abline(h=co1co2.dxy.avg,lty=3,col='red')

plot(cv1co1.dxy.mi7.10kb$window_pos_1,cv1co1.dxy.mi7.10kb$avg_dxy,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='dxy')
lines(cv1co1.dxy.mi7.100kb$window_pos_1,cv1co1.dxy.mi7.100kb$avg_dxy,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
cv1co1.dxy.mi7.avg <- mean(cv1co1.dxy.mi7.10kb$avg_dxy,na.rm=T)
abline(h=cv1co1.dxy.mi7.avg,lty=3,col='blue')
abline(h=cv1co1.dxy.avg,lty=3,col='red')

# Zoomed, 250 bp points, 10 kb windows
plot(cv1cv2.dxy.mi7.250bp$window_pos_1,cv1cv2.dxy.mi7.250bp$avg_dxy,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.002,0.045),ylab='dxy',xlab='Chromosome Position')
lines(cv1cv2.dxy.mi7.1kb$window_pos_1,cv1cv2.dxy.mi7.1kb$avg_dxy,type='l')
vert <- -0.0015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=cv1cv2.dxy.mi7.avg,lty=3,col='blue')
abline(h=cv1cv2.dxy.avg,lty=3,col='red')

plot(co1co2.dxy.mi7.250bp$window_pos_1,co1co2.dxy.mi7.250bp$avg_dxy,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.002,0.045),ylab='dxy',xlab='Chromosome Position')
lines(co1co2.dxy.mi7.1kb$window_pos_1,co1co2.dxy.mi7.1kb$avg_dxy,type='l')
vert <- -0.0015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=co1co2.dxy.mi7.avg,lty=3,col='blue')
abline(h=co1co2.dxy.avg,lty=3,col='red')

plot(cv1co1.dxy.mi7.250bp$window_pos_1,cv1co1.dxy.mi7.250bp$avg_dxy,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.002,0.045),ylab='dxy',xlab='Chromosome Position')
lines(cv1co1.dxy.mi7.1kb$window_pos_1,cv1co1.dxy.mi7.1kb$avg_dxy,type='l')
vert <- -0.0015
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.dxy.mi7.avg,lty=3,col='blue')
abline(h=cv1co1.dxy.avg,lty=3,col='red')

## fst
par(mfrow=c(3,1))
# 10 kb points, 100 kb windows
plot(cv1cv2.fst.mi7.10kb$window_pos_1,cv1cv2.fst.mi7.10kb$avg_wc_fst,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(cv1cv2.fst.mi7.100kb$window_pos_1,cv1cv2.fst.mi7.100kb$avg_wc_fst,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
cv1cv2.fst.mi7.avg <- mean(cv1cv2.fst.mi7.10kb$avg_wc_fst,na.rm=T)
abline(h=cv1cv2.fst.mi7.avg,lty=3,col='blue')
abline(h=cv1cv2.fst.avg,lty=3,col='red')

plot(co1co2.fst.mi7.10kb$window_pos_1,co1co2.fst.mi7.10kb$avg_wc_fst,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(co1co2.fst.mi7.100kb$window_pos_1,co1co2.fst.mi7.100kb$avg_wc_fst,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
co1co2.fst.mi7.avg <- mean(co1co2.fst.mi7.10kb$avg_wc_fst,na.rm=T)
abline(h=co1co2.fst.mi7.avg,lty=3,col='blue')
abline(h=co1co2.fst.avg,lty=3,col='red')

plot(cv1co1.fst.mi7.10kb$window_pos_1,cv1co1.fst.mi7.10kb$avg_wc_fst,pch=20,col=alpha('maroon',0.15),xlab='Chromosome Position (Mb)',ylab='fst')
lines(cv1co1.fst.mi7.100kb$window_pos_1,cv1co1.fst.mi7.100kb$avg_wc_fst,type='l')
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')
cv1co1.fst.mi7.avg <- mean(cv1co1.fst.mi7.10kb$avg_wc_fst,na.rm=T)
abline(h=cv1co1.fst.mi7.avg,lty=3,col='blue')
abline(h=cv1co1.fst.avg,lty=3,col='red')

# Zoomed, 1kb points, 10 kb windows
plot(cv1cv2.fst.mi7.250bp$window_pos_1,cv1cv2.fst.mi7.250bp$avg_wc_fst,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(cv1cv2.fst.mi7.1kb$window_pos_1,cv1cv2.fst.mi7.1kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=cv1cv2.fst.mi7.avg,lty=3,col='blue')
abline(h=cv1cv2.fst.avg,lty=3,col='red')

plot(co1co2.fst.mi7.250bp$window_pos_1,co1co2.fst.mi7.250bp$avg_wc_fst,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(co1co2.fst.mi7.1kb$window_pos_1,co1co2.fst.mi7.1kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=co1co2.fst.mi7.avg,lty=3,col='blue')
abline(h=co1co2.fst.avg,lty=3,col='red')

plot(cv1co1.fst.mi7.250bp$window_pos_1,cv1co1.fst.mi7.250bp$avg_wc_fst,pch=20,col=alpha('maroon',0.15),xlim=c(2980000,3080000),ylim=c(-0.1,1),ylab='fst',xlab='Chromosome Position')
lines(cv1co1.fst.mi7.1kb$window_pos_1,cv1co1.fst.mi7.1kb$avg_wc_fst,type='l')
vert <- -0.075
for (g in pla2) {
  segments(x0=pla2$V2,y0=vert,x1=pla2$V3,y1=vert,col='darkblue')
}
abline(h=cv1co1.fst.mi7.avg,lty=3,col='blue')
abline(h=cv1co1.fst.avg,lty=3,col='red')
