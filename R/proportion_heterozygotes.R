############################################################################
# Proportion of heterozygotes in venom gene regions
############################################################################

### Goal: quantify the proportion of heterozygotes in each population across
### venom gene regions, compared with the chromosome background. 

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & load dependencies-------------------------------

setwd('./heterozygosity/unphased/')

library(scales)
library(data.table)
library(dplyr)
library(Rmisc)

options('stringsAsFactors'=FALSE)

### Read in general coordinates---------------------------------------------

# Region locations
svmp.reg <- read.table('./resources/region_SVMP_scaffold-mi1.bed',header=F)
svsp.reg <- read.table('./resources/region_SVSP_scaffold-mi2.bed',header=F)
pla2.reg <- read.table('./resources/region_PLA2_scaffold-mi7.bed',header=F)

# Gene locations
svmp <- read.table('./resources/gene_SVMP.bed',header=F)
svsp <- read.table('./resources/gene_SVSP.bed',header=F)
pla2 <- read.table('./resources/gene_PLA2.bed',header=F)

### Read in data------------------------------------------------------------

cv1.mi1.100kb <- read.table('het.cv1.scaffold-mi1.100kb.txt',header=T)
cv1.mi1.10kb <- read.table('het.cv1.scaffold-mi1.10kb.txt',header=T)
cv1.mi1.1kb <- read.table('het.cv1.scaffold-mi1.1kb.txt',header=T)
cv1.mi1.all <- read.table('het.cv1.scaffold-mi1.txt',header=T)
cv1.mi2.100kb <- read.table('het.cv1.scaffold-mi2.100kb.txt',header=T)
cv1.mi2.10kb <- read.table('het.cv1.scaffold-mi2.10kb.txt',header=T)
cv1.mi2.1kb <- read.table('het.cv1.scaffold-mi2.1kb.txt',header=T)
cv1.mi2.all <- read.table('het.cv1.scaffold-mi2.txt',header=T)
cv1.mi7.100kb <- read.table('het.cv1.scaffold-mi7.100kb.txt',header=T)
cv1.mi7.10kb <- read.table('het.cv1.scaffold-mi7.10kb.txt',header=T)
cv1.mi7.1kb <- read.table('het.cv1.scaffold-mi7.1kb.txt',header=T)
cv1.mi7.all <- read.table('het.cv1.scaffold-mi7.txt',header=T)

cv2.mi1.100kb <- read.table('het.cv2.scaffold-mi1.100kb.txt',header=T)
cv2.mi1.10kb <- read.table('het.cv2.scaffold-mi1.10kb.txt',header=T)
cv2.mi1.1kb <- read.table('het.cv2.scaffold-mi1.1kb.txt',header=T)
cv2.mi1.all <- read.table('het.cv2.scaffold-mi1.txt',header=T)
cv2.mi2.100kb <- read.table('het.cv2.scaffold-mi2.100kb.txt',header=T)
cv2.mi2.10kb <- read.table('het.cv2.scaffold-mi2.10kb.txt',header=T)
cv2.mi2.1kb <- read.table('het.cv2.scaffold-mi2.1kb.txt',header=T)
cv2.mi2.all <- read.table('het.cv2.scaffold-mi2.txt',header=T)
cv2.mi7.100kb <- read.table('het.cv2.scaffold-mi7.100kb.txt',header=T)
cv2.mi7.10kb <- read.table('het.cv2.scaffold-mi7.10kb.txt',header=T)
cv2.mi7.1kb <- read.table('het.cv2.scaffold-mi7.1kb.txt',header=T)
cv2.mi7.all <- read.table('het.cv2.scaffold-mi7.txt',header=T)

co1.mi1.100kb <- read.table('het.co1.scaffold-mi1.100kb.txt',header=T)
co1.mi1.10kb <- read.table('het.co1.scaffold-mi1.10kb.cnvMask.txt',header=T)
co1.mi1.1kb <- read.table('het.co1.scaffold-mi1.1kb.cnvMask.txt',header=T)
co1.mi1.all <- read.table('het.co1.scaffold-mi1.txt',header=T)
co1.mi2.100kb <- read.table('het.co1.scaffold-mi2.100kb.txt',header=T)
co1.mi2.10kb <- read.table('het.co1.scaffold-mi2.10kb.cnvMask.txt',header=T)
co1.mi2.1kb <- read.table('het.co1.scaffold-mi2.1kb.cnvMask.txt',header=T)
co1.mi2.all <- read.table('het.co1.scaffold-mi2.txt',header=T)
co1.mi7.100kb <- read.table('het.co1.scaffold-mi7.100kb.txt',header=T)
co1.mi7.10kb <- read.table('het.co1.scaffold-mi7.10kb.cnvMask.txt',header=T)
co1.mi7.1kb <- read.table('het.co1.scaffold-mi7.1kb.cnvMask.txt',header=T)
co1.mi7.all <- read.table('het.co1.scaffold-mi7.txt',header=T)

co2.mi1.100kb <- read.table('het.co2.scaffold-mi1.100kb.txt',header=T)
co2.mi1.10kb <- read.table('het.co2.scaffold-mi1.10kb.cnvMask.txt',header=T)
co2.mi1.1kb <- read.table('het.co2.scaffold-mi1.1kb.cnvMask.txt',header=T)
co2.mi1.all <- read.table('het.co2.scaffold-mi1.txt',header=T)
co2.mi2.100kb <- read.table('het.co2.scaffold-mi2.100kb.txt',header=T)
co2.mi2.10kb <- read.table('het.co2.scaffold-mi2.10kb.cnvMask.txt',header=T)
co2.mi2.1kb <- read.table('het.co2.scaffold-mi2.1kb.cnvMask.txt',header=T)
co2.mi2.all <- read.table('het.co2.scaffold-mi2.txt',header=T)
co2.mi7.100kb <- read.table('het.co2.scaffold-mi7.100kb.txt',header=T)
co2.mi7.10kb <- read.table('het.co2.scaffold-mi7.10kb.cnvMask.txt',header=T)
co2.mi7.1kb <- read.table('het.co2.scaffold-mi7.1kb.cnvMask.txt',header=T)
co2.mi7.all <- read.table('het.co2.scaffold-mi7.txt',header=T)

### Parse venom from background---------------------------------------------

cv1.svmp.1kb <- cv1.mi1.1kb[which(cv1.mi1.1kb$start>=svmp.reg$V2 & cv1.mi1.1kb$start<=svmp.reg$V3),]
cv2.svmp.1kb <- cv2.mi1.1kb[which(cv2.mi1.1kb$start>=svmp.reg$V2 & cv2.mi1.1kb$start<=svmp.reg$V3),]
co1.svmp.1kb <- co1.mi1.1kb[which(co1.mi1.1kb$start>=svmp.reg$V2 & co1.mi1.1kb$start<=svmp.reg$V3),]
co2.svmp.1kb <- co2.mi1.1kb[which(co2.mi1.1kb$start>=svmp.reg$V2 & co2.mi1.1kb$start<=svmp.reg$V3),]

cv1.mi1.bk <- cv1.mi1.1kb[which(cv1.mi1.1kb$start<svmp.reg$V2 | cv1.mi1.1kb$start>svmp.reg$V3),]
cv2.mi1.bk <- cv2.mi1.1kb[which(cv2.mi1.1kb$start<svmp.reg$V2 | cv2.mi1.1kb$start>svmp.reg$V3),]
co1.mi1.bk <- co1.mi1.1kb[which(co1.mi1.1kb$start<svmp.reg$V2 | co1.mi1.1kb$start>svmp.reg$V3),]
co2.mi1.bk <- co2.mi1.1kb[which(co2.mi1.1kb$start<svmp.reg$V2 | co2.mi1.1kb$start>svmp.reg$V3),]

cv1.svsp.1kb <- cv1.mi2.1kb[which(cv1.mi2.1kb$start>=svsp.reg$V2 & cv1.mi2.1kb$start<=svsp.reg$V3),]
cv2.svsp.1kb <- cv2.mi2.1kb[which(cv2.mi2.1kb$start>=svsp.reg$V2 & cv2.mi2.1kb$start<=svsp.reg$V3),]
co1.svsp.1kb <- co1.mi2.1kb[which(co1.mi2.1kb$start>=svsp.reg$V2 & co1.mi2.1kb$start<=svsp.reg$V3),]
co2.svsp.1kb <- co2.mi2.1kb[which(co2.mi2.1kb$start>=svsp.reg$V2 & co2.mi2.1kb$start<=svsp.reg$V3),]

cv1.mi2.bk <- cv1.mi2.1kb[which(cv1.mi2.1kb$start<svsp.reg$V2 | cv1.mi2.1kb$start>svsp.reg$V3),]
cv2.mi2.bk <- cv2.mi2.1kb[which(cv2.mi2.1kb$start<svsp.reg$V2 | cv2.mi2.1kb$start>svsp.reg$V3),]
co1.mi2.bk <- co1.mi2.1kb[which(co1.mi2.1kb$start<svsp.reg$V2 | co1.mi2.1kb$start>svsp.reg$V3),]
co2.mi2.bk <- co2.mi2.1kb[which(co2.mi2.1kb$start<svsp.reg$V2 | co2.mi2.1kb$start>svsp.reg$V3),]

cv1.pla2.1kb <- cv1.mi7.1kb[which(cv1.mi7.1kb$start>=pla2.reg$V2 & cv1.mi7.1kb$start<=pla2.reg$V3),]
cv2.pla2.1kb <- cv2.mi7.1kb[which(cv2.mi7.1kb$start>=pla2.reg$V2 & cv2.mi7.1kb$start<=pla2.reg$V3),]
co1.pla2.1kb <- co1.mi7.1kb[which(co1.mi7.1kb$start>=pla2.reg$V2 & co1.mi7.1kb$start<=pla2.reg$V3),]
co2.pla2.1kb <- co2.mi7.1kb[which(co2.mi7.1kb$start>=pla2.reg$V2 & co2.mi7.1kb$start<=pla2.reg$V3),]

cv1.mi7.bk <- cv1.mi7.1kb[which(cv1.mi7.1kb$start<pla2.reg$V2 | cv1.mi7.1kb$start>pla2.reg$V3),]
cv2.mi7.bk <- cv2.mi7.1kb[which(cv2.mi7.1kb$start<pla2.reg$V2 | cv2.mi7.1kb$start>pla2.reg$V3),]
co1.mi7.bk <- co1.mi7.1kb[which(co1.mi7.1kb$start<pla2.reg$V2 | co1.mi7.1kb$start>pla2.reg$V3),]
co2.mi7.bk <- co2.mi7.1kb[which(co2.mi7.1kb$start<pla2.reg$V2 | co2.mi7.1kb$start>pla2.reg$V3),]

### Compare venom regions to background-------------------------------------

# SVMP
par(mfrow=c(1,4))
boxplot(as.numeric(cv1.svmp.1kb$prop_het),as.numeric(cv1.mi1.bk$prop_het))
boxplot(as.numeric(cv2.svmp.1kb$prop_het),as.numeric(cv2.mi1.bk$prop_het))
boxplot(as.numeric(co1.svmp.1kb$prop_het),as.numeric(co1.mi1.bk$prop_het))
boxplot(as.numeric(co2.svmp.1kb$prop_het),as.numeric(co2.mi1.bk$prop_het))

# SVSP
boxplot(as.numeric(cv1.svsp.1kb$prop_het),as.numeric(cv1.mi2.bk$prop_het))
boxplot(as.numeric(cv2.svsp.1kb$prop_het),as.numeric(cv2.mi2.bk$prop_het))
boxplot(as.numeric(co1.svsp.1kb$prop_het),as.numeric(co1.mi2.bk$prop_het))
boxplot(as.numeric(co2.svsp.1kb$prop_het),as.numeric(co2.mi2.bk$prop_het))

# PLA2
boxplot(as.numeric(cv1.pla2.1kb$prop_het),as.numeric(cv1.mi7.bk$prop_het))
boxplot(as.numeric(cv2.pla2.1kb$prop_het),as.numeric(cv2.mi7.bk$prop_het))
boxplot(as.numeric(co1.pla2.1kb$prop_het),as.numeric(co1.mi7.bk$prop_het))
boxplot(as.numeric(co2.pla2.1kb$prop_het),as.numeric(co2.mi7.bk$prop_het))

### Summary statistics------------------------------------------------------

mean(as.numeric(cv1.svmp.1kb$prop_het),na.rm=T)
mean(as.numeric(cv1.mi1.bk$prop_het),na.rm=T)
mean(as.numeric(cv2.svmp.1kb$prop_het),na.rm=T)
mean(as.numeric(cv2.mi1.bk$prop_het),na.rm=T)
mean(as.numeric(co1.svmp.1kb$prop_het),na.rm=T)
mean(as.numeric(co1.mi1.bk$prop_het),na.rm=T)
mean(as.numeric(co2.svmp.1kb$prop_het),na.rm=T)
mean(as.numeric(co2.mi1.bk$prop_het),na.rm=T)

mean(as.numeric(cv1.svsp.1kb$prop_het),na.rm=T)
mean(as.numeric(cv1.mi2.bk$prop_het),na.rm=T)
mean(as.numeric(cv2.svsp.1kb$prop_het),na.rm=T)
mean(as.numeric(cv2.mi2.bk$prop_het),na.rm=T)
mean(as.numeric(co1.svsp.1kb$prop_het),na.rm=T)
mean(as.numeric(co1.mi2.bk$prop_het),na.rm=T)
mean(as.numeric(co2.svsp.1kb$prop_het),na.rm=T)
mean(as.numeric(co2.mi2.bk$prop_het),na.rm=T)

mean(as.numeric(cv1.pla2.1kb$prop_het),na.rm=T)
mean(as.numeric(cv1.mi7.bk$prop_het),na.rm=T)
mean(as.numeric(cv2.pla2.1kb$prop_het),na.rm=T)
mean(as.numeric(cv2.mi7.bk$prop_het),na.rm=T)
mean(as.numeric(co1.pla2.1kb$prop_het),na.rm=T)
mean(as.numeric(co1.mi7.bk$prop_het),na.rm=T)
mean(as.numeric(co2.pla2.1kb$prop_het),na.rm=T)
mean(as.numeric(co2.mi7.bk$prop_het),na.rm=T)

wilcox.test(as.numeric(cv1.svmp.1kb$prop_het),as.numeric(cv1.mi1.bk$prop_het))
wilcox.test(as.numeric(cv2.svmp.1kb$prop_het),as.numeric(cv2.mi1.bk$prop_het))
wilcox.test(as.numeric(co1.svmp.1kb$prop_het),as.numeric(co1.mi1.bk$prop_het))
wilcox.test(as.numeric(co2.svmp.1kb$prop_het),as.numeric(co2.mi1.bk$prop_het))

wilcox.test(as.numeric(cv1.svsp.1kb$prop_het),as.numeric(cv1.mi2.bk$prop_het))
wilcox.test(as.numeric(cv2.svsp.1kb$prop_het),as.numeric(cv2.mi2.bk$prop_het))
wilcox.test(as.numeric(co1.svsp.1kb$prop_het),as.numeric(co1.mi2.bk$prop_het))
wilcox.test(as.numeric(co2.svsp.1kb$prop_het),as.numeric(co2.mi2.bk$prop_het))

wilcox.test(as.numeric(cv1.pla2.1kb$prop_het),as.numeric(cv1.mi7.bk$prop_het))
wilcox.test(as.numeric(cv2.pla2.1kb$prop_het),as.numeric(cv2.mi7.bk$prop_het))
wilcox.test(as.numeric(co1.pla2.1kb$prop_het),as.numeric(co1.mi7.bk$prop_het))
wilcox.test(as.numeric(co2.pla2.1kb$prop_het),as.numeric(co2.mi7.bk$prop_het))

### Test for enrichment of high prop. het-----------------------------------

# Set up data for comparisons
cv1.svmp.all <- read.table('het.beta.cv1.svmp.txt',header=T,row.names=NULL)
#cv2.svmp.all <- read.table('het.beta.cv2.svmp.txt',header=T,row.names=NULL)
co1.svmp.all <- read.table('het.beta.co1.svmp.txt',header=T,row.names=NULL)
#co2.svmp.all <- read.table('het.beta.co2.svmp.txt',header=T,row.names=NULL)

cv1.mi1.bk <- read.table('het.beta.cv1.scaffold-mi1.gene.txt',header=T,row.names=NULL)
#cv2.mi1.bk <- read.table('het.beta.cv2.scaffold-mi1.gene.txt',header=T,row.names=NULL)
co1.mi1.bk <- read.table('het.beta.co1.scaffold-mi1.gene.txt',header=T,row.names=NULL)
#co2.mi1.bk <- read.table('het.beta.co2.scaffold-mi1.gene.txt',header=T,row.names=NULL)

cv1.svsp.all <- read.table('het.beta.cv1.svsp.txt',header=T,row.names=NULL)
#cv2.svsp.all <- read.table('het.beta.cv2.svsp.txt',header=T,row.names=NULL)
co1.svsp.all <- read.table('het.beta.co1.svsp.txt',header=T,row.names=NULL)
#co2.svsp.all <- read.table('het.beta.co2.svsp.txt',header=T,row.names=NULL)

cv1.mi2.bk <- read.table('het.beta.cv1.scaffold-mi2.gene.txt',header=T,row.names=NULL)
#cv2.mi2.bk <- read.table('het.beta.cv2.scaffold-mi2.gene.txt',header=T,row.names=NULL)
co1.mi2.bk <- read.table('het.beta.co1.scaffold-mi2.gene.txt',header=T,row.names=NULL)
#co2.mi2.bk <- read.table('het.beta.co2.scaffold-mi2.gene.txt',header=T,row.names=NULL)

cv1.pla2.all <- read.table('het.beta.cv1.pla2.txt',header=T,row.names=NULL)
#cv2.pla2.all <- read.table('het.beta.cv2.pla2.txt',header=T,row.names=NULL)
co1.pla2.all <- read.table('het.beta.co1.pla2.txt',header=T,row.names=NULL)
#co2.pla2.all <- read.table('het.beta.co2.pla2.txt',header=T,row.names=NULL)

cv1.mi7.bk <- read.table('het.beta.cv1.scaffold-mi7.gene.txt',header=T,row.names=NULL)
#cv2.mi7.bk <- read.table('het.beta.cv2.scaffold-mi7.gene.txt',header=T,row.names=NULL)
co1.mi7.bk <- read.table('het.beta.co1.scaffold-mi7.gene.txt',header=T,row.names=NULL)
#co2.mi7.bk <- read.table('het.beta.co2.scaffold-mi7.gene.txt',header=T,row.names=NULL)

# Count numbers of sites with proportion of heterozygotes = 1

# SVMP
count(cv1.svmp.all$prop_het>=0.9)
count(cv1.mi1.bk$prop_het>=0.9)
svmp.cv1.comp <- matrix(c(226,969,30,3578),nrow = 2)

#count(cv2.svmp.all$prop_het==1)
#count(cv2.mi1.bk$prop_het==1)
#svmp.cv2.comp <- matrix(c(250,21560,160,240532),nrow = 2)

count(co1.svmp.all$prop_het>=0.9)
count(co1.mi1.bk$prop_het>=0.9)
svmp.co1.comp <- matrix(c(76,196,308,2895),nrow = 2)

#count(co2.svmp.all$prop_het==1)
#count(co2.mi1.bk$prop_het==1)
#svmp.co2.comp <- matrix(c(153,5682,645,251950),nrow = 2)

# SVSP
count(cv1.svsp.all$prop_het>=0.9)
count(cv1.mi2.bk$prop_het>=0.9)
svsp.cv1.comp <- matrix(c(51,295,94,3822),nrow = 2)

#count(cv2.svsp.all$prop_het==1)
#count(cv2.mi2.bk$prop_het==1)
#svsp.cv2.comp <- matrix(c(250,21560,71,206173),nrow = 2)

count(co1.svsp.all$prop_het>=0.9)
count(co1.mi2.bk$prop_het>=0.9)
svsp.co1.comp <- matrix(c(28,72,94,2402),nrow = 2)

#count(co2.svsp.all$prop_het==1)
#count(co2.mi2.bk$prop_het==1)
#svsp.co2.comp <- matrix(c(153,5682,136,206773),nrow = 2)

# PLA2
count(cv1.pla2.all$prop_het>=0.9)
count(cv1.mi7.bk$prop_het>=0.9)
pla2.cv1.comp <- matrix(c(0,25,43,2815),nrow = 2)

#count(cv2.pla2.all$prop_het==1)
#count(cv2.mi7.bk$prop_het==1)
#pla2.cv2.comp <- matrix(c(0,695,56,248823),nrow = 2)

count(co1.pla2.all$prop_het>=0.9)
count(co1.mi7.bk$prop_het>=0.9)
pla2.co1.comp <- matrix(c(0,6,39,2369),nrow = 2)

#count(co2.pla2.all$prop_het==1)
#count(co2.mi7.bk$prop_het==1)
#pla2.co2.comp <- matrix(c(4,194,129,250770),nrow = 2)

# Perform Fisher's exact tests
fisher.test(svmp.cv1.comp)
#fisher.test(svmp.cv2.comp)
fisher.test(svmp.co1.comp)
#fisher.test(svmp.co2.comp)

fisher.test(svsp.cv1.comp)
#fisher.test(svsp.cv2.comp)
fisher.test(svsp.co1.comp)
#fisher.test(svsp.co2.comp)

fisher.test(pla2.cv1.comp)
#fisher.test(pla2.cv2.comp)
fisher.test(pla2.co1.comp)
#fisher.test(pla2.co2.comp)

mean(cv1.svmp.all$prop_het)
mean(cv1.svsp.all$prop_het)
mean(cv1.pla2.all$prop_het)

mean(cv1.mi1.bk$prop_het)
mean(cv1.mi2.bk$prop_het)

mean(co1.svmp.all$prop_het)
mean(co1.svsp.all$prop_het)
mean(co1.pla2.all$prop_het)


### Plot scans--------------------------------------------------------------

# The zoomed out scans are largely useless, since they are not CNV-masked in
# CO1 and CO2

par(mfrow=c(3,1))
plot(cv1.mi1.10kb$start,cv1.mi1.10kb$prop_het,pch=20,col=alpha('lightblue3',0.25))
lines(cv1.mi1.100kb$start,cv1.mi1.100kb$prop_het)
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(cv1.mi2.10kb$start,cv1.mi2.10kb$prop_het,pch=20,col=alpha('aquamarine3',0.25))
lines(cv1.mi2.100kb$start,cv1.mi2.100kb$prop_het)
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(cv1.mi7.10kb$start,cv1.mi7.10kb$prop_het,pch=20,col=alpha('maroon',0.25))
lines(cv1.mi7.100kb$start,cv1.mi7.100kb$prop_het)
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')

plot(cv2.mi1.10kb$start,cv2.mi1.10kb$prop_het,pch=20,col=alpha('lightblue3',0.25))
lines(cv2.mi1.100kb$start,cv2.mi1.100kb$prop_het)
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(cv2.mi2.10kb$start,cv2.mi2.10kb$prop_het,pch=20,col=alpha('aquamarine3',0.25))
lines(cv2.mi2.100kb$start,cv2.mi2.100kb$prop_het)
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(cv2.mi7.10kb$start,cv2.mi7.10kb$prop_het,pch=20,col=alpha('maroon',0.25))
lines(cv2.mi7.100kb$start,cv2.mi7.100kb$prop_het)
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')

plot(co1.mi1.10kb$start,co1.mi1.10kb$prop_het,pch=20,col=alpha('lightblue3',0.25))
lines(co1.mi1.100kb$start,co1.mi1.100kb$prop_het)
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(co1.mi2.10kb$start,co1.mi2.10kb$prop_het,pch=20,col=alpha('aquamarine3',0.25))
lines(co1.mi2.100kb$start,co1.mi2.100kb$prop_het)
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(co1.mi7.10kb$start,co1.mi7.10kb$prop_het,pch=20,col=alpha('maroon',0.25))
lines(co1.mi7.100kb$start,co1.mi7.100kb$prop_het)
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')

plot(co2.mi1.10kb$start,co2.mi1.10kb$prop_het,pch=20,col=alpha('lightblue3',0.25))
lines(co2.mi1.100kb$start,co2.mi1.100kb$prop_het)
abline(v=svmp.reg$V2,col='grey4')
abline(v=svmp.reg$V3,col='grey4')
plot(co2.mi2.10kb$start,co2.mi2.10kb$prop_het,pch=20,col=alpha('aquamarine3',0.25))
lines(co2.mi2.100kb$start,co2.mi2.100kb$prop_het)
abline(v=svsp.reg$V2,col='grey4')
abline(v=svsp.reg$V3,col='grey4')
plot(co2.mi7.10kb$start,co2.mi7.10kb$prop_het,pch=20,col=alpha('maroon',0.25))
lines(co2.mi7.100kb$start,co2.mi7.100kb$prop_het)
abline(v=pla2.reg$V2,col='grey4')
abline(v=pla2.reg$V3,col='grey4')

# Zoomed, 1 kb, 10 kb windows
plot(cv1.mi1.1kb$start,cv1.mi1.1kb$prop_het,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.35),ylab='Heterozygosity',xlab='Chromosome Position')
lines(cv1.mi1.10kb$start,cv1.mi1.10kb$prop_het,type='l')
vert <- -0.02
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

plot(cv2.mi1.1kb$start,cv2.mi1.1kb$prop_het,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,0.35),ylab='Heterozygosity',xlab='Chromosome Position')
lines(cv2.mi1.10kb$start,cv2.mi1.10kb$prop_het,type='l')
vert <- -0.02
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

plot(co1.mi1.1kb$start,co1.mi1.1kb$prop_het,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,1),ylab='Heterozygosity',xlab='Chromosome Position')
lines(co1.mi1.10kb$start,co1.mi1.10kb$prop_het,type='l')
vert <- -0.02
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

plot(co2.mi1.1kb$start,co2.mi1.1kb$prop_het,pch=20,col=alpha('lightblue3',0.15),xlim=c(13600000,14700000),ylim=c(-0.025,1),ylab='Heterozygosity',xlab='Chromosome Position')
lines(co2.mi1.10kb$start,co2.mi1.10kb$prop_het,type='l')
vert <- -0.02
for (g in svmp) {
  segments(x0=svmp$V2,y0=vert,x1=svmp$V3,y1=vert,col='darkblue')
}

