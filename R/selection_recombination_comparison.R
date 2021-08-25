############################################################################
# Comparison of recombination rate with selection statistics
############################################################################

### Goal: determine if there are genome-wide and venom-region-specific links
### between recombination rate, iHS, and ß statistics

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & load dependencies-------------------------------

setwd('./')
library(scales)
library(data.table)
library(dplyr)

### Read in general coordinates---------------------------------------------

# Gene locations
svmp <- read.table('./resources/gene_SVMP.bed',header=F)
svsp <- read.table('./resources/gene_SVSP.bed',header=F)
pla2 <- read.table('./resources/gene_PLA2.bed',header=F)

# Region locations
svmp.reg <- read.table('./resources/region_SVMP_scaffold-mi1.bed',header=F)
svsp.reg <- read.table('./resources/region_SVSP_scaffold-mi2.bed',header=F)
pla2.reg <- read.table('./resources/region_PLA2_scaffold-mi7.bed',header=F)

### Read in π-corrected recombination rates---------------------------------

r.cv <- read.table('./recombination/sliding_windows/viridis.recomb.bpen10.windowed.10kb.centromereMask.piCorrected.txt',header=T)
r.co <- read.table('./recombination/sliding_windows/oreganus.recomb.bpen10.windowed.10kb.centromereMask.cnvMask.piCorrected.txt',header=T)

r.mi1.cv <- r.cv[r.cv$chrom=='scaffold-mi1',]
r.mi2.cv <- r.cv[r.cv$chrom=='scaffold-mi2',]
r.mi7.cv <- read.table('./recombination/sliding_windows/viridis.recomb.bpen10.scaffold-mi7.windowed.1kb.centromereMask.piCorrected.txt',header=T)
r.mi1.co <- r.co[r.co$chrom=='scaffold-mi1',]
r.mi2.co <- r.co[r.co$chrom=='scaffold-mi2',]
r.mi7.co <- read.table('./recombination/sliding_windows/oreganus.recomb.bpen10.scaffold-mi7.windowed.1kb.centromereMask.cnvMask.piCorrected.txt',header=T)

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

### Read in iHS values--------------------------------------------------------

cv1.ihs.10kb <- read.table('./rehh/cv.all_ihs.10kb.txt',header=T)
co1.ihs.10kb <- read.table('./rehh/co.all_ihs.10kb.cnvMask.txt',header=T)
cv1.ihs.1kb <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
co1.ihs.1kb <- read.table('./rehh/co.all_ihs.1kb.cnvMask.txt',header=T)

cv1.ihs.mi1.10kb <- cv1.ihs.10kb[which(cv1.ihs.10kb$chrom=='scaffold-mi1'),]
co1.ihs.mi1.10kb <- co1.ihs.10kb[which(co1.ihs.10kb$chrom=='scaffold-mi1'),]
cv1.ihs.mi2.10kb <- cv1.ihs.10kb[which(cv1.ihs.10kb$chrom=='scaffold-mi2'),]
co1.ihs.mi2.10kb <- co1.ihs.10kb[which(co1.ihs.10kb$chrom=='scaffold-mi2'),]
cv1.ihs.mi7.1kb <- cv1.ihs.1kb[which(cv1.ihs.1kb$chrom=='scaffold-mi7'),]
co1.ihs.mi7.1kb <- co1.ihs.1kb[which(co1.ihs.1kb$chrom=='scaffold-mi7'),]

cv1.ihs.svmp <- cv1.ihs.mi1.10kb[which(cv1.ihs.mi1.10kb$start>=svmp.reg$V2 & cv1.ihs.mi1.10kb$end<=svmp.reg$V3),]
co1.ihs.svmp <- co1.ihs.mi1.10kb[which(co1.ihs.mi1.10kb$start>=svmp.reg$V2 & co1.ihs.mi1.10kb$end<=svmp.reg$V3),]
cv1.ihs.svsp <- cv1.ihs.mi2.10kb[which(cv1.ihs.mi2.10kb$start>=svsp.reg$V2 & cv1.ihs.mi2.10kb$end<=svsp.reg$V3),]
co1.ihs.svsp <- co1.ihs.mi2.10kb[which(co1.ihs.mi2.10kb$start>=svsp.reg$V2 & co1.ihs.mi2.10kb$end<=svsp.reg$V3),]
cv1.ihs.pla2 <- cv1.ihs.mi7.1kb[which(cv1.ihs.mi7.1kb$start>=pla2.reg$V2 & cv1.ihs.mi7.1kb$end<=pla2.reg$V3),]
co1.ihs.pla2 <- co1.ihs.mi7.1kb[which(co1.ihs.mi7.1kb$start>=pla2.reg$V2 & co1.ihs.mi7.1kb$end<=pla2.reg$V3),]

cv1.ihs.mi1.bk <- cv1.ihs.mi1.10kb[which(cv1.ihs.mi1.10kb$start<svmp.reg$V2 & cv1.ihs.mi1.10kb$end<svmp.reg$V2 | cv1.ihs.mi1.10kb$start>svmp.reg$V3 & cv1.ihs.mi1.10kb$end>svmp.reg$V3),]
co1.ihs.mi1.bk <- co1.ihs.mi1.10kb[which(co1.ihs.mi1.10kb$start<svmp.reg$V2 & co1.ihs.mi1.10kb$end<svmp.reg$V2 | co1.ihs.mi1.10kb$start>svmp.reg$V3 & co1.ihs.mi1.10kb$end>svmp.reg$V3),]
cv1.ihs.mi2.bk <- cv1.ihs.mi2.10kb[which(cv1.ihs.mi2.10kb$start<svsp.reg$V2 & cv1.ihs.mi2.10kb$end<svsp.reg$V2 | cv1.ihs.mi2.10kb$start>svsp.reg$V3 & cv1.ihs.mi2.10kb$end>svsp.reg$V3),]
co1.ihs.mi2.bk <- co1.ihs.mi2.10kb[which(co1.ihs.mi2.10kb$start<svsp.reg$V2 & co1.ihs.mi2.10kb$end<svsp.reg$V2 | co1.ihs.mi2.10kb$start>svsp.reg$V3 & co1.ihs.mi2.10kb$end>svsp.reg$V3),]
cv1.ihs.mi7.bk <- cv1.ihs.mi7.10kb[which(cv1.ihs.mi7.10kb$start<pla2.reg$V2 & cv1.ihs.mi7.10kb$end<pla2.reg$V2 | cv1.ihs.mi7.10kb$start>pla2.reg$V3 & cv1.ihs.mi7.10kb$end>pla2.reg$V3),]
co1.ihs.mi7.bk <- co1.ihs.mi7.10kb[which(co1.ihs.mi7.10kb$start<pla2.reg$V2 & co1.ihs.mi7.10kb$end<pla2.reg$V2 | co1.ihs.mi7.10kb$start>pla2.reg$V3 & co1.ihs.mi7.10kb$end>pla2.reg$V3),]

### Read in ß values------------------------------------------------------------

cv1.beta.10kb <- read.table('./beta/results/cv1.phased.all.betascores.10kb.txt',header=T)
co1.beta.10kb <- read.table('./beta/results/co1.phased.all.betascores.10kb.cnvMask.txt',header=T)
cv1.beta.1kb <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
co1.beta.1kb <- read.table('./beta/results/co1.phased.all.betascores.1kb.cnvMask.txt',header=T)

cv1.beta.mi1.10kb <- cv1.beta.10kb[which(cv1.beta.10kb$chrom=='scaffold-mi1'),]
co1.beta.mi1.10kb <- co1.beta.10kb[which(co1.beta.10kb$chrom=='scaffold-mi1'),]
cv1.beta.mi2.10kb <- cv1.beta.10kb[which(cv1.beta.10kb$chrom=='scaffold-mi2'),]
co1.beta.mi2.10kb <- co1.beta.10kb[which(co1.beta.10kb$chrom=='scaffold-mi2'),]
cv1.beta.mi7.1kb <- cv1.beta.1kb[which(cv1.beta.1kb$chrom=='scaffold-mi7'),]
co1.beta.mi7.1kb <- co1.beta.1kb[which(co1.beta.1kb$chrom=='scaffold-mi7'),]

cv1.beta.svmp <- cv1.beta.mi1.10kb[which(cv1.beta.mi1.10kb$start>=svmp.reg$V2 & cv1.beta.mi1.10kb$end<=svmp.reg$V3),]
co1.beta.svmp <- co1.beta.mi1.10kb[which(co1.beta.mi1.10kb$start>=svmp.reg$V2 & co1.beta.mi1.10kb$end<=svmp.reg$V3),]
cv1.beta.svsp <- cv1.beta.mi2.10kb[which(cv1.beta.mi2.10kb$start>=svsp.reg$V2 & cv1.beta.mi2.10kb$end<=svsp.reg$V3),]
co1.beta.svsp <- co1.beta.mi2.10kb[which(co1.beta.mi2.10kb$start>=svsp.reg$V2 & co1.beta.mi2.10kb$end<=svsp.reg$V3),]
cv1.beta.pla2 <- cv1.beta.mi7.1kb[which(cv1.beta.mi7.1kb$start>=pla2.reg$V2 & cv1.beta.mi7.1kb$end<=pla2.reg$V3),]
co1.beta.pla2 <- co1.beta.mi7.1kb[which(co1.beta.mi7.1kb$start>=pla2.reg$V2 & co1.beta.mi7.1kb$end<=pla2.reg$V3),]

cv1.beta.mi1.bk <- cv1.beta.mi1.10kb[which(cv1.beta.mi1.10kb$start<svmp.reg$V2 & cv1.beta.mi1.10kb$end<svmp.reg$V2 | cv1.beta.mi1.10kb$start>svmp.reg$V3 & cv1.beta.mi1.10kb$end>svmp.reg$V3),]
co1.beta.mi1.bk <- co1.beta.mi1.10kb[which(co1.beta.mi1.10kb$start<svmp.reg$V2 & co1.beta.mi1.10kb$end<svmp.reg$V2 | co1.beta.mi1.10kb$start>svmp.reg$V3 & co1.beta.mi1.10kb$end>svmp.reg$V3),]
cv1.beta.mi2.bk <- cv1.beta.mi2.10kb[which(cv1.beta.mi2.10kb$start<svsp.reg$V2 & cv1.beta.mi2.10kb$end<svsp.reg$V2 | cv1.beta.mi2.10kb$start>svsp.reg$V3 & cv1.beta.mi2.10kb$end>svsp.reg$V3),]
co1.beta.mi2.bk <- co1.beta.mi2.10kb[which(co1.beta.mi2.10kb$start<svsp.reg$V2 & co1.beta.mi2.10kb$end<svsp.reg$V2 | co1.beta.mi2.10kb$start>svsp.reg$V3 & co1.beta.mi2.10kb$end>svsp.reg$V3),]
cv1.beta.mi7.bk <- cv1.beta.mi7.1kb[which(cv1.beta.mi7.1kb$start<pla2.reg$V2 & cv1.beta.mi7.1kb$end<pla2.reg$V2 | cv1.beta.mi7.1kb$start>pla2.reg$V3 & cv1.beta.mi7.1kb$end>pla2.reg$V3),]
co1.beta.mi7.bk <- co1.beta.mi7.1kb[which(co1.beta.mi7.1kb$start<pla2.reg$V2 & co1.beta.mi7.1kb$end<pla2.reg$V2 | co1.beta.mi7.1kb$start>pla2.reg$V3 & co1.beta.mi7.1kb$end>pla2.reg$V3),]

### Compare chromosome-wide statistics-------------------------------------------

# ρ versus iHS across chromosomes
plot(cv1.rho.10kb.mi1$mean,cv1.ihs.mi1.10kb$iHS,pch=20)
cor.test(cv1.rho.10kb.mi1$mean,as.numeric(cv1.ihs.mi1.10kb$iHS))

plot(cv1.rho.10kb.mi2$mean,cv1.ihs.mi2.10kb$iHS,pch=20)
cor.test(cv1.rho.10kb.mi2$mean,as.numeric(cv1.ihs.mi2.10kb$iHS))

plot(cv1.rho.10kb.mi7$mean,cv1.ihs.mi7.10kb$iHS,pch=20)
cor.test(cv1.rho.10kb.mi7$mean,as.numeric(cv1.ihs.mi7.10kb$iHS))

plot(co1.rho.10kb.mi1$mean,co1.ihs.mi1.10kb$iHS,pch=20)
cor.test(co1.rho.10kb.mi1$mean,as.numeric(co1.ihs.mi1.10kb$iHS))

plot(co1.rho.10kb.mi2$mean,co1.ihs.mi2.10kb$iHS,pch=20)
cor.test(co1.rho.10kb.mi2$mean,as.numeric(co1.ihs.mi2.10kb$iHS))

plot(co1.rho.10kb.mi7$mean,co1.ihs.mi7.10kb$iHS,pch=20)
cor.test(co1.rho.10kb.mi7$mean,as.numeric(co1.ihs.mi7.10kb$iHS))

# ρ versus ß across chromosomes
plot(cv1.rho.10kb.mi1$mean,cv1.beta.mi1.10kb$Beta1.,pch=20)
cor.test(cv1.rho.10kb.mi1$mean,as.numeric(cv1.beta.mi1.10kb$Beta1.))

plot(cv1.rho.10kb.mi2$mean,cv1.beta.mi2.10kb$Beta1.,pch=20)
cor.test(cv1.rho.10kb.mi2$mean,as.numeric(cv1.beta.mi2.10kb$Beta1.))

plot(cv1.rho.1kb.mi7$mean,cv1.beta.mi7.1kb$Beta1.,pch=20)
cor.test(cv1.rho.1kb.mi7$mean,as.numeric(cv1.beta.mi7.1kb$Beta1.))

# iHS versus ß across chromosomes
plot(cv1.ihs.mi1.10kb$iHS,cv1.beta.mi1.10kb$Beta1.,pch=20)
cor.test(as.numeric(cv1.ihs.mi1.10kb$iHS),as.numeric(cv1.beta.mi1.10kb$Beta1.))

plot(cv1.ihs.mi2.10kb$iHS,cv1.beta.mi2.10kb$Beta1.,pch=20)
cor.test(as.numeric(cv1.ihs.mi2.10kb$iHS),as.numeric(cv1.beta.mi2.10kb$Beta1.))

plot(cv1.ihs.mi7.10kb$iHS,cv1.beta.mi7.10kb$Beta1.,pch=20)
cor.test(as.numeric(cv1.ihs.mi7.10kb$iHS),as.numeric(cv1.beta.mi7.10kb$Beta1.))

# ρ versus iHS in venom regions
plot(r.svmp.cv$mean_corr[1:61],cv1.ihs.svmp$iHS,pch=20)
cor.test(r.svmp.cv$mean_corr[1:61],as.numeric(cv1.ihs.svmp$iHS),method='spearman')

plot(r.svsp.cv$mean_corr[1:51],cv1.ihs.svsp$iHS,pch=20)
cor.test(r.svsp.cv$mean_corr[1:51],as.numeric(cv1.ihs.svsp$iHS),method='spearman')

plot(r.pla2.cv$mean_corr[1:43],cv1.ihs.pla2$iHS,pch=20)
cor.test(as.numeric(r.pla2.cv$mean_corr[1:43]),as.numeric(cv1.ihs.pla2$iHS),method='spearman')

plot(r.svmp.co$mean_corr[1:61],co1.ihs.svmp$iHS,pch=20)
cor.test(r.svmp.co$mean_corr[1:61],as.numeric(co1.ihs.svmp$iHS),method='spearman')

plot(r.svsp.co$mean_corr[1:51],co1.ihs.svsp$iHS,pch=20)
cor.test(r.svsp.co$mean_corr[1:51],as.numeric(co1.ihs.svsp$iHS),method='spearman')

plot(r.pla2.co$mean_corr[1:43],co1.ihs.pla2$iHS,pch=20)
cor.test(as.numeric(r.pla2.co$mean_corr[1:43]),as.numeric(co1.ihs.pla2$iHS),method='spearman')

# ρ versus ß in venom regions
plot(r.svmp.cv$mean_corr[1:61],cv1.beta.svmp$Beta1.,pch=20)
cor.test(r.svmp.cv$mean_corr[1:61],as.numeric(cv1.beta.svmp$Beta1.))

plot(r.svsp.cv$mean_corr[1:51],cv1.beta.svsp$Beta1.,pch=20)
cor.test(r.svsp.cv$mean_corr[1:51],as.numeric(cv1.beta.svsp$Beta1.))

plot(r.pla2.cv$mean_corr[1:43],cv1.beta.pla2$Beta1.,pch=20)
cor.test(r.svsp.cv$mean_corr[1:43],as.numeric(cv1.beta.pla2$Beta1.))

plot(r.svmp.co$mean_corr[1:61],co1.beta.svmp$Beta1.,pch=20)
cor.test(r.svmp.co$mean_corr[1:61],as.numeric(co1.beta.svmp$Beta1.))

plot(r.svsp.co$mean_corr[1:51],co1.beta.svsp$Beta1.,pch=20)
cor.test(r.svsp.co$mean_corr[1:51],as.numeric(co1.beta.svsp$Beta1.))

plot(r.pla2.co$mean_corr[1:43],co1.beta.pla2$Beta1.,pch=20)
cor.test(r.svsp.co$mean_corr[1:43],as.numeric(co1.beta.pla2$Beta1.))


