############################################################################
# Compare π, dxy, and Fst values for venom genes to intergenic regions
############################################################################

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory and load dependencies-----------------------------

setwd('./gene_vs_intergenic/')

library(scales)
library(data.table)
library(dplyr)

### Read in data------------------------------------------------------------

svmp.pi.ig.cv1 <- read.table('svmp_intergenic.pi_mean.cv1.txt',header=T)
svmp.pi.ig.cv2 <- read.table('svmp_intergenic.pi_mean.cv2.txt',header=T)
svmp.pi.ig.co1 <- read.table('svmp_intergenic.pi_mean.co1.txt',header=T)
svmp.pi.ig.co2 <- read.table('svmp_intergenic.pi_mean.co2.txt',header=T)
svmp.pi.g.cv1 <- read.table('svmp_genic.pi_mean.cv1.txt',header=T)
svmp.pi.g.cv2 <- read.table('svmp_genic.pi_mean.cv2.txt',header=T)
svmp.pi.g.co1 <- read.table('svmp_genic.pi_mean.co1.txt',header=T)
svmp.pi.g.co2 <- read.table('svmp_genic.pi_mean.co2.txt',header=T)

svmp.dxy.ig.cv1co1 <- read.table('svmp_intergenic.dxy_mean.cv1co1.txt',header=T)
svmp.dxy.ig.cv1cv2 <- read.table('svmp_intergenic.dxy_mean.cv1cv2.txt',header=T)
svmp.dxy.ig.co1co2 <- read.table('svmp_intergenic.dxy_mean.co1co2.txt',header=T)
svmp.dxy.g.cv1co1 <- read.table('svmp_genic.dxy_mean.cv1co1.txt',header=T)
svmp.dxy.g.cv1cv2 <- read.table('svmp_genic.dxy_mean.cv1cv2.txt',header=T)
svmp.dxy.g.co1co2 <- read.table('svmp_genic.dxy_mean.co1co2.txt',header=T)

svmp.fst.ig.cv1co1 <- read.table('svmp_intergenic.fst_mean.cv1co1.txt',header=T)
svmp.fst.ig.cv1cv2 <- read.table('svmp_intergenic.fst_mean.cv1cv2.txt',header=T)
svmp.fst.ig.co1co2 <- read.table('svmp_intergenic.fst_mean.co1co2.txt',header=T)
svmp.fst.g.cv1co1 <- read.table('svmp_genic.fst_mean.cv1co1.txt',header=T)
svmp.fst.g.cv1cv2 <- read.table('svmp_genic.fst_mean.cv1cv2.txt',header=T)
svmp.fst.g.co1co2 <- read.table('svmp_genic.fst_mean.co1co2.txt',header=T)

svsp.pi.ig.cv1 <- read.table('svsp_intergenic.pi_mean.cv1.txt',header=T)
svsp.pi.ig.cv2 <- read.table('svsp_intergenic.pi_mean.cv2.txt',header=T)
svsp.pi.ig.co1 <- read.table('svsp_intergenic.pi_mean.co1.txt',header=T)
svsp.pi.ig.co2 <- read.table('svsp_intergenic.pi_mean.co2.txt',header=T)
svsp.pi.g.cv1 <- read.table('svsp_genic.pi_mean.cv1.txt',header=T)
svsp.pi.g.cv2 <- read.table('svsp_genic.pi_mean.cv2.txt',header=T)
svsp.pi.g.co1 <- read.table('svsp_genic.pi_mean.co1.txt',header=T)
svsp.pi.g.co2 <- read.table('svsp_genic.pi_mean.co2.txt',header=T)

svsp.dxy.ig.cv1co1 <- read.table('svsp_intergenic.dxy_mean.cv1co1.txt',header=T)
svsp.dxy.ig.cv1cv2 <- read.table('svsp_intergenic.dxy_mean.cv1cv2.txt',header=T)
svsp.dxy.ig.co1co2 <- read.table('svsp_intergenic.dxy_mean.co1co2.txt',header=T)
svsp.dxy.g.cv1co1 <- read.table('svsp_genic.dxy_mean.cv1co1.txt',header=T)
svsp.dxy.g.cv1cv2 <- read.table('svsp_genic.dxy_mean.cv1cv2.txt',header=T)
svsp.dxy.g.co1co2 <- read.table('svsp_genic.dxy_mean.co1co2.txt',header=T)

svsp.fst.ig.cv1co1 <- read.table('svsp_intergenic.fst_mean.cv1co1.txt',header=T)
svsp.fst.ig.cv1cv2 <- read.table('svsp_intergenic.fst_mean.cv1cv2.txt',header=T)
svsp.fst.ig.co1co2 <- read.table('svsp_intergenic.fst_mean.co1co2.txt',header=T)
svsp.fst.g.cv1co1 <- read.table('svsp_genic.fst_mean.cv1co1.txt',header=T)
svsp.fst.g.cv1cv2 <- read.table('svsp_genic.fst_mean.cv1cv2.txt',header=T)
svsp.fst.g.co1co2 <- read.table('svsp_genic.fst_mean.co1co2.txt',header=T)

pla2.pi.ig.cv1 <- read.table('pla2_intergenic.pi_mean.cv1.txt',header=T)
pla2.pi.ig.cv2 <- read.table('pla2_intergenic.pi_mean.cv2.txt',header=T)
pla2.pi.ig.co1 <- read.table('pla2_intergenic.pi_mean.co1.txt',header=T)
pla2.pi.ig.co2 <- read.table('pla2_intergenic.pi_mean.co2.txt',header=T)
pla2.pi.g.cv1 <- read.table('pla2_genic.pi_mean.cv1.txt',header=T)
pla2.pi.g.cv2 <- read.table('pla2_genic.pi_mean.cv2.txt',header=T)
pla2.pi.g.co1 <- read.table('pla2_genic.pi_mean.co1.txt',header=T)
pla2.pi.g.co2 <- read.table('pla2_genic.pi_mean.co2.txt',header=T)

pla2.dxy.ig.cv1co1 <- read.table('pla2_intergenic.dxy_mean.cv1co1.txt',header=T)
pla2.dxy.ig.cv1cv2 <- read.table('pla2_intergenic.dxy_mean.cv1cv2.txt',header=T)
pla2.dxy.ig.co1co2 <- read.table('pla2_intergenic.dxy_mean.co1co2.txt',header=T)
pla2.dxy.g.cv1co1 <- read.table('pla2_genic.dxy_mean.cv1co1.txt',header=T)
pla2.dxy.g.cv1cv2 <- read.table('pla2_genic.dxy_mean.cv1cv2.txt',header=T)
pla2.dxy.g.co1co2 <- read.table('pla2_genic.dxy_mean.co1co2.txt',header=T)

pla2.fst.ig.cv1co1 <- read.table('pla2_intergenic.fst_mean.cv1co1.txt',header=T)
pla2.fst.ig.cv1cv2 <- read.table('pla2_intergenic.fst_mean.cv1cv2.txt',header=T)
pla2.fst.ig.co1co2 <- read.table('pla2_intergenic.fst_mean.co1co2.txt',header=T)
pla2.fst.g.cv1co1 <- read.table('pla2_genic.fst_mean.cv1co1.txt',header=T)
pla2.fst.g.cv1cv2 <- read.table('pla2_genic.fst_mean.cv1cv2.txt',header=T)
pla2.fst.g.co1co2 <- read.table('pla2_genic.fst_mean.co1co2.txt',header=T)

### Compare genic vs intergenic distributions-------------------------------

boxplot(svmp.pi.g.cv1$pi_mean,svmp.pi.ig.cv1$pi_mean,svmp.pi.g.cv2$pi_mean,svmp.pi.ig.cv2$pi_mean,svmp.pi.g.co1$pi_mean,svmp.pi.ig.co1$pi_mean,svmp.pi.g.co2$pi_mean,svmp.pi.ig.co2$pi_mean,outline=F)
boxplot(svmp.dxy.g.cv1co1$dxy_mean,svmp.dxy.ig.cv1co1$dxy_mean,svmp.dxy.g.cv1cv2$dxy_mean,svmp.dxy.ig.cv1cv2$dxy_mean,svmp.dxy.g.co1co2$dxy_mean,svmp.dxy.ig.co1co2$dxy_mean,outline=F)
boxplot(svmp.fst.g.cv1co1$fst_mean,svmp.fst.ig.cv1co1$fst_mean,svmp.fst.g.cv1cv2$fst_mean,svmp.fst.ig.cv1cv2$fst_mean,svmp.fst.g.co1co2$fst_mean,svmp.fst.ig.co1co2$fst_mean,outline=F)

boxplot(svsp.pi.g.cv1$pi_mean,svsp.pi.ig.cv1$pi_mean,svsp.pi.g.cv2$pi_mean,svsp.pi.ig.cv2$pi_mean,svsp.pi.g.co1$pi_mean,svsp.pi.ig.co1$pi_mean,svsp.pi.g.co2$pi_mean,svsp.pi.ig.co2$pi_mean,outline=F)
boxplot(svsp.dxy.g.cv1co1$dxy_mean,svsp.dxy.ig.cv1co1$dxy_mean,svsp.dxy.g.cv1cv2$dxy_mean,svsp.dxy.ig.cv1cv2$dxy_mean,svsp.dxy.g.co1co2$dxy_mean,svsp.dxy.ig.co1co2$dxy_mean,outline=F)
boxplot(svsp.fst.g.cv1co1$fst_mean,svsp.fst.ig.cv1co1$fst_mean,svsp.fst.g.cv1cv2$fst_mean,svsp.fst.ig.cv1cv2$fst_mean,svsp.fst.g.co1co2$fst_mean,svsp.fst.ig.co1co2$fst_mean,outline=F)

boxplot(pla2.pi.g.cv1$pi_mean,pla2.pi.ig.cv1$pi_mean,pla2.pi.g.cv2$pi_mean,pla2.pi.ig.cv2$pi_mean,pla2.pi.g.co1$pi_mean,pla2.pi.ig.co1$pi_mean,pla2.pi.g.co2$pi_mean,pla2.pi.ig.co2$pi_mean,outline=F)
boxplot(pla2.dxy.g.cv1co1$dxy_mean,pla2.dxy.ig.cv1co1$dxy_mean,pla2.dxy.g.cv1cv2$dxy_mean,pla2.dxy.ig.cv1cv2$dxy_mean,pla2.dxy.g.co1co2$dxy_mean,pla2.dxy.ig.co1co2$dxy_mean,outline=F)
boxplot(pla2.fst.g.cv1co1$fst_mean,pla2.fst.ig.cv1co1$fst_mean,pla2.fst.g.cv1cv2$fst_mean,pla2.fst.ig.cv1cv2$fst_mean,pla2.fst.g.co1co2$fst_mean,pla2.fst.ig.co1co2$fst_mean,outline=F)

### Calculate summary statistics--------------------------------------------

mean(svmp.pi.g.cv1$pi_mean,na.rm=T)
mean(svmp.pi.g.cv2$pi_mean,na.rm=T)
mean(svmp.pi.g.co1$pi_mean,na.rm=T)
mean(svmp.pi.g.co2$pi_mean,na.rm=T)
mean(svmp.dxy.g.cv1co1$dxy_mean,na.rm=T)
mean(svmp.dxy.g.cv1cv2$dxy_mean,na.rm=T)
mean(svmp.dxy.g.co1co2$dxy_mean,na.rm=T)

mean(svsp.pi.g.cv1$pi_mean,na.rm=T)
mean(svsp.pi.g.cv2$pi_mean,na.rm=T)
mean(svsp.pi.g.co1$pi_mean,na.rm=T)
mean(svsp.pi.g.co2$pi_mean,na.rm=T)
mean(svsp.dxy.g.cv1co1$dxy_mean,na.rm=T)
mean(svsp.dxy.g.cv1cv2$dxy_mean,na.rm=T)
mean(svsp.dxy.g.co1co2$dxy_mean,na.rm=T)

mean(pla2.pi.g.cv1$pi_mean,na.rm=T)
mean(pla2.pi.g.cv2$pi_mean,na.rm=T)
mean(pla2.pi.g.co1$pi_mean,na.rm=T)
mean(pla2.pi.g.co2$pi_mean,na.rm=T)
mean(pla2.dxy.g.cv1co1$dxy_mean,na.rm=T)
mean(pla2.dxy.g.cv1cv2$dxy_mean,na.rm=T)
mean(pla2.dxy.g.co1co2$dxy_mean,na.rm=T)

sd(svmp.pi.g.cv1$pi_mean,na.rm=T)
sd(svmp.pi.g.cv2$pi_mean,na.rm=T)
sd(svmp.pi.g.co1$pi_mean,na.rm=T)
sd(svmp.pi.g.co2$pi_mean,na.rm=T)
sd(svmp.dxy.g.cv1co1$dxy_mean,na.rm=T)
sd(svmp.dxy.g.cv1cv2$dxy_mean,na.rm=T)
sd(svmp.dxy.g.co1co2$dxy_mean,na.rm=T)

sd(svsp.pi.g.cv1$pi_mean,na.rm=T)
sd(svsp.pi.g.cv2$pi_mean,na.rm=T)
sd(svsp.pi.g.co1$pi_mean,na.rm=T)
sd(svsp.pi.g.co2$pi_mean,na.rm=T)
sd(svsp.dxy.g.cv1co1$dxy_mean,na.rm=T)
sd(svsp.dxy.g.cv1cv2$dxy_mean,na.rm=T)
sd(svsp.dxy.g.co1co2$dxy_mean,na.rm=T)

sd(pla2.pi.g.cv1$pi_mean,na.rm=T)
sd(pla2.pi.g.cv2$pi_mean,na.rm=T)
sd(pla2.pi.g.co1$pi_mean,na.rm=T)
sd(pla2.pi.g.co2$pi_mean,na.rm=T)
sd(pla2.dxy.g.cv1co1$dxy_mean,na.rm=T)
sd(pla2.dxy.g.cv1cv2$dxy_mean,na.rm=T)
sd(pla2.dxy.g.co1co2$dxy_mean,na.rm=T)

mean(svmp.pi.ig.cv1$pi_mean,na.rm=T)
mean(svmp.pi.ig.cv2$pi_mean,na.rm=T)
mean(svmp.pi.ig.co1$pi_mean,na.rm=T)
mean(svmp.pi.ig.co2$pi_mean,na.rm=T)
mean(svmp.dxy.ig.cv1co1$dxy_mean,na.rm=T)
mean(svmp.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
mean(svmp.dxy.ig.co1co2$dxy_mean,na.rm=T)

mean(svsp.pi.ig.cv1$pi_mean,na.rm=T)
mean(svsp.pi.ig.cv2$pi_mean,na.rm=T)
mean(svsp.pi.ig.co1$pi_mean,na.rm=T)
mean(svsp.pi.ig.co2$pi_mean,na.rm=T)
mean(svsp.dxy.ig.cv1co1$dxy_mean,na.rm=T)
mean(svsp.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
mean(svsp.dxy.ig.co1co2$dxy_mean,na.rm=T)

mean(pla2.pi.ig.cv1$pi_mean,na.rm=T)
mean(pla2.pi.ig.cv2$pi_mean,na.rm=T)
mean(pla2.pi.ig.co1$pi_mean,na.rm=T)
mean(pla2.pi.ig.co2$pi_mean,na.rm=T)
mean(pla2.dxy.ig.cv1co1$dxy_mean,na.rm=T)
mean(pla2.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
mean(pla2.dxy.ig.co1co2$dxy_mean,na.rm=T)

sd(svmp.pi.ig.cv1$pi_mean,na.rm=T)
sd(svmp.pi.ig.cv2$pi_mean,na.rm=T)
sd(svmp.pi.ig.co1$pi_mean,na.rm=T)
sd(svmp.pi.ig.co2$pi_mean,na.rm=T)
sd(svmp.dxy.ig.cv1co1$dxy_mean,na.rm=T)
sd(svmp.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
sd(svmp.dxy.ig.co1co2$dxy_mean,na.rm=T)

sd(svsp.pi.ig.cv1$pi_mean,na.rm=T)
sd(svsp.pi.ig.cv2$pi_mean,na.rm=T)
sd(svsp.pi.ig.co1$pi_mean,na.rm=T)
sd(svsp.pi.ig.co2$pi_mean,na.rm=T)
sd(svsp.dxy.ig.cv1co1$dxy_mean,na.rm=T)
sd(svsp.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
sd(svsp.dxy.ig.co1co2$dxy_mean,na.rm=T)

sd(pla2.pi.ig.cv1$pi_mean,na.rm=T)
sd(pla2.pi.ig.cv2$pi_mean,na.rm=T)
sd(pla2.pi.ig.co1$pi_mean,na.rm=T)
sd(pla2.pi.ig.co2$pi_mean,na.rm=T)
sd(pla2.dxy.ig.cv1co1$dxy_mean,na.rm=T)
sd(pla2.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
sd(pla2.dxy.ig.co1co2$dxy_mean,na.rm=T)

### Mann-Whitney U statistics-----------------------------------------------

wilcox.test(svmp.pi.g.cv1$pi_mean,svmp.pi.ig.cv1$pi_mean)
wilcox.test(svmp.pi.g.cv2$pi_mean,svmp.pi.ig.cv2$pi_mean)
wilcox.test(svmp.pi.g.co1$pi_mean,svmp.pi.ig.co1$pi_mean)
wilcox.test(svmp.pi.g.co2$pi_mean,svmp.pi.ig.co2$pi_mean)

wilcox.test(svmp.dxy.g.cv1co1$dxy_mean,svmp.dxy.ig.cv1co1$dxy_mean)
wilcox.test(svmp.dxy.g.cv1cv2$dxy_mean,svmp.dxy.ig.cv1cv2$dxy_mean)
wilcox.test(svmp.dxy.g.co1co2$dxy_mean,svmp.dxy.ig.co1co2$dxy_mean)

wilcox.test(svsp.pi.g.cv1$pi_mean,svsp.pi.ig.cv1$pi_mean)
wilcox.test(svsp.pi.g.cv2$pi_mean,svsp.pi.ig.cv2$pi_mean)
wilcox.test(svsp.pi.g.co1$pi_mean,svsp.pi.ig.co1$pi_mean)
wilcox.test(svsp.pi.g.co2$pi_mean,svsp.pi.ig.co2$pi_mean)

wilcox.test(svsp.dxy.g.cv1co1$dxy_mean,svsp.dxy.ig.cv1co1$dxy_mean)
wilcox.test(svsp.dxy.g.cv1cv2$dxy_mean,svsp.dxy.ig.cv1cv2$dxy_mean)
wilcox.test(svsp.dxy.g.co1co2$dxy_mean,svsp.dxy.ig.co1co2$dxy_mean)

wilcox.test(pla2.pi.g.cv1$pi_mean,pla2.pi.ig.cv1$pi_mean)
wilcox.test(pla2.pi.g.cv2$pi_mean,pla2.pi.ig.cv2$pi_mean)
wilcox.test(pla2.pi.g.co1$pi_mean,pla2.pi.ig.co1$pi_mean)
wilcox.test(pla2.pi.g.co2$pi_mean,pla2.pi.ig.co2$pi_mean)

wilcox.test(pla2.dxy.g.cv1co1$dxy_mean,pla2.dxy.ig.cv1co1$dxy_mean)
wilcox.test(pla2.dxy.g.cv1cv2$dxy_mean,pla2.dxy.ig.cv1cv2$dxy_mean)
wilcox.test(pla2.dxy.g.co1co2$dxy_mean,pla2.dxy.ig.co1co2$dxy_mean)


mean(svmp.dxy.g.cv1co1$dxy_mean,na.rm=T)/mean(svmp.dxy.ig.cv1co1$dxy_mean,na.rm=T)
mean(svmp.dxy.g.cv1cv2$dxy_mean,na.rm=T)/mean(svmp.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
mean(svmp.dxy.g.co1co2$dxy_mean,na.rm=T)/mean(svmp.dxy.ig.co1co2$dxy_mean,na.rm=T)

mean(svsp.dxy.g.cv1co1$dxy_mean,na.rm=T)/mean(svsp.dxy.ig.cv1co1$dxy_mean,na.rm=T)
median(svsp.dxy.g.cv1cv2$dxy_mean,na.rm=T)/median(svsp.dxy.ig.cv1cv2$dxy_mean,na.rm=T)
mean(svsp.dxy.g.co1co2$dxy_mean,na.rm=T)/mean(svsp.dxy.ig.co1co2$dxy_mean,na.rm=T)
