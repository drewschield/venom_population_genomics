############################################################################
# Trans-species polymorphism in venom regions
############################################################################

### Goal: quantify trans-species polymorphisms inside/outside venom gene
### regions.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & load dependencies-------------------------------

setwd('./trans-species_polymorphism/')

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

poly.mi1 <- read.table('./mi1_MajorAlleleFreqs_27July21.txt',header=T)
poly.mi2 <- read.table('./mi2_MajorAlleleFreqs_27July21.txt',header=T)
poly.mi7 <- read.table('./mi7_MajorAlleleFreqs_27July21.txt',header=T)

### Parse venom region & background-----------------------------------------

poly.svmp <- poly.mi1[which(poly.mi1$Location>=svmp.reg$V2 & poly.mi1$Location<=svmp.reg$V3),]
poly.svsp <- poly.mi2[which(poly.mi2$Location>=svsp.reg$V2 & poly.mi2$Location<=svsp.reg$V3),]
poly.pla2 <- poly.mi7[which(poly.mi7$Location>=pla2.reg$V2 & poly.mi7$Location<=pla2.reg$V3),]

poly.mi1.bk <- poly.mi1[which(poly.mi1$Location<svmp.reg$V2 | poly.mi1$Location>svmp.reg$V3),]
poly.mi2.bk <- poly.mi2[which(poly.mi2$Location<svsp.reg$V2 | poly.mi2$Location>svsp.reg$V3),]
poly.mi7.bk <- poly.mi7[which(poly.mi7$Location<pla2.reg$V2 | poly.mi7$Location>pla2.reg$V3),]

# If we want to focus on PLA2 A1, which has no evidence of CNV between species
poly.pla2 <- poly.mi7[which(poly.mi7$Location>=pla2$V2[5] & poly.mi7$Location<=pla2$V3[5]),]
poly.mi7.bk <- poly.mi7[which(poly.mi7$Location<pla2$V2[5] | poly.mi7$Location>pla2$V3[5]),]

### Quantify trans-species polymorphisms------------------------------------

# Count fixed differences
count(poly.svmp$Pop1==1 & poly.svmp$Pop2==1)
count(poly.mi1.bk$Pop1==1 & poly.mi1.bk$Pop2==1)

count(poly.svsp$Pop1==1 & poly.svsp$Pop2==1)
count(poly.mi2.bk$Pop1==1 & poly.mi2.bk$Pop2==1)

count(poly.pla2$Pop1==1 & poly.pla2$Pop2==1)
count(poly.mi7.bk$Pop1==1 & poly.mi7.bk$Pop2==1)

# Count trans-species polymorphisms in venom and background regions
count((1-poly.svmp$Pop1)>0 & 0.5 >=(1-poly.svmp$Pop1) & (1-poly.svmp$Pop2)>0 & 0.5>=(1-poly.svmp$Pop2))
count((1-poly.mi1.bk$Pop1)>0 & 0.5 >=(1-poly.mi1.bk$Pop1) & (1-poly.mi1.bk$Pop2)>0 & 0.5>=(1-poly.mi1.bk$Pop2))

count((1-poly.svsp$Pop1)>0 & 0.5 >=(1-poly.svsp$Pop1) & (1-poly.svsp$Pop2)>0 & 0.5>=(1-poly.svsp$Pop2))
count((1-poly.mi2.bk$Pop1)>0 & 0.5 >=(1-poly.mi2.bk$Pop1) & (1-poly.mi2.bk$Pop2)>0 & 0.5>=(1-poly.mi2.bk$Pop2))

count((1-poly.pla2$Pop1)>0 & 0.5 >=(1-poly.pla2$Pop1) & (1-poly.pla2$Pop2)>0 & 0.5>=(1-poly.pla2$Pop2))
count((1-poly.mi7.bk$Pop1)>0 & 0.5 >=(1-poly.mi7.bk$Pop1) & (1-poly.mi7.bk$Pop2)>0 & 0.5>=(1-poly.mi7.bk$Pop2))

# Count private polymorphisms per species
count((1-poly.svmp$Pop2)==0 & (1-poly.svmp$Pop1>0))
count((1-poly.svmp$Pop1)==0 & (1-poly.svmp$Pop2>0))

count((1-poly.mi1.bk$Pop2)==0 & (1-poly.mi1.bk$Pop1>0))
count((1-poly.mi1.bk$Pop1)==0 & (1-poly.mi1.bk$Pop2>0))

count((1-poly.svsp$Pop2)==0 & (1-poly.svsp$Pop1>0))
count((1-poly.svsp$Pop1)==0 & (1-poly.svsp$Pop2>0))

count((1-poly.mi2.bk$Pop2)==0 & (1-poly.mi2.bk$Pop1>0))
count((1-poly.mi2.bk$Pop1)==0 & (1-poly.mi2.bk$Pop2>0))

count((1-poly.pla2$Pop2)==0 & (1-poly.pla2$Pop1>0))
count((1-poly.pla2$Pop1)==0 & (1-poly.pla2$Pop2>0))

count((1-poly.mi7.bk$Pop2)==0 & (1-poly.mi7.bk$Pop1>0))
count((1-poly.mi7.bk$Pop1)==0 & (1-poly.mi7.bk$Pop2>0))

### Pie charts!-------------------------------------------------------------

par(mfrow=c(3,2))
# Simple Pie Chart
slices <- c(280,1203,4058,4649)
lbls <- c("Fixed","TSP", "CV1", "CO1")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
pie(slices, labels = lbls, main="SVMP",radius=1)

slices.bk <- c(19677,18478,120921,100755)
lbls <- c("Fixed","TSP", "CV1", "CO1")
pct <- round(slices.bk/sum(slices.bk)*100)
lbls <- paste(lbls, pct) # add percents to labels
pie(slices.bk, labels = lbls, main="Chromosome 9",radius=1)

slices <- c(242,781,3084,2568)
lbls <- c("Fixed","TSP", "CV1", "CO1")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
pie(slices, labels = lbls, main="SVSP",radius=1)

slices.bk <- c(11399,21559,96246,72207)
lbls <- c("Fixed","TSP", "CV1", "CO1")
pct <- round(slices.bk/sum(slices.bk)*100)
lbls <- paste(lbls, pct) # add percents to labels
pie(slices.bk, labels = lbls, main="Chromosome 10",radius=1)

slices <- c(4,20,47,42)
lbls <- c("Fixed","TSP", "CV1", "CO1")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
pie(slices, labels = lbls, main="PLA2",radius=1)

slices.bk <- c(9647,13787,79692,60870)
lbls <- c("Fixed","TSP", "CV1", "CO1")
pct <- round(slices.bk/sum(slices.bk)*100)
lbls <- paste(lbls, pct) # add percents to labels
pie(slices.bk, labels = lbls, main="Chromosome 15",radius=1)

### Fisher's exact tests----------------------------------------------------

# Set up matrices
svmp.tsp.comp <- matrix(c(1203,10190,18478,259831),nrow = 2)
svsp.tsp.comp <- matrix(c(781,6675,21559,201411),nrow = 2)
pla2.tsp.comp <- matrix(c(20,113,13857,164702),nrow = 2)

# Perform Fisher's exact tests
fisher.test(svmp.tsp.comp)
fisher.test(svsp.tsp.comp)
fisher.test(pla2.tsp.comp)

### Read in amino acid results----------------------------------------------

aa.svmp <- read.table('./amino_acid/SVMP_all.txt')
aa.svsp <- read.table('./amino_acid/SVSP_all.txt')
aa.pla2 <- read.table('./amino_acid/PLA2_all.txt')

### Quantify trans-species polymorphisms------------------------------------

# Count trans-species polymorphisms in venom and background regions
count((1-aa.svmp$V2)>0.1 & 0.5>=(1-aa.svmp$V2) & (1-aa.svmp$V3)>0.1 & 0.5>=(1-aa.svmp$V3))
count((1-aa.svsp$V2)>0.1 & 0.5>=(1-aa.svsp$V2) & (1-aa.svsp$V3)>0.1 & 0.5>=(1-aa.svsp$V3))
count((1-aa.pla2$V2)>0.1 & 0.5>=(1-aa.pla2$V2) & (1-aa.pla2$V3)>0.1 & 0.5>=(1-aa.pla2$V3))

# Count private polymorphisms per species
count((1-aa.svmp$V3)==0 & (1-aa.svmp$V2>0.1))
count((1-aa.svmp$V2)==0 & (1-aa.svmp$V3>0.1))

count((1-aa.svsp$V3)==0 & (1-aa.svsp$V2>0.1))
count((1-aa.svsp$V2)==0 & (1-aa.svsp$V3>0.1))

count((1-aa.pla2$V3)==0 & (1-aa.pla2$V2>0.1))
count((1-aa.pla2$V2)==0 & (1-aa.pla2$V3>0.1))
