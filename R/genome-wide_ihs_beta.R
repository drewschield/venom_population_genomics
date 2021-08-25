############################################################################
# Genome-wide patterns of selection statistics
############################################################################

### Goal: examine genomic scans of iHS and ß in rattlesnake
### species/populations across the genome.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('./')

library(scales)
library(data.table)
library(tidyverse)
library(stringr)

### Read in data------------------------------------------------------------

cv1.ihs.1mb <- read.table('./rehh/cv.all_ihs.1Mb.txt',header=T)
cv1.ihs.100kb <- read.table('./rehh/cv.all_ihs.100kb.txt',header=T)
co1.ihs.1mb <- read.table('./rehh/co.all_ihs.1Mb.txt',header=T)
co1.ihs.100kb <- read.table('./rehh/co.all_ihs.100kb.txt',header=T)

cv1.beta.1mb <- read.table('./beta/results/cv1.phased.all.betascores.1Mb.txt',header=T)
cv1.beta.100kb <- read.table('./beta/results/cv1.phased.all.betascores.100kb.txt',header=T)
co1.beta.1mb <- read.table('./beta/results/co1.phased.all.betascores.1Mb.txt',header=T)
co1.beta.100kb <- read.table('./beta/results/co1.phased.all.betascores.100kb.txt',header=T)

### Plot scans---------------------------------------------------------------

color1 = rep(NA, length=length(cv1.ihs.100kb$chrom))
color1[which(cv1.ihs.100kb$chrom=="scaffold-ma1")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-ma2")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-ma3")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-ma4")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-ma5")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-ma6")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-ma7")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-Z")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi1")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi2")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi3")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi4")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi5")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi6")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi7")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi8")] = "goldenrod"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi9")] = "darkorange"
color1[which(cv1.ihs.100kb$chrom=="scaffold-mi10")] = "goldenrod"

# Multiply 1 Mb rows by 9.959 because 1303 windows * 10 = 13030; 12977 100 kb windows exist, so 12977/13030 = 0.9959
# Export as 10 x 8
par(mfrow=c(4,1))
plot(rownames(cv1.ihs.100kb),cv1.ihs.100kb$iHS,pch=20,cex=0.9,col=alpha(color1,0.05),ylab='iHS',xlab='Chromosome Position (Mb)')
lines(as.integer(rownames(cv1.ihs.1mb))*9.959,cv1.ihs.1mb$iHS,lwd=0.85)
plot(rownames(co1.ihs.100kb),co1.ihs.100kb$iHS,pch=20,cex=0.9,col=alpha(color1,0.05),ylab='iHS',xlab='Chromosome Position (Mb)')
lines(as.integer(rownames(co1.ihs.1mb))*9.959,co1.ihs.1mb$iHS,lwd=0.85)
plot(rownames(cv1.beta.100kb),as.numeric(cv1.beta.100kb$Beta1.),pch=20,cex=0.9,col=alpha(color1,0.05),ylim=c(-1,5),ylab='ß',xlab='Chromosome Position (Mb)')
lines(as.integer(rownames(cv1.beta.1mb))*9.959,as.numeric(cv1.beta.1mb$Beta1.),lwd=0.85)
plot(rownames(co1.beta.100kb),co1.beta.100kb$Beta1.,pch=20,cex=0.9,col=alpha(color1,0.05),ylim=c(-1,5),ylab='ß',xlab='Chromosome Position (Mb)')
lines(as.integer(rownames(co1.beta.1mb))*9.959,co1.beta.1mb$Beta1.,lwd=0.85)


