############################################################################
# Genome-wide patterns of genetic diversity and divergence
############################################################################

### Goal: examine genomic scans of pi, dxy, and fst between rattlesnake
### species/populations across the genome and make comparisons to characterize
### the genomic landscape of divergence.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('./pixy/')

library(scales)
library(data.table)
library(tidyverse)
library(stringr)

### Read in data------------------------------------------------------------

pi.100kb <- read.table('pixy.all.100kb_pi.txt',header=T)
dxy.100kb <- read.table('pixy.all.100kb_dxy.txt',header=T)
fst.100kb <- read.table('pixy.all.100kb_fst.txt',header=T)

pi.10kb <- read.table('pixy.all.10kb_pi.txt',header=T)
dxy.10kb <- read.table('pixy.all.10kb_dxy.txt',header=T)
fst.10kb <- read.table('pixy.all.10kb_fst.txt',header=T)

rho.100kb.cv1 <- read.table('../recombination/sliding_windows/viridis.recomb.bpen10.windowed.100kb.centromereMask.txt',header=T)
rho.100kb.co1 <- read.table('../recombination/sliding_windows/oreganus.recomb.bpen10.windowed.100kb.centromereMask.txt',header=T)

### Parse data---------------------------------------------------------------

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

### Write 100 kb data to files-----------------------------------------------

write.table(pi.100kb.cv1,'./pi.100kb.cv1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(pi.100kb.cv2,'./pi.100kb.cv2.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(pi.100kb.co1,'./pi.100kb.co1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(pi.100kb.co2,'./pi.100kb.co2.txt',quote=FALSE,row.names=FALSE,sep = "\t")

write.table(dxy.100kb.cv1co1,'./dxy.100kb.cv1co1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(dxy.100kb.cv1cv2,'./dxy.100kb.cv1cv2.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(dxy.100kb.co1co2,'./dxy.100kb.co1co2.txt',quote=FALSE,row.names=FALSE,sep = "\t")

write.table(fst.100kb.cv1co1,'./fst.100kb.cv1co1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(fst.100kb.cv1cv2,'./fst.100kb.cv1cv2.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(fst.100kb.co1co2,'./fst.100kb.co1co2.txt',quote=FALSE,row.names=FALSE,sep = "\t")

### Write 10 kb data to files-----------------------------------------------

write.table(pi.10kb.cv1,'./pi.10kb.cv1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(pi.10kb.cv2,'./pi.10kb.cv2.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(pi.10kb.co1,'./pi.10kb.co1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(pi.10kb.co2,'./pi.10kb.co2.txt',quote=FALSE,row.names=FALSE,sep = "\t")

write.table(dxy.10kb.cv1co1,'./dxy.10kb.cv1co1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(dxy.10kb.cv1cv2,'./dxy.10kb.cv1cv2.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(dxy.10kb.co1co2,'./dxy.10kb.co1co2.txt',quote=FALSE,row.names=FALSE,sep = "\t")

write.table(fst.10kb.cv1co1,'./fst.10kb.cv1co1.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(fst.10kb.cv1cv2,'./fst.10kb.cv1cv2.txt',quote=FALSE,row.names=FALSE,sep = "\t")
write.table(fst.10kb.co1co2,'./fst.10kb.co1co2.txt',quote=FALSE,row.names=FALSE,sep = "\t")

### Read in 1 Mb windowed averages-------------------------------------------

pi.1mb.cv1 <- read.table('pi.1Mb.cv1.txt',header=T)
pi.1mb.cv2 <- read.table('pi.1Mb.cv2.txt',header=T)
pi.1mb.co1 <- read.table('pi.1Mb.co1.txt',header=T)
pi.1mb.co2 <- read.table('pi.1Mb.co2.txt',header=T)

dxy.1mb.cv1co1 <- read.table('dxy.1Mb.cv1co1.txt',header=T)
dxy.1mb.cv1cv2 <- read.table('dxy.1Mb.cv1cv2.txt',header=T)
dxy.1mb.co1co2 <- read.table('dxy.1Mb.co1co2.txt',header=T)

fst.1mb.cv1co1 <- read.table('fst.1Mb.cv1co1.txt',header=T)
fst.1mb.cv1cv2 <- read.table('fst.1Mb.cv1cv2.txt',header=T)
fst.1mb.co1co2 <- read.table('fst.1Mb.co1co2.txt',header=T)

### Plot scans---------------------------------------------------------------

color = rep(NA, length=length(pi.100kb.cv1$chromosome))
color[which(pi.10kb.cv1$chromosome=="scaffold-ma1")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-ma2")] = "grey"
color[which(pi.10kb.cv1$chromosome=="scaffold-ma3")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-ma4")] = "grey"
color[which(pi.10kb.cv1$chromosome=="scaffold-ma5")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-ma6")] = "grey"
color[which(pi.10kb.cv1$chromosome=="scaffold-ma7")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-Z")] = "seagreen"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi1")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi2")] = "grey"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi3")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi4")] = "grey"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi5")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi6")] = "grey"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi7")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi8")] = "grey"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi9")] = "black"
color[which(pi.10kb.cv1$chromosome=="scaffold-mi10")] = "grey"

color1 = rep(NA, length=length(pi.100kb.cv1$chromosome))
color1[which(pi.100kb.cv1$chromosome=="scaffold-ma1")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-ma2")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-ma3")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-ma4")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-ma5")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-ma6")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-ma7")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-Z")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi1")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi2")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi3")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi4")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi5")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi6")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi7")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi8")] = "goldenrod"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi9")] = "darkorange"
color1[which(pi.100kb.cv1$chromosome=="scaffold-mi10")] = "goldenrod"

color2 = rep(NA, length=length(pi.100kb.cv1$chromosome))
color2[which(pi.100kb.cv1$chromosome=="scaffold-ma1")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-ma2")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-ma3")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-ma4")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-ma5")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-ma6")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-ma7")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-Z")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi1")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi2")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi3")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi4")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi5")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi6")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi7")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi8")] = "lightblue"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi9")] = "navy"
color2[which(pi.100kb.cv1$chromosome=="scaffold-mi10")] = "lightblue"

color3 = rep(NA, length=length(pi.100kb.cv1$chromosome))
color3[which(pi.100kb.cv1$chromosome=="scaffold-ma1")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-ma2")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-ma3")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-ma4")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-ma5")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-ma6")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-ma7")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-Z")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi1")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi2")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi3")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi4")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi5")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi6")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi7")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi8")] = "purple4"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi9")] = "purple"
color3[which(pi.100kb.cv1$chromosome=="scaffold-mi10")] = "purple4"

# Multiply 1 Mb rows by 9.959 because 1303 windows * 10 = 13030; 12977 100 kb windows exist, so 12977/13030 = 0.9959
# Export as 10 x 8
par(mfrow=c(4,1))
plot(rownames(pi.100kb.cv1),pi.100kb.cv1$avg_pi,pch=20,cex=0.9,col=alpha(color1,0.05))
lines(as.integer(rownames(pi.1mb.cv1))*9.959,pi.1mb.cv1$avg_pi,lwd=0.85)
plot(rownames(pi.100kb.cv2),pi.100kb.cv2$avg_pi,pch=20,cex=0.9,col=alpha(color1,0.05))
lines(as.integer(rownames(pi.1mb.cv2))*9.959,pi.1mb.cv2$avg_pi,lwd=0.85)
plot(rownames(pi.100kb.co1),pi.100kb.co1$avg_pi,pch=20,cex=0.9,col=alpha(color2,0.05))
lines(as.integer(rownames(pi.1mb.co1))*9.959,pi.1mb.co1$avg_pi,lwd=0.85)
plot(rownames(pi.100kb.co2),pi.100kb.co2$avg_pi,pch=20,cex=0.9,col=alpha(color2,0.05))
lines(as.integer(rownames(pi.1mb.co2))*9.959,pi.1mb.co2$avg_pi,lwd=0.85)

plot(rownames(dxy.100kb.cv1co1),dxy.100kb.cv1co1$avg_dxy,pch=20,cex=0.9,col=alpha(color3,0.05)) # dummy scan for size consistency
plot(rownames(dxy.100kb.cv1co1),dxy.100kb.cv1co1$avg_dxy,pch=20,cex=0.9,col=alpha(color3,0.05))
lines(as.integer(rownames(dxy.1mb.cv1co1))*9.959,dxy.1mb.cv1co1$avg_dxy,lwd=0.85)
plot(rownames(dxy.100kb.cv1cv2),dxy.100kb.cv1cv2$avg_dxy,pch=20,cex=0.9,col=alpha(color1,0.05))
lines(as.integer(rownames(dxy.1mb.cv1cv2))*9.959,dxy.1mb.cv1cv2$avg_dxy,lwd=0.85)
plot(rownames(dxy.100kb.co1co2),dxy.100kb.co1co2$avg_dxy,pch=20,cex=0.9,col=alpha(color2,0.05))
lines(as.integer(rownames(dxy.1mb.co1co2))*9.959,dxy.1mb.co1co2$avg_dxy,lwd=0.85)

plot(rownames(fst.100kb.cv1co1),fst.100kb.cv1co1$avg_wc_fst,pch=20,cex=0.9,col=alpha(color3,0.05)) # dummy scan for size consistency
plot(rownames(fst.100kb.cv1co1),fst.100kb.cv1co1$avg_wc_fst,pch=20,cex=0.9,col=alpha(color3,0.05))
lines(as.integer(rownames(fst.1mb.cv1co1))*9.959,fst.1mb.cv1co1$avg_wc_fst,lwd=0.85)
plot(rownames(fst.100kb.cv1cv2),fst.100kb.cv1cv2$avg_wc_fst,pch=20,cex=0.9,col=alpha(color1,0.05))
lines(as.integer(rownames(fst.1mb.cv1cv2))*9.959,fst.1mb.cv1cv2$avg_wc_fst,lwd=0.85)
plot(rownames(fst.100kb.co1co2),fst.100kb.co1co2$avg_wc_fst,pch=20,cex=0.9,col=alpha(color2,0.05))
lines(as.integer(rownames(fst.1mb.co1co2))*9.959,fst.1mb.co1co2$avg_wc_fst,lwd=0.85)

### Summary statistics-------------------------------------------------------

pi.100kb.cv1.ma <- pi.100kb.cv1 %>% filter(str_detect(chromosome, "scaffold-ma"))
pi.100kb.cv1.mi <- pi.100kb.cv1 %>% filter(str_detect(chromosome, "scaffold-mi"))
pi.100kb.cv1.Z <- pi.100kb.cv1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
pi.100kb.cv1.par <- pi.100kb.cv1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

pi.100kb.cv2.ma <- pi.100kb.cv2 %>% filter(str_detect(chromosome, "scaffold-ma"))
pi.100kb.cv2.mi <- pi.100kb.cv2 %>% filter(str_detect(chromosome, "scaffold-mi"))
pi.100kb.cv2.Z <- pi.100kb.cv2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
pi.100kb.cv2.par <- pi.100kb.cv2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

pi.100kb.co1.ma <- pi.100kb.co1 %>% filter(str_detect(chromosome, "scaffold-ma"))
pi.100kb.co1.mi <- pi.100kb.co1 %>% filter(str_detect(chromosome, "scaffold-mi"))
pi.100kb.co1.Z <- pi.100kb.co1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
pi.100kb.co1.par <- pi.100kb.co1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

pi.100kb.co2.ma <- pi.100kb.co2 %>% filter(str_detect(chromosome, "scaffold-ma"))
pi.100kb.co2.mi <- pi.100kb.co2 %>% filter(str_detect(chromosome, "scaffold-mi"))
pi.100kb.co2.Z <- pi.100kb.co2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
pi.100kb.co2.par <- pi.100kb.co2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

dxy.100kb.cv1co1.ma <- dxy.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-ma"))
dxy.100kb.cv1co1.mi <- dxy.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-mi"))
dxy.100kb.cv1co1.Z <- dxy.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
dxy.100kb.cv1co1.par <- dxy.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

dxy.100kb.cv1cv2.ma <- dxy.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-ma"))
dxy.100kb.cv1cv2.mi <- dxy.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-mi"))
dxy.100kb.cv1cv2.Z <- dxy.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
dxy.100kb.cv1cv2.par <- dxy.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

dxy.100kb.co1co2.ma <- dxy.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-ma"))
dxy.100kb.co1co2.mi <- dxy.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-mi"))
dxy.100kb.co1co2.Z <- dxy.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
dxy.100kb.co1co2.par <- dxy.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

fst.100kb.cv1co1.ma <- fst.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-ma"))
fst.100kb.cv1co1.mi <- fst.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-mi"))
fst.100kb.cv1co1.Z <- fst.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
fst.100kb.cv1co1.par <- fst.100kb.cv1co1 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

fst.100kb.cv1cv2.ma <- fst.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-ma"))
fst.100kb.cv1cv2.mi <- fst.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-mi"))
fst.100kb.cv1cv2.Z <- fst.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
fst.100kb.cv1cv2.par <- fst.100kb.cv1cv2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

fst.100kb.co1co2.ma <- fst.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-ma"))
fst.100kb.co1co2.mi <- fst.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-mi"))
fst.100kb.co1co2.Z <- fst.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 < 106800000)
fst.100kb.co1co2.par <- fst.100kb.co1co2 %>% filter(str_detect(chromosome, "scaffold-Z") & window_pos_1 >= 106800000)

## Genome-wide
# pi
mean(pi.100kb.cv1$avg_pi)
sd(pi.100kb.cv1$avg_pi)
mean(pi.100kb.cv2$avg_pi)
sd(pi.100kb.cv2$avg_pi)
mean(pi.100kb.co1$avg_pi)
sd(pi.100kb.co1$avg_pi)
mean(pi.100kb.co2$avg_pi)
sd(pi.100kb.co2$avg_pi)

# dxy
mean(dxy.100kb.cv1co1$avg_dxy)
sd(dxy.100kb.cv1co1$avg_dxy)
mean(dxy.100kb.cv1cv2$avg_dxy)
sd(dxy.100kb.cv1cv2$avg_dxy)
mean(dxy.100kb.co1co2$avg_dxy)
sd(dxy.100kb.co1co2$avg_dxy)

# fst
mean(fst.100kb.cv1co1$avg_wc_fst)
sd(fst.100kb.cv1co1$avg_wc_fst)
mean(fst.100kb.cv1cv2$avg_wc_fst)
sd(fst.100kb.cv1cv2$avg_wc_fst)
mean(fst.100kb.co1co2$avg_wc_fst)
sd(fst.100kb.co1co2$avg_wc_fst)

## Macrochromosome
# pi
mean(pi.100kb.cv1.ma$avg_pi)
sd(pi.100kb.cv1.ma$avg_pi)
mean(pi.100kb.cv2.ma$avg_pi)
sd(pi.100kb.cv2.ma$avg_pi)
mean(pi.100kb.co1.ma$avg_pi)
sd(pi.100kb.co1.ma$avg_pi)
mean(pi.100kb.co2.ma$avg_pi)
sd(pi.100kb.co2.ma$avg_pi)

# dxy
mean(dxy.100kb.cv1co1.ma$avg_dxy)
sd(dxy.100kb.cv1co1.ma$avg_dxy)
mean(dxy.100kb.cv1cv2.ma$avg_dxy)
sd(dxy.100kb.cv1cv2.ma$avg_dxy)
mean(dxy.100kb.co1co2.ma$avg_dxy)
sd(dxy.100kb.co1co2.ma$avg_dxy)

# fst
mean(fst.100kb.cv1co1.ma$avg_wc_fst)
sd(fst.100kb.cv1co1.ma$avg_wc_fst)
mean(fst.100kb.cv1cv2.ma$avg_wc_fst)
sd(fst.100kb.cv1cv2.ma$avg_wc_fst)
mean(fst.100kb.co1co2.ma$avg_wc_fst)
sd(fst.100kb.co1co2.ma$avg_wc_fst)

## Microchromosome
# pi
mean(pi.100kb.cv1.mi$avg_pi)
sd(pi.100kb.cv1.mi$avg_pi)
mean(pi.100kb.cv2.mi$avg_pi)
sd(pi.100kb.cv2.mi$avg_pi)
mean(pi.100kb.co1.mi$avg_pi)
sd(pi.100kb.co1.mi$avg_pi)
mean(pi.100kb.co2.mi$avg_pi)
sd(pi.100kb.co2.mi$avg_pi)

# dxy
mean(dxy.100kb.cv1co1.mi$avg_dxy)
sd(dxy.100kb.cv1co1.mi$avg_dxy)
mean(dxy.100kb.cv1cv2.mi$avg_dxy)
sd(dxy.100kb.cv1cv2.mi$avg_dxy)
mean(dxy.100kb.co1co2.mi$avg_dxy)
sd(dxy.100kb.co1co2.mi$avg_dxy)

# fst
mean(fst.100kb.cv1co1.mi$avg_wc_fst)
sd(fst.100kb.cv1co1.mi$avg_wc_fst)
mean(fst.100kb.cv1cv2.mi$avg_wc_fst)
sd(fst.100kb.cv1cv2.mi$avg_wc_fst)
mean(fst.100kb.co1co2.mi$avg_wc_fst)
sd(fst.100kb.co1co2.mi$avg_wc_fst)

## Z chromosome
# pi
mean(pi.100kb.cv1.Z$avg_pi)
sd(pi.100kb.cv1.Z$avg_pi)
mean(pi.100kb.cv2.Z$avg_pi)
sd(pi.100kb.cv2.Z$avg_pi)
mean(pi.100kb.co1.Z$avg_pi)
sd(pi.100kb.co1.Z$avg_pi)
mean(pi.100kb.co2.Z$avg_pi)
sd(pi.100kb.co2.Z$avg_pi)

# dxy
mean(dxy.100kb.cv1co1.Z$avg_dxy)
sd(dxy.100kb.cv1co1.Z$avg_dxy)
mean(dxy.100kb.cv1cv2.Z$avg_dxy)
sd(dxy.100kb.cv1cv2.Z$avg_dxy)
mean(dxy.100kb.co1co2.Z$avg_dxy)
sd(dxy.100kb.co1co2.Z$avg_dxy)

# fst
mean(fst.100kb.cv1co1.Z$avg_wc_fst)
sd(fst.100kb.cv1co1.Z$avg_wc_fst)
mean(fst.100kb.cv1cv2.Z$avg_wc_fst)
sd(fst.100kb.cv1cv2.Z$avg_wc_fst)
mean(fst.100kb.co1co2.Z$avg_wc_fst)
sd(fst.100kb.co1co2.Z$avg_wc_fst)

## PAR
# pi
mean(pi.100kb.cv1.par$avg_pi)
sd(pi.100kb.cv1.par$avg_pi)
mean(pi.100kb.cv2.par$avg_pi)
sd(pi.100kb.cv2.par$avg_pi)
mean(pi.100kb.co1.par$avg_pi)
sd(pi.100kb.co1.par$avg_pi)
mean(pi.100kb.co2.par$avg_pi)
sd(pi.100kb.co2.par$avg_pi)

# dxy
mean(dxy.100kb.cv1co1.par$avg_dxy)
sd(dxy.100kb.cv1co1.par$avg_dxy)
mean(dxy.100kb.cv1cv2.par$avg_dxy)
sd(dxy.100kb.cv1cv2.par$avg_dxy)
mean(dxy.100kb.co1co2.par$avg_dxy)
sd(dxy.100kb.co1co2.par$avg_dxy)

# fst
mean(fst.100kb.cv1co1.par$avg_wc_fst)
sd(fst.100kb.cv1co1.par$avg_wc_fst)
mean(fst.100kb.cv1cv2.par$avg_wc_fst)
sd(fst.100kb.cv1cv2.par$avg_wc_fst)
mean(fst.100kb.co1co2.par$avg_wc_fst)
sd(fst.100kb.co1co2.par$avg_wc_fst)

### t-tests------------------------------------------------------------------

t.test(pi.100kb.cv1.ma$avg_pi,pi.100kb.cv1.mi$avg_pi)
t.test(pi.100kb.cv2.ma$avg_pi,pi.100kb.cv2.mi$avg_pi)
t.test(pi.100kb.co1.ma$avg_pi,pi.100kb.co1.mi$avg_pi)
t.test(pi.100kb.co2.ma$avg_pi,pi.100kb.co2.mi$avg_pi)

### Boxplots of chromosome distributions-------------------------------------

par(mfrow=c(1,4))
boxplot(pi.100kb.cv1.ma$avg_pi,pi.100kb.cv1.mi$avg_pi,pi.100kb.cv1.Z$avg_pi,pi.100kb.cv1.par$avg_pi,outline=F,ylim=c(0,0.0115),col=c('grey','lightblue','seagreen','aquamarine4'),ylab='pi CV1')
boxplot(pi.100kb.cv2.ma$avg_pi,pi.100kb.cv2.mi$avg_pi,pi.100kb.cv2.Z$avg_pi,pi.100kb.cv2.par$avg_pi,outline=F,ylim=c(0,0.0115),col=c('grey','lightblue','seagreen','aquamarine4'))
boxplot(pi.100kb.co1.ma$avg_pi,pi.100kb.co1.mi$avg_pi,pi.100kb.co1.Z$avg_pi,pi.100kb.co1.par$avg_pi,outline=F,ylim=c(0,0.0115),col=c('grey','lightblue','seagreen','aquamarine4'),ylab='pi CO1')
boxplot(pi.100kb.co2.ma$avg_pi,pi.100kb.co2.mi$avg_pi,pi.100kb.co2.Z$avg_pi,pi.100kb.co2.par$avg_pi,outline=F,ylim=c(0,0.0115),col=c('grey','lightblue','seagreen','aquamarine4'))

par(mfrow=c(1,3))
boxplot(dxy.100kb.cv1co1.ma$avg_dxy,dxy.100kb.cv1co1.mi$avg_dxy,dxy.100kb.cv1co1.Z$avg_dxy,dxy.100kb.cv1co1.par$avg_dxy,outline=F,ylim=c(0,0.014),col=c('grey','lightblue','seagreen','aquamarine4'),ylab='dxy CV1-CO1')
boxplot(dxy.100kb.cv1cv2.ma$avg_dxy,dxy.100kb.cv1cv2.mi$avg_dxy,dxy.100kb.cv1cv2.Z$avg_dxy,dxy.100kb.cv1cv2.par$avg_dxy,outline=F,ylim=c(0,0.014),col=c('grey','lightblue','seagreen','aquamarine4'))
boxplot(dxy.100kb.co1co2.ma$avg_dxy,dxy.100kb.co1co2.mi$avg_dxy,dxy.100kb.co1co2.Z$avg_dxy,dxy.100kb.co1co2.par$avg_dxy,outline=F,ylim=c(0,0.014),col=c('grey','lightblue','seagreen','aquamarine4'))

par(mfrow=c(1,3))
boxplot(fst.100kb.cv1co1.ma$avg_wc_fst,fst.100kb.cv1co1.mi$avg_wc_fst,fst.100kb.cv1co1.Z$avg_wc_fst,fst.100kb.cv1co1.par$avg_wc_fst,outline=F,col=c('grey','lightblue','seagreen','aquamarine4'),ylab='fst CV1-CO1')
boxplot(fst.100kb.cv1cv2.ma$avg_wc_fst,fst.100kb.cv1cv2.mi$avg_wc_fst,fst.100kb.cv1cv2.Z$avg_wc_fst,fst.100kb.cv1cv2.par$avg_wc_fst,outline=F,col=c('grey','lightblue','seagreen','aquamarine4'))
boxplot(fst.100kb.co1co2.ma$avg_wc_fst,fst.100kb.co1co2.mi$avg_wc_fst,fst.100kb.co1co2.Z$avg_wc_fst,fst.100kb.co1co2.par$avg_wc_fst,outline=F,col=c('grey','lightblue','seagreen','aquamarine4'))

### Correlation coefficients-------------------------------------------------

# Compare within statistics
cor.test(pi.100kb.cv1$avg_pi,pi.100kb.co1$avg_pi,method='spearman')
cor.test(pi.100kb.cv1$avg_pi,pi.100kb.cv2$avg_pi,method='spearman')
cor.test(pi.100kb.cv1$avg_pi,pi.100kb.co2$avg_pi,method='spearman')
cor.test(pi.100kb.co1$avg_pi,pi.100kb.co2$avg_pi,method='spearman')
cor.test(pi.100kb.co1$avg_pi,pi.100kb.cv2$avg_pi,method='spearman')
cor.test(pi.100kb.cv2$avg_pi,pi.100kb.co2$avg_pi,method='spearman')

cor.test(dxy.100kb.cv1co1$avg_dxy,dxy.100kb.cv1cv2$avg_dxy,method='spearman')
cor.test(dxy.100kb.cv1co1$avg_dxy,dxy.100kb.co1co2$avg_dxy,method='spearman')
cor.test(dxy.100kb.cv1cv2$avg_dxy,dxy.100kb.co1co2$avg_dxy,method='spearman')

cor.test(fst.100kb.cv1co1$avg_wc_fst,fst.100kb.cv1cv2$avg_wc_fst,method='spearman')
cor.test(fst.100kb.cv1co1$avg_wc_fst,fst.100kb.co1co2$avg_wc_fst,method='spearman')
cor.test(fst.100kb.cv1cv2$avg_wc_fst,fst.100kb.co1co2$avg_wc_fst,method='spearman')

cor.test(rho.100kb.cv1$mean,rho.100kb.co1$mean,method='spearman')

# Compare between statistics
cor.test(pi.100kb.cv1$avg_pi,dxy.100kb.cv1cv2$avg_dxy)
cor.test(pi.100kb.cv1$avg_pi,dxy.100kb.cv1co1$avg_dxy)
cor.test(pi.100kb.cv2$avg_pi,dxy.100kb.cv1cv2$avg_dxy)
cor.test(pi.100kb.co1$avg_pi,dxy.100kb.co1co2$avg_dxy)
cor.test(pi.100kb.co1$avg_pi,dxy.100kb.cv1co1$avg_dxy)
cor.test(pi.100kb.co2$avg_pi,dxy.100kb.co1co2$avg_dxy)

cor.test(pi.100kb.cv1$avg_pi,fst.100kb.cv1cv2$avg_wc_fst)
cor.test(pi.100kb.cv1$avg_pi,fst.100kb.cv1co1$avg_wc_fst)
cor.test(pi.100kb.cv2$avg_pi,fst.100kb.cv1cv2$avg_wc_fst)
cor.test(pi.100kb.co1$avg_pi,fst.100kb.co1co2$avg_wc_fst)
cor.test(pi.100kb.co1$avg_pi,fst.100kb.cv1co1$avg_wc_fst)
cor.test(pi.100kb.co2$avg_pi,fst.100kb.co1co2$avg_wc_fst)

cor.test(dxy.100kb.cv1co1$avg_dxy,fst.100kb.cv1co1$avg_wc_fst)
cor.test(dxy.100kb.cv1cv2$avg_dxy,fst.100kb.cv1cv2$avg_wc_fst)
cor.test(dxy.100kb.co1co2$avg_dxy,fst.100kb.co1co2$avg_wc_fst)

par(mfrow=c(2,3))
# Export as 5.25 x 7
plot(pi.100kb.cv1$avg_pi,dxy.100kb.cv1co1$avg_dxy,pch=20,cex=0.8,col=alpha('black',0.1))
plot(pi.100kb.cv1$avg_pi,fst.100kb.cv1co1$avg_wc_fst,pch=20,cex=0.8,col=alpha('black',0.1))
plot(dxy.100kb.cv1co1$avg_dxy,fst.100kb.cv1co1$avg_wc_fst,pch=20,cex=0.8,col=alpha('black',0.1))
plot(pi.100kb.co1$avg_pi,dxy.100kb.cv1co1$avg_dxy,pch=20,cex=0.8,col=alpha('black',0.1))
plot(pi.100kb.co1$avg_pi,fst.100kb.cv1co1$avg_wc_fst,pch=20,cex=0.8,col=alpha('black',0.1))
plot(pi.100kb.cv1$avg_pi,pi.100kb.co1$avg_pi,pch=20,cex=0.8,col=alpha('black',0.1))

# Compare with recombination rates
cor.test(pi.100kb.cv1$avg_pi,rho.100kb.cv1$mean,method='spearman')
cor.test(pi.100kb.co1$avg_pi,rho.100kb.co1$mean,method='spearman')
cor.test(pi.100kb.cv2$avg_pi,rho.100kb.cv1$mean,method='spearman')
cor.test(pi.100kb.co2$avg_pi,rho.100kb.co1$mean,method='spearman')

cor.test(dxy.100kb.cv1co1$avg_dxy,rho.100kb.cv1$mean,method='spearman')
cor.test(dxy.100kb.cv1co1$avg_dxy,rho.100kb.co1$mean,method='spearman')
cor.test(dxy.100kb.cv1cv2$avg_dxy,rho.100kb.cv1$mean,method='spearman')
cor.test(dxy.100kb.co1co2$avg_dxy,rho.100kb.co1$mean,method='spearman')

cor.test(fst.100kb.cv1co1$avg_wc_fst,rho.100kb.cv1$mean,method='spearman')
cor.test(fst.100kb.cv1co1$avg_wc_fst,rho.100kb.co1$mean,method='spearman')
cor.test(fst.100kb.cv1cv2$avg_wc_fst,rho.100kb.cv1$mean,method='spearman')
cor.test(fst.100kb.co1co2$avg_wc_fst,rho.100kb.co1$mean,method='spearman')

plot(pi.100kb.cv1$avg_pi,rho.100kb.cv1$mean,pch=20,cex=0.8,col=alpha('black',0.1))
plot(pi.100kb.co1$avg_pi,rho.100kb.co1$mean,pch=20,cex=0.8,col=alpha('black',0.1),ylim=c(0,0.2))
plot(dxy.100kb.cv1co1$avg_dxy,rho.100kb.cv1$mean,pch=20,cex=0.8,col=alpha('black',0.1))
plot(fst.100kb.cv1co1$avg_wc_fst,rho.100kb.cv1$mean,pch=20,cex=0.8,col=alpha('black',0.1))
