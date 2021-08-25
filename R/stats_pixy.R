############################################################################
# Summary statistics and tests on genetic diversity and differnetiation
# in venom gene regions
############################################################################

### Goal: compare estimates of pi, dxy, and fst in venom gene regions to the
### genome/chromosomal background, and calculate summary statistics per gene

### i. Clear environment----------------------------------------------------

rm(list = ls())

### ii. Set working directory & dependencies--------------------------------

setwd('./pixy/')

### 1. Read in general data locations---------------------------------------

# Gene locations
svmp <- read.table('./resources/gene_SVMP.bed',header=F)
svsp <- read.table('./resources/gene_SVSP.bed',header=F)
pla2 <- read.table('./resources/gene_PLA2.bed',header=F)
# Region locations
svmp.reg <- read.table('./resources/region_SVMP_scaffold-mi1.bed',header=F)
svsp.reg <- read.table('./resources/region_SVSP_scaffold-mi2.bed',header=F)
pla2.reg <- read.table('./resources/region_PLA2_scaffold-mi7.bed',header=F)

### 2. Calculate genome-wide means------------------------------------------

pi.10kb <- read.table('pixy.all.10kb_pi.txt',header=T)
dxy.10kb <- read.table('pixy.all.10kb_dxy.txt',header=T)
fst.10kb <- read.table('pixy.all.10kb_fst.txt',header=T)

pi.10kb.cv1 <- pi.10kb %>% filter(str_detect(pop,'CV1'))
pi.10kb.cv2 <- pi.10kb %>% filter(str_detect(pop,'CV2'))
pi.10kb.co1 <- pi.10kb %>% filter(str_detect(pop,'CO1'))
pi.10kb.co2 <- pi.10kb %>% filter(str_detect(pop,'CO2'))

dxy.10kb.cv1co1 <- dxy.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CO1'))
dxy.10kb.cv1cv2 <- dxy.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CV2'))
dxy.10kb.co1co2 <- dxy.10kb %>% filter(str_detect(pop1, 'CO1') & str_detect(pop2, 'CO2'))

fst.10kb.cv1co1 <- fst.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CO1'))
fst.10kb.cv1cv2 <- fst.10kb %>% filter(str_detect(pop1, 'CV1') & str_detect(pop2, 'CV2'))
fst.10kb.co1co2 <- fst.10kb %>% filter(str_detect(pop1, 'CO1') & str_detect(pop2, 'CO2'))

cv1pi.avg <- mean(pi.10kb.cv1$avg_pi)
cv2pi.avg <- mean(pi.10kb.cv2$avg_pi)
co1pi.avg <- mean(pi.10kb.co1$avg_pi)
co2pi.avg <- mean(pi.10kb.co2$avg_pi)

cv1cv2.dxy.avg <- mean(dxy.10kb.cv1cv2$avg_dxy,na.rm=T)
cv1co2.dxy.avg <- mean(dxy.10kb.co1co2$avg_dxy,na.rm=T)
cv1co1.dxy.avg <- mean(dxy.10kb.cv1co1$avg_dxy,na.rm=T)

cv1cv2.fst.avg <- mean(fst.10kb.cv1cv2$avg_wc_fst,na.rm=T)
cv1co2.fst.avg <- mean(fst.10kb.co1co2$avg_wc_fst,na.rm=T)
cv1co1.fst.avg <- mean(fst.10kb.cv1co1$avg_wc_fst,na.rm=T)

### 3. Read in data (1 kb windows)------------------------------------------

pi.mi1.1kb <- read.table('pixy.scaffold-mi1.1kb_pi.txt',header=T)
dxy.mi1.1kb <- read.table('pixy.scaffold-mi1.1kb_dxy.txt',header=T)
fst.mi1.1kb <- read.table('pixy.scaffold-mi1.1kb_fst.txt',header=T)

pi.mi2.1kb <- read.table('pixy.scaffold-mi2.1kb_pi.txt',header=T)
dxy.mi2.1kb <- read.table('pixy.scaffold-mi2.1kb_dxy.txt',header=T)
fst.mi2.1kb <- read.table('pixy.scaffold-mi2.1kb_fst.txt',header=T)

pi.mi7.1kb <- read.table('pixy.scaffold-mi7.1kb_pi.txt',header=T)
dxy.mi7.1kb <- read.table('pixy.scaffold-mi7.1kb_dxy.txt',header=T)
fst.mi7.1kb <- read.table('pixy.scaffold-mi7.1kb_fst.txt',header=T)

### 3.1 Read in non-venom paralog values (10kb and 1kb windows)-------------

cv1.pi.nv.svmp <- read.table('./non-venom/CV1.non-venom_SVMP.10kb_pi.txt',header=T)
cv2.pi.nv.svmp <- read.table('./non-venom/CV2.non-venom_SVMP.10kb_pi.txt',header=T)
co1.pi.nv.svmp <- read.table('./non-venom/CO1.non-venom_SVMP.10kb_pi.txt',header=T)
co2.pi.nv.svmp <- read.table('./non-venom/CO2.non-venom_SVMP.10kb_pi.txt',header=T)

cv1.pi.nv.svsp <- read.table('./non-venom/CV1.non-venom_SVSP.10kb_pi.txt',header=T)
cv2.pi.nv.svsp <- read.table('./non-venom/CV2.non-venom_SVSP.10kb_pi.txt',header=T)
co1.pi.nv.svsp <- read.table('./non-venom/CO1.non-venom_SVSP.10kb_pi.txt',header=T)
co2.pi.nv.svsp <- read.table('./non-venom/CO2.non-venom_SVSP.10kb_pi.txt',header=T)

cv1.pi.nv.pla2 <- read.table('./non-venom/CV1.non-venom_PLA2.1kb_pi.txt',header=T)
cv2.pi.nv.pla2 <- read.table('./non-venom/CV2.non-venom_PLA2.1kb_pi.txt',header=T)
co1.pi.nv.pla2 <- read.table('./non-venom/CO1.non-venom_PLA2.1kb_pi.txt',header=T)
co2.pi.nv.pla2 <- read.table('./non-venom/CO2.non-venom_PLA2.1kb_pi.txt',header=T)

cv1co1.dxy.nv.svmp <- read.table('./non-venom/cv1co1.non-venom_SVMP.10kb_dxy.txt',header=T)
cv1cv2.dxy.nv.svmp <- read.table('./non-venom/cv1cv2.non-venom_SVMP.10kb_dxy.txt',header=T)
co1co2.dxy.nv.svmp <- read.table('./non-venom/co1co2.non-venom_SVMP.10kb_dxy.txt',header=T)

cv1co1.dxy.nv.svsp <- read.table('./non-venom/cv1co1.non-venom_SVSP.10kb_dxy.txt',header=T)
cv1cv2.dxy.nv.svsp <- read.table('./non-venom/cv1cv2.non-venom_SVSP.10kb_dxy.txt',header=T)
co1co2.dxy.nv.svsp <- read.table('./non-venom/co1co2.non-venom_SVSP.10kb_dxy.txt',header=T)

cv1co1.dxy.nv.pla2 <- read.table('./non-venom/cv1co1.non-venom_PLA2.1kb_dxy.txt',header=T)
cv1cv2.dxy.nv.pla2 <- read.table('./non-venom/cv1cv2.non-venom_PLA2.1kb_dxy.txt',header=T)
co1co2.dxy.nv.pla2 <- read.table('./non-venom/co1co2.non-venom_PLA2.1kb_dxy.txt',header=T)

cv1co1.fst.nv.svmp <- read.table('./non-venom/cv1co1.non-venom_SVMP.10kb_fst.txt',header=T)
cv1cv2.fst.nv.svmp <- read.table('./non-venom/cv1cv2.non-venom_SVMP.10kb_fst.txt',header=T)
co1co2.fst.nv.svmp <- read.table('./non-venom/co1co2.non-venom_SVMP.10kb_fst.txt',header=T)

cv1co1.fst.nv.svsp <- read.table('./non-venom/cv1co1.non-venom_SVSP.10kb_fst.txt',header=T)
cv1cv2.fst.nv.svsp <- read.table('./non-venom/cv1cv2.non-venom_SVSP.10kb_fst.txt',header=T)
co1co2.fst.nv.svsp <- read.table('./non-venom/co1co2.non-venom_SVSP.10kb_fst.txt',header=T)

cv1co1.fst.nv.pla2 <- read.table('./non-venom/cv1co1.non-venom_PLA2.1kb_fst.txt',header=T)
cv1cv2.fst.nv.pla2 <- read.table('./non-venom/cv1cv2.non-venom_PLA2.1kb_fst.txt',header=T)
co1co2.fst.nv.pla2 <- read.table('./non-venom/co1co2.non-venom_PLA2.1kb_fst.txt',header=T)

### 4. Parse data and read in cnv-masked data for C. oreganus---------------

# chromosome 9
cv1.pi.mi1 <- pi.mi1.1kb[which(pi.mi1.1kb$pop=='CV1'),]
cv2.pi.mi1 <- pi.mi1.1kb[which(pi.mi1.1kb$pop=='CV2'),]
co1.pi.mi1 <- read.table('pi.scaffold-mi1.1kb.co1.cnvMask.txt',header=T)
co2.pi.mi1 <- read.table('pi.scaffold-mi1.1kb.co2.cnvMask.txt',header=T)

cv1cv2.dxy.mi1 <- dxy.mi1.1kb[which(dxy.mi1.1kb$pop1=='CV1' & dxy.mi1.1kb$pop2=='CV2'),]
cv1co1.dxy.mi1 <- read.table('dxy.scaffold-mi1.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi1 <- read.table('dxy.scaffold-mi1.1kb.co1co2.cnvMask.txt',header=T)

cv1cv2.fst.mi1 <- fst.mi1.1kb[which(fst.mi1.1kb$pop1=='CV1' & fst.mi1.1kb$pop2=='CV2'),]
cv1co1.fst.mi1 <- read.table('fst.scaffold-mi1.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi1 <- read.table('fst.scaffold-mi1.1kb.co1co2.cnvMask.txt',header=T)

# chromosome 10
cv1.pi.mi2 <- pi.mi2.1kb[which(pi.mi2.1kb$pop=='CV1'),]
cv2.pi.mi2 <- pi.mi2.1kb[which(pi.mi2.1kb$pop=='CV2'),]
co1.pi.mi2 <- read.table('pi.scaffold-mi2.1kb.co1.cnvMask.txt',header=T)
co2.pi.mi2 <- read.table('pi.scaffold-mi2.1kb.co2.cnvMask.txt',header=T)
cv1cv2.dxy.mi2 <- dxy.mi2.1kb[which(dxy.mi2.1kb$pop1=='CV1' & dxy.mi2.1kb$pop2=='CV2'),]
cv1co1.dxy.mi2 <- read.table('dxy.scaffold-mi2.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi2 <- read.table('dxy.scaffold-mi2.1kb.co1co2.cnvMask.txt',header=T)
cv1cv2.fst.mi2 <- fst.mi2.1kb[which(fst.mi2.1kb$pop1=='CV1' & fst.mi2.1kb$pop2=='CV2'),]
cv1co1.fst.mi2 <- read.table('fst.scaffold-mi2.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi2 <- read.table('fst.scaffold-mi2.1kb.co1co2.cnvMask.txt',header=T)

# chromosome 15
cv1.pi.mi7 <- pi.mi7.1kb[which(pi.mi7.1kb$pop=='CV1'),]
cv2.pi.mi7 <- pi.mi7.1kb[which(pi.mi7.1kb$pop=='CV2'),]
co1.pi.mi7 <- read.table('pi.scaffold-mi7.1kb.co1.cnvMask.txt',header=T)
co2.pi.mi7 <- read.table('pi.scaffold-mi7.1kb.co2.cnvMask.txt',header=T)
cv1cv2.dxy.mi7 <- dxy.mi7.1kb[which(dxy.mi7.1kb$pop1=='CV1' & dxy.mi7.1kb$pop2=='CV2'),]
cv1co1.dxy.mi7 <- read.table('dxy.scaffold-mi7.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.dxy.mi7 <- read.table('dxy.scaffold-mi7.1kb.co1co2.cnvMask.txt',header=T)
cv1cv2.fst.mi7 <- fst.mi7.1kb[which(fst.mi7.1kb$pop1=='CV1' & fst.mi7.1kb$pop2=='CV2'),]
cv1co1.fst.mi7 <- read.table('fst.scaffold-mi7.1kb.cv1co1.cnvMask.txt',header=T)
co1co2.fst.mi7 <- read.table('fst.scaffold-mi7.1kb.co1co2.cnvMask.txt',header=T)

### 5. Extract estimates in/out of venom regions----------------------------

# SVMP and chromosome 9 background
cv1.pi.svmp <- cv1.pi.mi1[which(cv1.pi.mi1$window_pos_1>=svmp.reg$V2 & cv1.pi.mi1$window_pos_2 <= svmp.reg$V3),]
cv1.pi.mi1.bk <- cv1.pi.mi1[which(cv1.pi.mi1$window_pos_1<svmp.reg$V2 & cv1.pi.mi1$window_pos_2 <= svmp.reg$V2 | cv1.pi.mi1$window_pos_1>svmp.reg$V3 & cv1.pi.mi1$window_pos_2 >= svmp.reg$V3),]
cv2.pi.svmp <- cv2.pi.mi1[which(cv2.pi.mi1$window_pos_1>=svmp.reg$V2 & cv2.pi.mi1$window_pos_2 <= svmp.reg$V3),]
cv2.pi.mi1.bk <- cv2.pi.mi1[which(cv2.pi.mi1$window_pos_1<svmp.reg$V2 & cv2.pi.mi1$window_pos_2 <= svmp.reg$V2 | cv2.pi.mi1$window_pos_1>svmp.reg$V3 & cv2.pi.mi1$window_pos_2 >= svmp.reg$V3),]
co1.pi.svmp <- co1.pi.mi1[which(co1.pi.mi1$window_pos_1>=svmp.reg$V2 & co1.pi.mi1$window_pos_2 <= svmp.reg$V3),]
co1.pi.mi1.bk <- co1.pi.mi1[which(co1.pi.mi1$window_pos_1<svmp.reg$V2 & co1.pi.mi1$window_pos_2 <= svmp.reg$V2 | co1.pi.mi1$window_pos_1>svmp.reg$V3 & co1.pi.mi1$window_pos_2 >= svmp.reg$V3),]
co2.pi.svmp <- co2.pi.mi1[which(co2.pi.mi1$window_pos_1>=svmp.reg$V2 & co2.pi.mi1$window_pos_2 <= svmp.reg$V3),]
co2.pi.mi1.bk <- co2.pi.mi1[which(co2.pi.mi1$window_pos_1<svmp.reg$V2 & co2.pi.mi1$window_pos_2 <= svmp.reg$V2 | co2.pi.mi1$window_pos_1>svmp.reg$V3 & co2.pi.mi1$window_pos_2 >= svmp.reg$V3),]

cv1cv2.dxy.svmp <- cv1cv2.dxy.mi1[which(cv1cv2.dxy.mi1$window_pos_1>=svmp.reg$V2 & cv1cv2.dxy.mi1$window_pos_2 <= svmp.reg$V3),]
cv1cv2.dxy.mi1.bk <- cv1cv2.dxy.mi1[which(cv1cv2.dxy.mi1$window_pos_1<svmp.reg$V2 & cv1cv2.dxy.mi1$window_pos_2 <= svmp.reg$V2 | cv1cv2.dxy.mi1$window_pos_1>svmp.reg$V3 & cv1cv2.dxy.mi1$window_pos_2 >= svmp.reg$V3),]
cv1co1.dxy.svmp <- cv1co1.dxy.mi1[which(cv1co1.dxy.mi1$window_pos_1>=svmp.reg$V2 & cv1co1.dxy.mi1$window_pos_2 <= svmp.reg$V3),]
cv1co1.dxy.mi1.bk <- cv1co1.dxy.mi1[which(cv1co1.dxy.mi1$window_pos_1<svmp.reg$V2 & cv1co1.dxy.mi1$window_pos_2 <= svmp.reg$V2 | cv1co1.dxy.mi1$window_pos_1>svmp.reg$V3 & cv1co1.dxy.mi1$window_pos_2 >= svmp.reg$V3),]
co1co2.dxy.svmp <- co1co2.dxy.mi1[which(co1co2.dxy.mi1$window_pos_1>=svmp.reg$V2 & co1co2.dxy.mi1$window_pos_2 <= svmp.reg$V3),]
co1co2.dxy.mi1.bk <- co1co2.dxy.mi1[which(co1co2.dxy.mi1$window_pos_1<svmp.reg$V2 & co1co2.dxy.mi1$window_pos_2 <= svmp.reg$V2 | co1co2.dxy.mi1$window_pos_1>svmp.reg$V3 & co1co2.dxy.mi1$window_pos_2 >= svmp.reg$V3),]

cv1cv2.fst.svmp <- cv1cv2.fst.mi1[which(cv1cv2.fst.mi1$window_pos_1>=svmp.reg$V2 & cv1cv2.fst.mi1$window_pos_2 <= svmp.reg$V3),]
cv1cv2.fst.mi1.bk <- cv1cv2.fst.mi1[which(cv1cv2.fst.mi1$window_pos_1<svmp.reg$V2 & cv1cv2.fst.mi1$window_pos_2 <= svmp.reg$V2 | cv1cv2.fst.mi1$window_pos_1>svmp.reg$V3 & cv1cv2.fst.mi1$window_pos_2 >= svmp.reg$V3),]
cv1co1.fst.svmp <- cv1co1.fst.mi1[which(cv1co1.fst.mi1$window_pos_1>=svmp.reg$V2 & cv1co1.fst.mi1$window_pos_2 <= svmp.reg$V3),]
cv1co1.fst.mi1.bk <- cv1co1.fst.mi1[which(cv1co1.fst.mi1$window_pos_1<svmp.reg$V2 & cv1co1.fst.mi1$window_pos_2 <= svmp.reg$V2 | cv1co1.fst.mi1$window_pos_1>svmp.reg$V3 & cv1co1.fst.mi1$window_pos_2 >= svmp.reg$V3),]
co1co2.fst.svmp <- co1co2.fst.mi1[which(co1co2.fst.mi1$window_pos_1>=svmp.reg$V2 & co1co2.fst.mi1$window_pos_2 <= svmp.reg$V3),]
co1co2.fst.mi1.bk <- co1co2.fst.mi1[which(co1co2.fst.mi1$window_pos_1<svmp.reg$V2 & co1co2.fst.mi1$window_pos_2 <= svmp.reg$V2 | co1co2.fst.mi1$window_pos_1>svmp.reg$V3 & co1co2.fst.mi1$window_pos_2 >= svmp.reg$V3),]

# SVSP and chromosome 10 background
cv1.pi.svsp <- cv1.pi.mi2[which(cv1.pi.mi2$window_pos_1>=svsp.reg$V2 & cv1.pi.mi2$window_pos_2 <= svsp.reg$V3),]
cv1.pi.mi2.bk <- cv1.pi.mi2[which(cv1.pi.mi2$window_pos_1<svsp.reg$V2 & cv1.pi.mi2$window_pos_2 <= svsp.reg$V2 | cv1.pi.mi2$window_pos_1>svsp.reg$V3 & cv1.pi.mi2$window_pos_2 >= svsp.reg$V3),]
cv2.pi.svsp <- cv2.pi.mi2[which(cv2.pi.mi2$window_pos_1>=svsp.reg$V2 & cv2.pi.mi2$window_pos_2 <= svsp.reg$V3),]
cv2.pi.mi2.bk <- cv2.pi.mi2[which(cv2.pi.mi2$window_pos_1<svsp.reg$V2 & cv2.pi.mi2$window_pos_2 <= svsp.reg$V2 | cv2.pi.mi2$window_pos_1>svsp.reg$V3 & cv2.pi.mi2$window_pos_2 >= svsp.reg$V3),]
co1.pi.svsp <- co1.pi.mi2[which(co1.pi.mi2$window_pos_1>=svsp.reg$V2 & co1.pi.mi2$window_pos_2 <= svsp.reg$V3),]
co1.pi.mi2.bk <- co1.pi.mi2[which(co1.pi.mi2$window_pos_1<svsp.reg$V2 & co1.pi.mi2$window_pos_2 <= svsp.reg$V2 | co1.pi.mi2$window_pos_1>svsp.reg$V3 & co1.pi.mi2$window_pos_2 >= svsp.reg$V3),]
co2.pi.svsp <- co2.pi.mi2[which(co2.pi.mi2$window_pos_1>=svsp.reg$V2 & co2.pi.mi2$window_pos_2 <= svsp.reg$V3),]
co2.pi.mi2.bk <- co2.pi.mi2[which(co2.pi.mi2$window_pos_1<svsp.reg$V2 & co2.pi.mi2$window_pos_2 <= svsp.reg$V2 | co2.pi.mi2$window_pos_1>svsp.reg$V3 & co2.pi.mi2$window_pos_2 >= svsp.reg$V3),]

cv1cv2.dxy.svsp <- cv1cv2.dxy.mi2[which(cv1cv2.dxy.mi2$window_pos_1>=svsp.reg$V2 & cv1cv2.dxy.mi2$window_pos_2 <= svsp.reg$V3),]
cv1cv2.dxy.mi2.bk <- cv1cv2.dxy.mi2[which(cv1cv2.dxy.mi2$window_pos_1<svsp.reg$V2 & cv1cv2.dxy.mi2$window_pos_2 <= svsp.reg$V2 | cv1cv2.dxy.mi2$window_pos_1>svsp.reg$V3 & cv1cv2.dxy.mi2$window_pos_2 >= svsp.reg$V3),]
cv1co1.dxy.svsp <- cv1co1.dxy.mi2[which(cv1co1.dxy.mi2$window_pos_1>=svsp.reg$V2 & cv1co1.dxy.mi2$window_pos_2 <= svsp.reg$V3),]
cv1co1.dxy.mi2.bk <- cv1co1.dxy.mi2[which(cv1co1.dxy.mi2$window_pos_1<svsp.reg$V2 & cv1co1.dxy.mi2$window_pos_2 <= svsp.reg$V2 | cv1co1.dxy.mi2$window_pos_1>svsp.reg$V3 & cv1co1.dxy.mi2$window_pos_2 >= svsp.reg$V3),]
co1co2.dxy.svsp <- co1co2.dxy.mi2[which(co1co2.dxy.mi2$window_pos_1>=svsp.reg$V2 & co1co2.dxy.mi2$window_pos_2 <= svsp.reg$V3),]
co1co2.dxy.mi2.bk <- co1co2.dxy.mi2[which(co1co2.dxy.mi2$window_pos_1<svsp.reg$V2 & co1co2.dxy.mi2$window_pos_2 <= svsp.reg$V2 | co1co2.dxy.mi2$window_pos_1>svsp.reg$V3 & co1co2.dxy.mi2$window_pos_2 >= svsp.reg$V3),]

cv1cv2.fst.svsp <- cv1cv2.fst.mi2[which(cv1cv2.fst.mi2$window_pos_1>=svsp.reg$V2 & cv1cv2.fst.mi2$window_pos_2 <= svsp.reg$V3),]
cv1cv2.fst.mi2.bk <- cv1cv2.fst.mi2[which(cv1cv2.fst.mi2$window_pos_1<svsp.reg$V2 & cv1cv2.fst.mi2$window_pos_2 <= svsp.reg$V2 | cv1cv2.fst.mi2$window_pos_1>svsp.reg$V3 & cv1cv2.fst.mi2$window_pos_2 >= svsp.reg$V3),]
cv1co1.fst.svsp <- cv1co1.fst.mi2[which(cv1co1.fst.mi2$window_pos_1>=svsp.reg$V2 & cv1co1.fst.mi2$window_pos_2 <= svsp.reg$V3),]
cv1co1.fst.mi2.bk <- cv1co1.fst.mi2[which(cv1co1.fst.mi2$window_pos_1<svsp.reg$V2 & cv1co1.fst.mi2$window_pos_2 <= svsp.reg$V2 | cv1co1.fst.mi2$window_pos_1>svsp.reg$V3 & cv1co1.fst.mi2$window_pos_2 >= svsp.reg$V3),]
co1co2.fst.svsp <- co1co2.fst.mi2[which(co1co2.fst.mi2$window_pos_1>=svsp.reg$V2 & co1co2.fst.mi2$window_pos_2 <= svsp.reg$V3),]
co1co2.fst.mi2.bk <- co1co2.fst.mi2[which(co1co2.fst.mi2$window_pos_1<svsp.reg$V2 & co1co2.fst.mi2$window_pos_2 <= svsp.reg$V2 | co1co2.fst.mi2$window_pos_1>svsp.reg$V3 & co1co2.fst.mi2$window_pos_2 >= svsp.reg$V3),]

# PLA2 and chromosome 15 background
cv1.pi.pla2 <- cv1.pi.mi7[which(cv1.pi.mi7$window_pos_1>=pla2.reg$V2 & cv1.pi.mi7$window_pos_2 <= pla2.reg$V3),]
cv1.pi.mi7.bk <- cv1.pi.mi7[which(cv1.pi.mi7$window_pos_1<pla2.reg$V2 & cv1.pi.mi7$window_pos_2 <= pla2.reg$V2 | cv1.pi.mi7$window_pos_1>pla2.reg$V3 & cv1.pi.mi7$window_pos_2 >= pla2.reg$V3),]
cv2.pi.pla2 <- cv2.pi.mi7[which(cv2.pi.mi7$window_pos_1>=pla2.reg$V2 & cv2.pi.mi7$window_pos_2 <= pla2.reg$V3),]
cv2.pi.mi7.bk <- cv2.pi.mi7[which(cv2.pi.mi7$window_pos_1<pla2.reg$V2 & cv2.pi.mi7$window_pos_2 <= pla2.reg$V2 | cv2.pi.mi7$window_pos_1>pla2.reg$V3 & cv2.pi.mi7$window_pos_2 >= pla2.reg$V3),]
co1.pi.pla2 <- co1.pi.mi7[which(co1.pi.mi7$window_pos_1>=pla2.reg$V2 & co1.pi.mi7$window_pos_2 <= pla2.reg$V3),]
co1.pi.mi7.bk <- co1.pi.mi7[which(co1.pi.mi7$window_pos_1<pla2.reg$V2 & co1.pi.mi7$window_pos_2 <= pla2.reg$V2 | co1.pi.mi7$window_pos_1>pla2.reg$V3 & co1.pi.mi7$window_pos_2 >= pla2.reg$V3),]
co2.pi.pla2 <- co2.pi.mi7[which(co2.pi.mi7$window_pos_1>=pla2.reg$V2 & co2.pi.mi7$window_pos_2 <= pla2.reg$V3),]
co2.pi.mi7.bk <- co2.pi.mi7[which(co2.pi.mi7$window_pos_1<pla2.reg$V2 & co2.pi.mi7$window_pos_2 <= pla2.reg$V2 | co2.pi.mi7$window_pos_1>pla2.reg$V3 & co2.pi.mi7$window_pos_2 >= pla2.reg$V3),]

cv1cv2.dxy.pla2 <- cv1cv2.dxy.mi7[which(cv1cv2.dxy.mi7$window_pos_1>=pla2.reg$V2 & cv1cv2.dxy.mi7$window_pos_2 <= pla2.reg$V3),]
cv1cv2.dxy.mi7.bk <- cv1cv2.dxy.mi7[which(cv1cv2.dxy.mi7$window_pos_1<pla2.reg$V2 & cv1cv2.dxy.mi7$window_pos_2 <= pla2.reg$V2 | cv1cv2.dxy.mi7$window_pos_1>pla2.reg$V3 & cv1cv2.dxy.mi7$window_pos_2 >= pla2.reg$V3),]
cv1co1.dxy.pla2 <- cv1co1.dxy.mi7[which(cv1co1.dxy.mi7$window_pos_1>=pla2.reg$V2 & cv1co1.dxy.mi7$window_pos_2 <= pla2.reg$V3),]
cv1co1.dxy.mi7.bk <- cv1co1.dxy.mi7[which(cv1co1.dxy.mi7$window_pos_1<pla2.reg$V2 & cv1co1.dxy.mi7$window_pos_2 <= pla2.reg$V2 | cv1co1.dxy.mi7$window_pos_1>pla2.reg$V3 & cv1co1.dxy.mi7$window_pos_2 >= pla2.reg$V3),]
co1co2.dxy.pla2 <- co1co2.dxy.mi7[which(co1co2.dxy.mi7$window_pos_1>=pla2.reg$V2 & co1co2.dxy.mi7$window_pos_2 <= pla2.reg$V3),]
co1co2.dxy.mi7.bk <- co1co2.dxy.mi7[which(co1co2.dxy.mi7$window_pos_1<pla2.reg$V2 & co1co2.dxy.mi7$window_pos_2 <= pla2.reg$V2 | co1co2.dxy.mi7$window_pos_1>pla2.reg$V3 & co1co2.dxy.mi7$window_pos_2 >= pla2.reg$V3),]

cv1cv2.fst.pla2 <- cv1cv2.fst.mi7[which(cv1cv2.fst.mi7$window_pos_1>=pla2.reg$V2 & cv1cv2.fst.mi7$window_pos_2 <= pla2.reg$V3),]
cv1cv2.fst.mi7.bk <- cv1cv2.fst.mi7[which(cv1cv2.fst.mi7$window_pos_1<pla2.reg$V2 & cv1cv2.fst.mi7$window_pos_2 <= pla2.reg$V2 | cv1cv2.fst.mi7$window_pos_1>pla2.reg$V3 & cv1cv2.fst.mi7$window_pos_2 >= pla2.reg$V3),]
cv1co1.fst.pla2 <- cv1co1.fst.mi7[which(cv1co1.fst.mi7$window_pos_1>=pla2.reg$V2 & cv1co1.fst.mi7$window_pos_2 <= pla2.reg$V3),]
cv1co1.fst.mi7.bk <- cv1co1.fst.mi7[which(cv1co1.fst.mi7$window_pos_1<pla2.reg$V2 & cv1co1.fst.mi7$window_pos_2 <= pla2.reg$V2 | cv1co1.fst.mi7$window_pos_1>pla2.reg$V3 & cv1co1.fst.mi7$window_pos_2 >= pla2.reg$V3),]
co1co2.fst.pla2 <- co1co2.fst.mi7[which(co1co2.fst.mi7$window_pos_1>=pla2.reg$V2 & co1co2.fst.mi7$window_pos_2 <= pla2.reg$V3),]
co1co2.fst.mi7.bk <- co1co2.fst.mi7[which(co1co2.fst.mi7$window_pos_1<pla2.reg$V2 & co1co2.fst.mi7$window_pos_2 <= pla2.reg$V2 | co1co2.fst.mi7$window_pos_1>pla2.reg$V3 & co1co2.fst.mi7$window_pos_2 >= pla2.reg$V3),]

### 6. Summary statistics for venom gene regions----------------------------

# SVMP
mean(cv1.pi.svmp$avg_pi,na.rm=T)
mean(cv2.pi.svmp$avg_pi,na.rm=T)
mean(co1.pi.svmp$avg_pi,na.rm=T)
mean(co2.pi.svmp$avg_pi,na.rm=T)
sd(cv1.pi.svmp$avg_pi,na.rm=T)
sd(cv2.pi.svmp$avg_pi,na.rm=T)
sd(co1.pi.svmp$avg_pi,na.rm=T)
sd(co2.pi.svmp$avg_pi,na.rm=T)

mean(cv1co1.dxy.svmp$avg_dxy,na.rm=T)
mean(cv1cv2.dxy.svmp$avg_dxy,na.rm=T)
mean(co1co2.dxy.svmp$avg_dxy,na.rm=T)
sd(cv1co1.dxy.svmp$avg_dxy,na.rm=T)
sd(cv1cv2.dxy.svmp$avg_dxy,na.rm=T)
sd(co1co2.dxy.svmp$avg_dxy,na.rm=T)

mean(cv1co1.fst.svmp$avg_wc_fst,na.rm=T)
mean(cv1cv2.fst.svmp$avg_wc_fst,na.rm=T)
mean(co1co2.fst.svmp$avg_wc_fst,na.rm=T)
sd(cv1co1.fst.svmp$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.svmp$avg_wc_fst,na.rm=T)
sd(co1co2.fst.svmp$avg_wc_fst,na.rm=T)

# SVSP
mean(cv1.pi.svsp$avg_pi,na.rm=T)
mean(cv2.pi.svsp$avg_pi,na.rm=T)
mean(co1.pi.svsp$avg_pi,na.rm=T)
mean(co2.pi.svsp$avg_pi,na.rm=T)
sd(cv1.pi.svsp$avg_pi,na.rm=T)
sd(cv2.pi.svsp$avg_pi,na.rm=T)
sd(co1.pi.svsp$avg_pi,na.rm=T)
sd(co2.pi.svsp$avg_pi,na.rm=T)

mean(cv1co1.dxy.svsp$avg_dxy,na.rm=T)
mean(cv1cv2.dxy.svsp$avg_dxy,na.rm=T)
mean(co1co2.dxy.svsp$avg_dxy,na.rm=T)
sd(cv1co1.dxy.svsp$avg_dxy,na.rm=T)
sd(cv1cv2.dxy.svsp$avg_dxy,na.rm=T)
sd(co1co2.dxy.svsp$avg_dxy,na.rm=T)

mean(cv1co1.fst.svsp$avg_wc_fst,na.rm=T)
mean(cv1cv2.fst.svsp$avg_wc_fst,na.rm=T)
mean(co1co2.fst.svsp$avg_wc_fst,na.rm=T)
sd(cv1co1.fst.svsp$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.svsp$avg_wc_fst,na.rm=T)
sd(co1co2.fst.svsp$avg_wc_fst,na.rm=T)

# PLA2
mean(cv1.pi.pla2$avg_pi,na.rm=T)
mean(cv2.pi.pla2$avg_pi,na.rm=T)
mean(co1.pi.pla2$avg_pi,na.rm=T)
mean(co2.pi.pla2$avg_pi,na.rm=T)
sd(cv1.pi.pla2$avg_pi,na.rm=T)
sd(cv2.pi.pla2$avg_pi,na.rm=T)
sd(co1.pi.pla2$avg_pi,na.rm=T)
sd(co2.pi.pla2$avg_pi,na.rm=T)

mean(cv1co1.dxy.pla2$avg_dxy,na.rm=T)
mean(cv1cv2.dxy.pla2$avg_dxy,na.rm=T)
mean(co1co2.dxy.pla2$avg_dxy,na.rm=T)
sd(cv1co1.dxy.pla2$avg_dxy,na.rm=T)
sd(cv1cv2.dxy.pla2$avg_dxy,na.rm=T)
sd(co1co2.dxy.pla2$avg_dxy,na.rm=T)

mean(cv1co1.fst.pla2$avg_wc_fst,na.rm=T)
mean(cv1cv2.fst.pla2$avg_wc_fst,na.rm=T)
mean(co1co2.fst.pla2$avg_wc_fst,na.rm=T)
sd(cv1co1.fst.pla2$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.pla2$avg_wc_fst,na.rm=T)
sd(co1co2.fst.pla2$avg_wc_fst,na.rm=T)

### 7. Summary statistics for chromosome-specific backgrounds---------------

# Chromosome 9
mean(cv1.pi.mi1.bk$avg_pi,na.rm=T)
mean(cv2.pi.mi1.bk$avg_pi,na.rm=T)
mean(co1.pi.mi1.bk$avg_pi,na.rm=T)
mean(co2.pi.mi1.bk$avg_pi,na.rm=T)
sd(cv1.pi.mi1.bk$avg_pi,na.rm=T)
sd(cv2.pi.mi1.bk$avg_pi,na.rm=T)
sd(co1.pi.mi1.bk$avg_pi,na.rm=T)
sd(co2.pi.mi1.bk$avg_pi,na.rm=T)

mean(cv1co1.dxy.mi1.bk$avg_dxy,na.rm=T)
mean(cv1cv2.dxy.mi1.bk$avg_dxy,na.rm=T)
mean(co1co2.dxy.mi1.bk$avg_dxy,na.rm=T)
sd(cv1co1.dxy.mi1.bk$avg_dxy,na.rm=T)
sd(cv1cv2.dxy.mi1.bk$avg_dxy,na.rm=T)
sd(co1co2.dxy.mi1.bk$avg_dxy,na.rm=T)

mean(cv1co1.fst.mi1.bk$avg_wc_fst,na.rm=T)
mean(cv1cv2.fst.mi1.bk$avg_wc_fst,na.rm=T)
mean(co1co2.fst.mi1.bk$avg_wc_fst,na.rm=T)
sd(cv1co1.fst.mi1.bk$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.mi1.bk$avg_wc_fst,na.rm=T)
sd(co1co2.fst.mi1.bk$avg_wc_fst,na.rm=T)

# Chromosome 10
mean(cv1.pi.mi2.bk$avg_pi,na.rm=T)
mean(cv2.pi.mi2.bk$avg_pi,na.rm=T)
mean(co1.pi.mi2.bk$avg_pi,na.rm=T)
mean(co2.pi.mi2.bk$avg_pi,na.rm=T)
sd(cv1.pi.mi2.bk$avg_pi,na.rm=T)
sd(cv2.pi.mi2.bk$avg_pi,na.rm=T)
sd(co1.pi.mi2.bk$avg_pi,na.rm=T)
sd(co2.pi.mi2.bk$avg_pi,na.rm=T)

mean(cv1co1.dxy.mi2.bk$avg_dxy,na.rm=T)
mean(cv1cv2.dxy.mi2.bk$avg_dxy,na.rm=T)
mean(co1co2.dxy.mi2.bk$avg_dxy,na.rm=T)
sd(cv1co1.dxy.mi2.bk$avg_dxy,na.rm=T)
sd(cv1cv2.dxy.mi2.bk$avg_dxy,na.rm=T)
sd(co1co2.dxy.mi2.bk$avg_dxy,na.rm=T)

mean(cv1co1.fst.mi2.bk$avg_wc_fst,na.rm=T)
mean(cv1cv2.fst.mi2.bk$avg_wc_fst,na.rm=T)
mean(co1co2.fst.mi2.bk$avg_wc_fst,na.rm=T)
sd(cv1co1.fst.mi2.bk$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.mi2.bk$avg_wc_fst,na.rm=T)
sd(co1co2.fst.mi2.bk$avg_wc_fst,na.rm=T)

# Chromosome 15
mean(cv1.pi.mi7.bk$avg_pi,na.rm=T)
mean(cv2.pi.mi7.bk$avg_pi,na.rm=T)
mean(co1.pi.mi7.bk$avg_pi,na.rm=T)
mean(co2.pi.mi7.bk$avg_pi,na.rm=T)
sd(cv1.pi.mi7.bk$avg_pi,na.rm=T)
sd(cv2.pi.mi7.bk$avg_pi,na.rm=T)
sd(co1.pi.mi7.bk$avg_pi,na.rm=T)
sd(co2.pi.mi7.bk$avg_pi,na.rm=T)

mean(cv1co1.dxy.mi7.bk$avg_dxy,na.rm=T)
mean(cv1cv2.dxy.mi7.bk$avg_dxy,na.rm=T)
mean(co1co2.dxy.mi7.bk$avg_dxy,na.rm=T)
sd(cv1co1.dxy.mi7.bk$avg_dxy,na.rm=T)
sd(cv1cv2.dxy.mi7.bk$avg_dxy,na.rm=T)
sd(co1co2.dxy.mi7.bk$avg_dxy,na.rm=T)

mean(cv1co1.fst.mi7.bk$avg_wc_fst,na.rm=T)
mean(cv1cv2.fst.mi7.bk$avg_wc_fst,na.rm=T)
mean(co1co2.fst.mi7.bk$avg_wc_fst,na.rm=T)
sd(cv1co1.fst.mi7.bk$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.mi7.bk$avg_wc_fst,na.rm=T)
sd(co1co2.fst.mi7.bk$avg_wc_fst,na.rm=T)

### 7.1 Summary statistics for non-venom homologs---------------------------

mean(cv1.pi.nv.svmp$avg_pi)
sd(cv1.pi.nv.svmp$avg_pi)
mean(cv2.pi.nv.svmp$avg_pi)
sd(cv2.pi.nv.svmp$avg_pi)
mean(co1.pi.nv.svmp$avg_pi)
sd(co1.pi.nv.svmp$avg_pi)
mean(co2.pi.nv.svmp$avg_pi)
sd(co2.pi.nv.svmp$avg_pi)

mean(cv1.pi.nv.svsp$avg_pi)
sd(cv1.pi.nv.svsp$avg_pi)
mean(cv2.pi.nv.svsp$avg_pi)
sd(cv2.pi.nv.svsp$avg_pi)
mean(co1.pi.nv.svsp$avg_pi)
sd(co1.pi.nv.svsp$avg_pi)
mean(co2.pi.nv.svsp$avg_pi)
sd(co2.pi.nv.svsp$avg_pi)

mean(cv1.pi.nv.pla2$avg_pi)
sd(cv1.pi.nv.pla2$avg_pi)
mean(cv2.pi.nv.pla2$avg_pi)
sd(cv2.pi.nv.pla2$avg_pi)
mean(co1.pi.nv.pla2$avg_pi)
sd(co1.pi.nv.pla2$avg_pi)
mean(co2.pi.nv.pla2$avg_pi)
sd(co2.pi.nv.pla2$avg_pi)

mean(cv1co1.dxy.nv.svmp$avg_dxy)
sd(cv1co1.dxy.nv.svmp$avg_dxy)
mean(cv1cv2.dxy.nv.svmp$avg_dxy)
sd(cv1cv2.dxy.nv.svmp$avg_dxy)
mean(co1co2.dxy.nv.svmp$avg_dxy)
sd(co1co2.dxy.nv.svmp$avg_dxy)

mean(cv1co1.dxy.nv.svsp$avg_dxy)
sd(cv1co1.dxy.nv.svsp$avg_dxy)
mean(cv1cv2.dxy.nv.svsp$avg_dxy)
sd(cv1cv2.dxy.nv.svsp$avg_dxy)
mean(co1co2.dxy.nv.svsp$avg_dxy)
sd(co1co2.dxy.nv.svsp$avg_dxy)

mean(cv1co1.dxy.nv.pla2$avg_dxy)
sd(cv1co1.dxy.nv.pla2$avg_dxy)
mean(cv1cv2.dxy.nv.pla2$avg_dxy)
sd(cv1cv2.dxy.nv.pla2$avg_dxy)
mean(co1co2.dxy.nv.pla2$avg_dxy)
sd(co1co2.dxy.nv.pla2$avg_dxy)

mean(cv1co1.fst.nv.svmp$avg_wc_fst)
sd(cv1co1.fst.nv.svmp$avg_wc_fst)
mean(cv1cv2.fst.nv.svmp$avg_wc_fst)
sd(cv1cv2.fst.nv.svmp$avg_wc_fst)
mean(co1co2.fst.nv.svmp$avg_wc_fst)
sd(co1co2.fst.nv.svmp$avg_wc_fst)

mean(cv1co1.fst.nv.svsp$avg_wc_fst)
sd(cv1co1.fst.nv.svsp$avg_wc_fst)
mean(cv1cv2.fst.nv.svsp$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.nv.svsp$avg_wc_fst,na.rm=T)
mean(co1co2.fst.nv.svsp$avg_wc_fst)
sd(co1co2.fst.nv.svsp$avg_wc_fst)

mean(cv1co1.fst.nv.pla2$avg_wc_fst)
sd(cv1co1.fst.nv.pla2$avg_wc_fst)
mean(cv1cv2.fst.nv.pla2$avg_wc_fst,na.rm=T)
sd(cv1cv2.fst.nv.pla2$avg_wc_fst,na.rm=T)
mean(co1co2.fst.nv.pla2$avg_wc_fst)
sd(co1co2.fst.nv.pla2$avg_wc_fst)

### 8. Mann-Whitney U tests comparing venom regions to backgrounds----------

# Pi and dxy distributions are definitely not normally-distributed, so we'll
# avoid two-sample t-tests, and use wilcox tests instead.

# SVMP/chromosome 9
wilcox.test(cv1.pi.svmp$avg_pi,cv1.pi.mi1.bk$avg_pi)
wilcox.test(cv2.pi.svmp$avg_pi,cv2.pi.mi1.bk$avg_pi)
wilcox.test(co1.pi.svmp$avg_pi,co1.pi.mi1.bk$avg_pi)
wilcox.test(co2.pi.svmp$avg_pi,co2.pi.mi1.bk$avg_pi)

wilcox.test(cv1co1.dxy.svmp$avg_dxy,cv1co1.dxy.mi1.bk$avg_dxy)
wilcox.test(cv1cv2.dxy.svmp$avg_dxy,cv1cv2.dxy.mi1.bk$avg_dxy)
wilcox.test(co1co2.dxy.svmp$avg_dxy,co1co2.dxy.mi1.bk$avg_dxy)

wilcox.test(cv1co1.fst.svmp$avg_wc_fst,cv1co1.fst.mi1.bk$avg_wc_fst)
wilcox.test(cv1cv2.fst.svmp$avg_wc_fst,cv1cv2.fst.mi1.bk$avg_wc_fst)
wilcox.test(co1co2.fst.svmp$avg_wc_fst,co1co2.fst.mi1.bk$avg_wc_fst)

wilcox.test(cv1.pi.svmp$avg_pi,cv1.pi.nv.svmp$avg_pi)
wilcox.test(cv2.pi.svmp$avg_pi,cv2.pi.nv.svmp$avg_pi)
wilcox.test(co1.pi.svmp$avg_pi,co1.pi.nv.svmp$avg_pi)
wilcox.test(co2.pi.svmp$avg_pi,co2.pi.nv.svmp$avg_pi)

wilcox.test(cv1co1.dxy.svmp$avg_dxy,cv1co1.dxy.nv.svmp$avg_dxy)
wilcox.test(cv1cv2.dxy.svmp$avg_dxy,cv1cv2.dxy.nv.svmp$avg_dxy)
wilcox.test(co1co2.dxy.svmp$avg_dxy,co1co2.dxy.nv.svmp$avg_dxy)

wilcox.test(cv1co1.fst.svmp$avg_wc_fst,cv1co1.fst.nv.svmp$avg_wc_fst)
wilcox.test(cv1cv2.fst.svmp$avg_wc_fst,cv1cv2.fst.nv.svmp$avg_wc_fst)
wilcox.test(co1co2.fst.svmp$avg_wc_fst,co1co2.fst.nv.svmp$avg_wc_fst)

# SVSP/chromosome 10
wilcox.test(cv1.pi.svsp$avg_pi,cv1.pi.mi2.bk$avg_pi)
wilcox.test(cv2.pi.svsp$avg_pi,cv2.pi.mi2.bk$avg_pi)
wilcox.test(co1.pi.svsp$avg_pi,co1.pi.mi2.bk$avg_pi)
wilcox.test(co2.pi.svsp$avg_pi,co2.pi.mi2.bk$avg_pi)

wilcox.test(cv1co1.dxy.svsp$avg_dxy,cv1co1.dxy.mi2.bk$avg_dxy)
wilcox.test(cv1cv2.dxy.svsp$avg_dxy,cv1cv2.dxy.mi2.bk$avg_dxy)
wilcox.test(co1co2.dxy.svsp$avg_dxy,co1co2.dxy.mi2.bk$avg_dxy)

wilcox.test(cv1co1.fst.svsp$avg_wc_fst,cv1co1.fst.mi2.bk$avg_wc_fst)
wilcox.test(cv1cv2.fst.svsp$avg_wc_fst,cv1cv2.fst.mi2.bk$avg_wc_fst)
wilcox.test(co1co2.fst.svsp$avg_wc_fst,co1co2.fst.mi2.bk$avg_wc_fst)

wilcox.test(cv1.pi.svsp$avg_pi,cv1.pi.nv.svsp$avg_pi)
wilcox.test(cv2.pi.svsp$avg_pi,cv2.pi.nv.svsp$avg_pi)
wilcox.test(co1.pi.svsp$avg_pi,co1.pi.nv.svsp$avg_pi)
wilcox.test(co2.pi.svsp$avg_pi,co2.pi.nv.svsp$avg_pi)

wilcox.test(cv1co1.dxy.svsp$avg_dxy,cv1co1.dxy.nv.svsp$avg_dxy)
wilcox.test(cv1cv2.dxy.svsp$avg_dxy,cv1cv2.dxy.nv.svsp$avg_dxy)
wilcox.test(co1co2.dxy.svsp$avg_dxy,co1co2.dxy.nv.svsp$avg_dxy)

wilcox.test(cv1co1.fst.svsp$avg_wc_fst,cv1co1.fst.nv.svsp$avg_wc_fst)
wilcox.test(cv1cv2.fst.svsp$avg_wc_fst,cv1cv2.fst.nv.svsp$avg_wc_fst)
wilcox.test(co1co2.fst.svsp$avg_wc_fst,co1co2.fst.nv.svsp$avg_wc_fst)

# PLA2/chromosome 15
wilcox.test(cv1.pi.pla2$avg_pi,cv1.pi.mi7.bk$avg_pi)
wilcox.test(cv2.pi.pla2$avg_pi,cv2.pi.mi7.bk$avg_pi)
wilcox.test(co1.pi.pla2$avg_pi,co1.pi.mi7.bk$avg_pi)
wilcox.test(co2.pi.pla2$avg_pi,co2.pi.mi7.bk$avg_pi)

wilcox.test(cv1co1.dxy.pla2$avg_dxy,cv1co1.dxy.mi7.bk$avg_dxy)
wilcox.test(cv1cv2.dxy.pla2$avg_dxy,cv1cv2.dxy.mi7.bk$avg_dxy)
wilcox.test(co1co2.dxy.pla2$avg_dxy,co1co2.dxy.mi7.bk$avg_dxy)

wilcox.test(cv1co1.fst.pla2$avg_wc_fst,cv1co1.fst.mi7.bk$avg_wc_fst)
wilcox.test(cv1cv2.fst.pla2$avg_wc_fst,cv1cv2.fst.mi7.bk$avg_wc_fst)
wilcox.test(co1co2.fst.pla2$avg_wc_fst,co1co2.fst.mi7.bk$avg_wc_fst)

wilcox.test(cv1.pi.pla2$avg_pi,cv1.pi.nv.pla2$avg_pi)
wilcox.test(cv2.pi.pla2$avg_pi,cv2.pi.nv.pla2$avg_pi)
wilcox.test(co1.pi.pla2$avg_pi,co1.pi.nv.pla2$avg_pi)
wilcox.test(co2.pi.pla2$avg_pi,co2.pi.nv.pla2$avg_pi)

wilcox.test(cv1co1.dxy.pla2$avg_dxy,cv1co1.dxy.nv.pla2$avg_dxy)
wilcox.test(cv1cv2.dxy.pla2$avg_dxy,cv1cv2.dxy.nv.pla2$avg_dxy)
wilcox.test(co1co2.dxy.pla2$avg_dxy,co1co2.dxy.nv.pla2$avg_dxy)

wilcox.test(cv1co1.fst.pla2$avg_wc_fst,cv1co1.fst.nv.pla2$avg_wc_fst)
wilcox.test(cv1cv2.fst.pla2$avg_wc_fst,cv1cv2.fst.nv.pla2$avg_wc_fst)
wilcox.test(co1co2.fst.pla2$avg_wc_fst,co1co2.fst.nv.pla2$avg_wc_fst)

### 9. Boxplots comparing gene-means to chromosome backgrounds--------------

# Make vectors of gene means
cv1.pi.svmp.m <- c(0.0056,0.0048,0.0058,0.0043,0.0045,0.0051,0.0034,0.0028,0.0028,0.0032,0.0024)
cv2.pi.svmp.m <- c(0.0041,0.0036,0.0042,0.0038,0.0041,0.0046,0.0036,0.0041,0.0034,0.0034,0.0018)
co1.pi.svmp.m <- c(0.0071,0.0113,0.0112,0.0038)
co2.pi.svmp.m <- c(0.0077,0.0131,0.0114,0.0043)

cv1.pi.svsp.m <- c(0.0076,0.0071,0.0049,0.0053,0.0168,0.0126,0.0101,0.0030,0.0046,0.0078,0.0042)
cv2.pi.svsp.m <- c(0.0038,0.0044,0.0030,0.0025,0.0144,0.0089,0.0061,0.0029,0.0035,0.0079,0.0037)
co1.pi.svsp.m <- c(0.0099,0.0055,0.0063,0.0068,0.0040,0.0054,0.0064)
co2.pi.svsp.m <- c(0.0126,0.0068,0.0106,0.0083,0.0064,0.0074,0.0103)

cv1.pi.pla2.m <- c(0.0024,0.0008,0.0053,0.0133)
cv2.pi.pla2.m <- c(0.0016,0.0009,0.0064,0.0174)
co1.pi.pla2.m <- c(0.0121)
co2.pi.pla2.m <- c(0.0024)

cv1co1.dxy.svmp.m <- c(0.0087,0.0172,0.0305,0.0060)
cv1cv2.dxy.svmp.m <- c(0.0052,0.0043,0.0053,0.0041,0.0046,0.0051,0.0036,0.0039,0.0036,0.0036,0.0022)
co1co2.dxy.svmp.m <- c(0.0079,0.0123,0.0113,0.0050)

cv1co1.dxy.svsp.m <- c(0.0157,0.0087,0.0137,0.0125,0.0052,0.0097,0.0122)
cv1cv2.dxy.svsp.m <- c(0.0059,0.0059,0.0039,0.0040,0.0153,0.0109,0.0082,0.0030,0.0041,0.0078,0.0039)
co1co2.dxy.svsp.m <- c(0.0111,0.0062,0.0086,0.0077,0.0058,0.0077,0.0099)

cv1co1.dxy.pla2.m <- c(0.0241)
cv1cv2.dxy.pla2.m <- c(0.0030,0.0009,0.0061,0.0165)
co1co2.dxy.pla2.m <- c(0.0156)

cv1co1.fst.svmp.m <- c(0.3065,0.6018,0.5793,0.5239)
cv1cv2.fst.svmp.m <- c(0.0881,0.0640,0.1024,0.0720,0.0668,0.0546,0.0533,0.0966,0.1744,0.0813,0.0448)
co1co2.fst.svmp.m <- c(0.0905,0.0606,0.0652,0.1963)

cv1co1.fst.svsp.m <- c(0.3705,0.3031,0.2677,0.2982,0.2943,0.5704,0.4209)
cv1cv2.fst.svsp.m <- c(0.0374,0.0512,0.0917,0.1006,0.0025,0.0387,0.0344,0.0019,0.0455,0.0141,0.0098)
co1co2.fst.svsp.m <- c(0.0192,0.0238,0.1026,0.0624,0.1766,0.2873,0.1708)

cv1co1.fst.pla2.m <- c(0.5316)
cv1cv2.fst.pla2.m <- c(0.3783,0.0407,0.0545,0.0831)
co1co2.fst.pla2.m <- c(0.5246)

# Plot boxplots with gene-specific points
par(mfrow=c(3,3))
boxplot(cv1.pi.mi1.bk$avg_pi,cv2.pi.mi1.bk$avg_pi,co1.pi.mi1.bk$avg_pi,co2.pi.mi1.bk$avg_pi,col='lightblue3',outline=F,ylim=c(0,0.015),ylab='pi',names=c('CV1','CV2','CO1','CO2'))
for (p in cv1.pi.svmp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv2.pi.svmp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.pi.svmp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co2.pi.svmp.m){
  points(4,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1.pi.mi2.bk$avg_pi,cv2.pi.mi2.bk$avg_pi,co1.pi.mi2.bk$avg_pi,co2.pi.mi2.bk$avg_pi,col='aquamarine3',outline=F,ylim=c(0,0.015),ylab='pi',names=c('CV1','CV2','CO1','CO2'))
for (p in cv1.pi.svsp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv2.pi.svsp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.pi.svsp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co2.pi.svsp.m){
  points(4,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1.pi.mi7.bk$avg_pi,cv2.pi.mi7.bk$avg_pi,co1.pi.mi7.bk$avg_pi,co2.pi.mi7.bk$avg_pi,col='maroon',outline=F,ylim=c(0,0.015),ylab='pi',names=c('CV1','CV2','CO1','CO2'))
for (p in cv1.pi.pla2.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv2.pi.pla2.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1.pi.pla2.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co2.pi.pla2.m){
  points(4,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.dxy.mi1.bk$avg_dxy,cv1cv2.dxy.mi1.bk$avg_dxy,co1co2.dxy.mi1.bk$avg_dxy,col='lightblue3',outline=F,ylim=c(0,0.035),ylab='dxy',names=c('CV1-CO1','CV1-CV2','CO1-CO2'))
for (p in cv1co1.dxy.svmp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv1cv2.dxy.svmp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1co2.dxy.svmp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.dxy.mi2.bk$avg_dxy,cv1cv2.dxy.mi2.bk$avg_dxy,co1co2.dxy.mi2.bk$avg_dxy,col='aquamarine3',outline=F,ylim=c(0,0.035),ylab='dxy',names=c('CV1-CO1','CV1-CV2','CO1-CO2'))
for (p in cv1co1.dxy.svsp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv1cv2.dxy.svsp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1co2.dxy.svsp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.dxy.mi7.bk$avg_dxy,cv1cv2.dxy.mi7.bk$avg_dxy,co1co2.dxy.mi7.bk$avg_dxy,col='maroon',outline=F,ylim=c(0,0.035),ylab='dxy',names=c('CV1-CO1','CV1-CV2','CO1-CO2'))
for (p in cv1co1.dxy.pla2.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv1cv2.dxy.pla2.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1co2.dxy.pla2.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.fst.mi1.bk$avg_wc_fst,cv1cv2.fst.mi1.bk$avg_wc_fst,co1co2.fst.mi1.bk$avg_wc_fst,col='lightblue3',outline=F,ylab='fst',names=c('CV1-CO1','CV1-CV2','CO1-CO2'))
for (p in cv1co1.fst.svmp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv1cv2.fst.svmp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1co2.fst.svmp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.fst.mi2.bk$avg_wc_fst,cv1cv2.fst.mi2.bk$avg_wc_fst,co1co2.fst.mi2.bk$avg_wc_fst,col='aquamarine3',outline=F,ylab='fst',names=c('CV1-CO1','CV1-CV2','CO1-CO2'))
for (p in cv1co1.fst.svsp.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv1cv2.fst.svsp.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1co2.fst.svsp.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 

boxplot(cv1co1.fst.mi7.bk$avg_wc_fst,cv1cv2.fst.mi7.bk$avg_wc_fst,co1co2.fst.mi7.bk$avg_wc_fst,col='maroon',outline=F,ylab='fst',names=c('CV1-CO1','CV1-CV2','CO1-CO2'))
for (p in cv1co1.fst.pla2.m){
  points(1,p,pch=24,lwd=1.5,col='navy')
} 
for (p in cv1cv2.fst.pla2.m){
  points(2,p,pch=24,lwd=1.5,col='navy')
} 
for (p in co1co2.fst.pla2.m){
  points(3,p,pch=24,lwd=1.5,col='navy')
} 
