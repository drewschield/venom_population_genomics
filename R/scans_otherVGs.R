############################################################################
# Genomic scans over 'other' venom genes
############################################################################

### Goal: study variation in population genetic parameters across venom genes
### other than the 3 major venom gene tandem arrays. The 'other' genes are
### those with significant gene expression in the C. viridis venom gland. These
### include CRISPs, C-type lectin, LAAOs, and vespryn.

### 'Other' venom gene genomic coordinates:

# C-type lectin
# scaffold-mi5	11650747	11653723	C-type lectin	maker-scaffold-mi5-augustus-gene-38.4	crovir-transcript-564	CTL

# CRISPs
# scaffold-ma1	169423774	169434684	Cysteine-rich secretory protein	maker-scaffold-ma1-augustus-gene-564.1	crovir-transcript-8573	CRISP_1
# scaffold-ma1	169434958	169437996	Cysteine-rich secretory protein	maker-scaffold-ma1-augustus-gene-564.2	crovir-transcript-8572	CRISP_2

# Exonuclease
# scaffold-ma6	12590208	12591465	Exonuclease	maker-scaffold-ma6-augustus-gene-42.10	crovir-transcript-11720	Exonuclease_1
# scaffold-mi3	10271502	10274220	Exonuclease	maker-scaffold-mi3-augustus-gene-34.12	crovir-transcript-11292	Exonuclease_2
# scaffold-mi7	8097114	8103411	Exonuclease	maker-scaffold-mi7-augustus-gene-27.16	crovir-transcript-1466	Exonuclease_3

# Glutaminyl cyclase
# scaffold-ma1	256551622	256564040	Glutaminyl cyclase	maker-scaffold-ma1-augustus-gene-855.11	crovir-transcript-9447	vQC_1
# scaffold-mi7	5091107	5094268	Glutaminyl cyclase	maker-scaffold-mi7-augustus-gene-17.22	crovir-transcript-1322	vQC_2

# LAAOs
# scaffold-ma2	4654769	4658293	L-amino acid oxidase	maker-scaffold-ma2-augustus-gene-15.17	crovir-transcript-13249	LAAO_1
# scaffold-ma2	4658599	4661642	L-amino acid oxidase	maker-scaffold-ma2-augustus-gene-15.18	crovir-transcript-13246	LAAO_2
# scaffold-ma4	85461961	85468906	L-amino acid oxidase	maker-scaffold-ma4-augustus-gene-284.10	crovir-transcript-3893	LAAO_3

# VEGFs
# scaffold-ma1	260248287	260272500	Vascular endothelial growth factor	maker-scaffold-ma1-augustus-gene-867.2	crovir-transcript-9439	VEGF_1
# scaffold-ma7	40288572	40327884	Vascular endothelial growth factor	maker-scaffold-ma7-augustus-gene-134.2	crovir-transcript-4609	VEGF_2

# Vespryn
# scaffold-ma2	4377779	4385668	Vespryn/Ohanin	maker-scaffold-ma2-augustus-gene-14.15	crovir-transcript-13215	Vespryn

### Clear environment-------------------------------------------------------

rm(list = ls())

### Load dependencies-------------------------------------------------------

library(scales)
library(data.table)
library(dplyr)

options('stringsAsFactors'=FALSE)

### C-type lectin-----------------------------------------------------------

# Coordinates:
# scaffold-mi5	11650747	11653723

# Read in data
ctl.1kb.pi.cv <- read.table('./pixy/pixy.scaffold-mi5.1kb_pi.txt',header=T)
ctl.1kb.pi.cv <- ctl.1kb.pi.cv[which(ctl.1kb.pi.cv$pop=='CV1'),]

ctl.1kb.pi.co <- read.table('./pixy/pixy.scaffold-mi5.1kb_pi.txt',header=T)
ctl.1kb.pi.co <- ctl.1kb.pi.co[which(ctl.1kb.pi.co$pop=='CO1'),]

ctl.1kb.dxy <- read.table('./pixy/pixy.scaffold-mi5.1kb_dxy.txt',header=T)
ctl.1kb.dxy <- ctl.1kb.dxy[which(ctl.1kb.dxy$pop1=='CV1' & ctl.1kb.dxy$pop2=='CO1'),]

ctl.1kb.fst <- read.table('./pixy/pixy.scaffold-mi5.1kb_fst.txt',header=T)
ctl.1kb.fst <- ctl.1kb.fst[which(ctl.1kb.fst$pop1=='CV1' & ctl.1kb.fst$pop2=='CO1'),]

ctl.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
ctl.1kb.taj.cv <- ctl.1kb.taj.cv[which(ctl.1kb.taj.cv$CHROM=='scaffold-mi5'),]

ctl.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
ctl.1kb.taj.co <- ctl.1kb.taj.co[which(ctl.1kb.taj.co$CHROM=='scaffold-mi5'),]

ctl.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-mi5.cv1.co1.txt',header=T)

ctl.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
ctl.1kb.ihs.cv <- ctl.1kb.ihs.cv[which(ctl.1kb.ihs.cv$chrom=='scaffold-mi5'),]

ctl.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
ctl.1kb.ihs.co <- ctl.1kb.ihs.co[which(ctl.1kb.ihs.co$chrom=='scaffold-mi5'),]

ctl.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
ctl.1kb.beta.cv <- ctl.1kb.beta.cv[which(ctl.1kb.beta.cv$chrom=='scaffold-mi5'),]

ctl.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
ctl.1kb.beta.co <- ctl.1kb.beta.co[which(ctl.1kb.beta.co$chrom=='scaffold-mi5'),]

# Plot scans for CTL
par(mfrow=c(7,1))
plot(ctl.1kb.pi.cv$window_pos_1,ctl.1kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(11600000,11700000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')
plot(ctl.1kb.pi.co$window_pos_1,ctl.1kb.pi.co$avg_pi,type='l',col='navy',xlim=c(11600000,11700000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')

plot(ctl.1kb.dxy$window_pos_1,ctl.1kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(11600000,11700000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')

plot(ctl.1kb.fst$window_pos_1,ctl.1kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(11600000,11700000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')

plot(ctl.1kb.taj.cv$BIN_START,ctl.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(11600000,11700000),ylim=c(-2,2),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')
plot(ctl.1kb.taj.co$BIN_START,ctl.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(11600000,11700000),ylim=c(-2,2),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')

plot(ctl.1kb.df$start,ctl.1kb.df$df,type='l',col='darkgreen',xlim=c(11600000,11700000),ylim=c(-0.05,0.2),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')

plot(ctl.1kb.ihs.cv$start,ctl.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(11600000,11700000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')
plot(ctl.1kb.ihs.co$start,ctl.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(11600000,11700000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')

plot(ctl.1kb.beta.cv$start,ctl.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(11600000,11700000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')
plot(ctl.1kb.beta.co$start,ctl.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(11600000,11700000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=11650747,y0=vert,x1=11653723,y1=vert,col='darkblue')

### CRISPs------------------------------------------------------------------

# Coordinates:
# scaffold-ma1	169423774	169434684	Cysteine-rich secretory protein	maker-scaffold-ma1-augustus-gene-564.1	crovir-transcript-8573	CRISP_1
# scaffold-ma1	169434958	169437996	Cysteine-rich secretory protein	maker-scaffold-ma1-augustus-gene-564.2	crovir-transcript-8572	CRISP_2

# Read in data
csp.10kb.pi.cv <- read.table('./pixy/pi.10kb.cv1.txt',header=T)
csp.10kb.pi.cv <- csp.10kb.pi.cv[which(csp.10kb.pi.cv$chromosome=='scaffold-ma1'),]

csp.10kb.pi.co <- read.table('./pixy/pi.10kb.co1.txt',header=T)
csp.10kb.pi.co <- csp.10kb.pi.co[which(csp.10kb.pi.co$chromosome=='scaffold-ma1'),]

csp.10kb.dxy <- read.table('./pixy/dxy.10kb.cv1co1.txt',header=T)
csp.10kb.dxy <- csp.10kb.dxy[which(csp.10kb.dxy$chromosome=='scaffold-ma1'),]

csp.10kb.fst <- read.table('./pixy/fst.10kb.cv1co1.txt',header=T)
csp.10kb.fst <- csp.10kb.fst[which(csp.10kb.fst$chromosome=='scaffold-ma1'),]

csp.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
csp.1kb.taj.cv <- csp.1kb.taj.cv[which(csp.1kb.taj.cv$CHROM=='scaffold-ma1'),]

csp.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
csp.1kb.taj.co <- csp.1kb.taj.co[which(csp.1kb.taj.co$CHROM=='scaffold-ma1'),]

csp.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma1.cv1.co1.txt',header=T)

csp.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
csp.1kb.ihs.cv <- csp.1kb.ihs.cv[which(csp.1kb.ihs.cv$chrom=='scaffold-ma1'),]

csp.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
csp.1kb.ihs.co <- csp.1kb.ihs.co[which(csp.1kb.ihs.co$chrom=='scaffold-ma1'),]

csp.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
csp.1kb.beta.cv <- csp.1kb.beta.cv[which(csp.1kb.beta.cv$chrom=='scaffold-ma1'),]

csp.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
csp.1kb.beta.co <- csp.1kb.beta.co[which(csp.1kb.beta.co$chrom=='scaffold-ma1'),]

# Plot scans for CRISPs

# scaffold-ma1	169423774	169434684	Cysteine-rich secretory protein	maker-scaffold-ma1-augustus-gene-564.1	crovir-transcript-8573	CRISP_1
# scaffold-ma1	169434958	169437996	Cysteine-rich secretory protein	maker-scaffold-ma1-augustus-gene-564.2	crovir-transcript-8572	CRISP_2

par(mfrow=c(7,1))
plot(csp.10kb.pi.cv$window_pos_1,csp.10kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(169390000,169465000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')
plot(csp.10kb.pi.co$window_pos_1,csp.10kb.pi.co$avg_pi,type='l',col='navy',xlim=c(169390000,169465000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')

plot(csp.10kb.dxy$window_pos_1,csp.10kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(169390000,169465000),ylim=c(-0.004,0.03),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')

plot(csp.10kb.fst$window_pos_1,csp.10kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(169390000,169465000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')

plot(csp.1kb.taj.cv$BIN_START,csp.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(169390000,169465000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')
plot(csp.1kb.taj.co$BIN_START,csp.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(169390000,169465000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')

plot(csp.1kb.df$start,csp.1kb.df$df,type='l',col='darkgreen',xlim=c(169390000,169465000),ylim=c(-0.05,0.4),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')

plot(csp.1kb.ihs.cv$start,csp.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(169390000,169465000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')
plot(csp.1kb.ihs.co$start,csp.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(169390000,169465000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')

plot(csp.1kb.beta.cv$start,csp.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(169390000,169465000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')
plot(csp.1kb.beta.co$start,csp.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(169390000,169465000),ylim=c(-3,20),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=169423774,y0=vert,x1=169434684,y1=vert,col='darkblue')
segments(x0=169434958,y0=vert,x1=169437996,y1=vert,col='darkblue')

### Exonuclease-------------------------------------------------------------

# Coordinates:
# scaffold-ma6	12590208	12591465	Exonuclease	maker-scaffold-ma6-augustus-gene-42.10	crovir-transcript-11720	Exonuclease_1
# scaffold-mi3	10271502	10274220	Exonuclease	maker-scaffold-mi3-augustus-gene-34.12	crovir-transcript-11292	Exonuclease_2
# scaffold-mi7	8097114	8103411	Exonuclease	maker-scaffold-mi7-augustus-gene-27.16	crovir-transcript-1466	Exonuclease_3

# Read in data for exonuclease 1
exo1.1kb.pi.cv <- read.table('./pixy/pixy.scaffold-ma6.1kb_pi.txt',header=T)
exo1.1kb.pi.cv <- exo1.1kb.pi.cv[which(exo1.1kb.pi.cv$pop=='CV1'),]

exo1.1kb.pi.co <- read.table('./pixy/pixy.scaffold-ma6.1kb_pi.txt',header=T)
exo1.1kb.pi.co <- exo1.1kb.pi.co[which(exo1.1kb.pi.co$pop=='CO1'),]

exo1.1kb.dxy <- read.table('./pixy/pixy.scaffold-ma6.1kb_dxy.txt',header=T)
exo1.1kb.dxy <- exo1.1kb.dxy[which(exo1.1kb.dxy$pop1=='CV1' & exo1.1kb.dxy$pop2=='CO1'),]

exo1.1kb.fst <- read.table('./pixy/pixy.scaffold-ma6.1kb_fst.txt',header=T)
exo1.1kb.fst <- exo1.1kb.fst[which(exo1.1kb.fst$pop1=='CV1' & exo1.1kb.fst$pop2=='CO1'),]

exo1.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
exo1.1kb.taj.cv <- exo1.1kb.taj.cv[which(exo1.1kb.taj.cv$CHROM=='scaffold-ma6'),]

exo1.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
exo1.1kb.taj.co <- exo1.1kb.taj.co[which(exo1.1kb.taj.co$CHROM=='scaffold-ma6'),]

exo1.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma6.cv1.co1.txt',header=T)

exo1.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
exo1.1kb.ihs.cv <- exo1.1kb.ihs.cv[which(exo1.1kb.ihs.cv$chrom=='scaffold-ma6'),]

exo1.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
exo1.1kb.ihs.co <- exo1.1kb.ihs.co[which(exo1.1kb.ihs.co$chrom=='scaffold-ma6'),]

exo1.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
exo1.1kb.beta.cv <- exo1.1kb.beta.cv[which(exo1.1kb.beta.cv$chrom=='scaffold-ma6'),]

exo1.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
exo1.1kb.beta.co <- exo1.1kb.beta.co[which(exo1.1kb.beta.co$chrom=='scaffold-ma6'),]

# Read in data for exonuclease 2
exo2.1kb.pi.cv <- read.table('./pixy/pixy.scaffold-mi3.1kb_pi.txt',header=T)
exo2.1kb.pi.cv <- exo2.1kb.pi.cv[which(exo2.1kb.pi.cv$pop=='CV1'),]

exo2.1kb.pi.co <- read.table('./pixy/pixy.scaffold-mi3.1kb_pi.txt',header=T)
exo2.1kb.pi.co <- exo2.1kb.pi.co[which(exo2.1kb.pi.co$pop=='CO1'),]

exo2.1kb.dxy <- read.table('./pixy/pixy.scaffold-mi3.1kb_dxy.txt',header=T)
exo2.1kb.dxy <- exo2.1kb.dxy[which(exo2.1kb.dxy$pop1=='CV1' & exo2.1kb.dxy$pop2=='CO1'),]

exo2.1kb.fst <- read.table('./pixy/pixy.scaffold-mi3.1kb_fst.txt',header=T)
exo2.1kb.fst <- exo2.1kb.fst[which(exo2.1kb.fst$pop1=='CV1' & exo2.1kb.fst$pop2=='CO1'),]

exo2.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
exo2.1kb.taj.cv <- exo2.1kb.taj.cv[which(exo2.1kb.taj.cv$CHROM=='scaffold-mi3'),]

exo2.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
exo2.1kb.taj.co <- exo2.1kb.taj.co[which(exo2.1kb.taj.co$CHROM=='scaffold-mi3'),]

exo2.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-mi3.cv1.co1.txt',header=T)

exo2.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
exo2.1kb.ihs.cv <- exo2.1kb.ihs.cv[which(exo2.1kb.ihs.cv$chrom=='scaffold-mi3'),]

exo2.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
exo2.1kb.ihs.co <- exo2.1kb.ihs.co[which(exo2.1kb.ihs.co$chrom=='scaffold-mi3'),]

exo2.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
exo2.1kb.beta.cv <- exo2.1kb.beta.cv[which(exo2.1kb.beta.cv$chrom=='scaffold-mi3'),]

exo2.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
exo2.1kb.beta.co <- exo2.1kb.beta.co[which(exo2.1kb.beta.co$chrom=='scaffold-mi3'),]

# Read in data for exonuclease 3
exo3.1kb.pi.cv <- read.table('./pixy/pixy.scaffold-mi7.1kb_pi.txt',header=T)
exo3.1kb.pi.cv <- exo3.1kb.pi.cv[which(exo3.1kb.pi.cv$pop=='CV1'),]

exo3.1kb.pi.co <- read.table('./pixy/pixy.scaffold-mi7.1kb_pi.txt',header=T)
exo3.1kb.pi.co <- exo3.1kb.pi.co[which(exo3.1kb.pi.co$pop=='CO1'),]

exo3.1kb.dxy <- read.table('./pixy/pixy.scaffold-mi7.1kb_dxy.txt',header=T)
exo3.1kb.dxy <- exo3.1kb.dxy[which(exo3.1kb.dxy$pop1=='CV1' & exo3.1kb.dxy$pop2=='CO1'),]

exo3.1kb.fst <- read.table('./pixy/pixy.scaffold-mi7.1kb_fst.txt',header=T)
exo3.1kb.fst <- exo3.1kb.fst[which(exo3.1kb.fst$pop1=='CV1' & exo3.1kb.fst$pop2=='CO1'),]

exo3.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
exo3.1kb.taj.cv <- exo3.1kb.taj.cv[which(exo3.1kb.taj.cv$CHROM=='scaffold-mi7'),]

exo3.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
exo3.1kb.taj.co <- exo3.1kb.taj.co[which(exo3.1kb.taj.co$CHROM=='scaffold-mi7'),]

exo3.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-mi7.cv1.co1.txt',header=T)

exo3.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
exo3.1kb.ihs.cv <- exo3.1kb.ihs.cv[which(exo3.1kb.ihs.cv$chrom=='scaffold-mi7'),]

exo3.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
exo3.1kb.ihs.co <- exo3.1kb.ihs.co[which(exo3.1kb.ihs.co$chrom=='scaffold-mi7'),]

exo3.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
exo3.1kb.beta.cv <- exo3.1kb.beta.cv[which(exo3.1kb.beta.cv$chrom=='scaffold-mi7'),]

exo3.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
exo3.1kb.beta.co <- exo3.1kb.beta.co[which(exo3.1kb.beta.co$chrom=='scaffold-mi7'),]

# Plot scans for Exonuclease 1

par(mfrow=c(7,1))
plot(exo1.1kb.pi.cv$window_pos_1,exo1.1kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(12560000,12620000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')
plot(exo1.1kb.pi.co$window_pos_1,exo1.1kb.pi.co$avg_pi,type='l',col='navy',xlim=c(12560000,12620000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')

plot(exo1.1kb.dxy$window_pos_1,exo1.1kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(12560000,12620000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')

plot(exo1.1kb.fst$window_pos_1,exo1.1kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(12560000,12620000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')

plot(exo1.1kb.taj.cv$BIN_START,exo1.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(12560000,12620000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')
plot(exo1.1kb.taj.co$BIN_START,exo1.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(12560000,12620000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')

plot(exo1.1kb.df$start,exo1.1kb.df$df,type='l',col='darkgreen',xlim=c(12560000,12620000),ylim=c(-0.05,0.2),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')

plot(exo1.1kb.ihs.cv$start,exo1.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(12560000,12620000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')
plot(exo1.1kb.ihs.co$start,exo1.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(12560000,12620000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')

plot(exo1.1kb.beta.cv$start,exo1.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(12560000,12620000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')
plot(exo1.1kb.beta.co$start,exo1.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(12560000,12620000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=12590208,y0=vert,x1=12591465,y1=vert,col='darkblue')

# Plot scans for Exonuclease 2

# scaffold-mi3	10271502	10274220	Exonuclease	maker-scaffold-mi3-augustus-gene-34.12	crovir-transcript-11292	Exonuclease_2

par(mfrow=c(7,1))
plot(exo2.1kb.pi.cv$window_pos_1,exo2.1kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(10250000,10300000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')
plot(exo2.1kb.pi.co$window_pos_1,exo2.1kb.pi.co$avg_pi,type='l',col='navy',xlim=c(10250000,10300000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')

plot(exo2.1kb.dxy$window_pos_1,exo2.1kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(10250000,10300000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')

plot(exo2.1kb.fst$window_pos_1,exo2.1kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(10250000,10300000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')

plot(exo2.1kb.taj.cv$BIN_START,exo2.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(10250000,10300000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')
plot(exo2.1kb.taj.co$BIN_START,exo2.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(10250000,10300000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')

plot(exo2.1kb.df$start,exo2.1kb.df$df,type='l',col='darkgreen',xlim=c(10250000,10300000),ylim=c(-0.05,0.2),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')

plot(exo2.1kb.ihs.cv$start,exo2.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(10250000,10300000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')
plot(exo2.1kb.ihs.co$start,exo2.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(10250000,10300000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')

plot(exo2.1kb.beta.cv$start,exo2.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(10250000,10300000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')
plot(exo2.1kb.beta.co$start,exo2.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(10250000,10300000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=10271502,y0=vert,x1=10274220,y1=vert,col='darkblue')

# Plot scans for Exonuclease 3

# scaffold-mi7	8097114	8103411	Exonuclease	maker-scaffold-mi7-augustus-gene-27.16	crovir-transcript-1466	Exonuclease_3

par(mfrow=c(7,1))
plot(exo3.1kb.pi.cv$window_pos_1,exo3.1kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(8080000,8117000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')
plot(exo3.1kb.pi.co$window_pos_1,exo3.1kb.pi.co$avg_pi,type='l',col='navy',xlim=c(8080000,8117000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')

plot(exo3.1kb.dxy$window_pos_1,exo3.1kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(8080000,8117000),ylim=c(-0.004,0.025),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')

plot(exo3.1kb.fst$window_pos_1,exo3.1kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(8080000,8117000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')

plot(exo3.1kb.taj.cv$BIN_START,exo3.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(8080000,8117000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')
plot(exo3.1kb.taj.co$BIN_START,exo3.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(8080000,8117000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')

plot(exo3.1kb.df$start,exo3.1kb.df$df,type='l',col='darkgreen',xlim=c(8080000,8117000),ylim=c(-0.05,0.2),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')

plot(exo3.1kb.ihs.cv$start,exo3.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(8080000,8117000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')
plot(exo3.1kb.ihs.co$start,exo3.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(8080000,8117000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')

plot(exo3.1kb.beta.cv$start,exo3.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(8080000,8117000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')
plot(exo3.1kb.beta.co$start,exo3.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(8080000,8117000),ylim=c(-1,4),xlab='Chromosome Position',ylab="ß")
vert <- -0.5
segments(x0=8097114,y0=vert,x1=8103411,y1=vert,col='darkblue')

### Glutaminyl cyclases-----------------------------------------------------

# Coordinates:
# scaffold-ma1	256551622	256564040	Glutaminyl cyclase	maker-scaffold-ma1-augustus-gene-855.11	crovir-transcript-9447	vQC_1
# scaffold-mi7	5091107	5094268	Glutaminyl cyclase	maker-scaffold-mi7-augustus-gene-17.22	crovir-transcript-1322	vQC_2

# Read in data for GC1
gc1.10kb.pi.cv <- read.table('./pixy/pi.10kb.cv1.txt',header=T)
gc1.10kb.pi.cv <- gc1.10kb.pi.cv[which(gc1.10kb.pi.cv$chromosome=='scaffold-ma1'),]

gc1.10kb.pi.co <- read.table('./pixy/pi.10kb.co1.txt',header=T)
gc1.10kb.pi.co <- gc1.10kb.pi.co[which(gc1.10kb.pi.co$chromosome=='scaffold-ma1'),]

gc1.10kb.dxy <- read.table('./pixy/dxy.10kb.cv1co1.txt',header=T)
gc1.10kb.dxy <- gc1.10kb.dxy[which(gc1.10kb.dxy$chromosome=='scaffold-ma1'),]

gc1.10kb.fst <- read.table('./pixy/fst.10kb.cv1co1.txt',header=T)
gc1.10kb.fst <- gc1.10kb.fst[which(gc1.10kb.fst$chromosome=='scaffold-ma1'),]

gc1.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
gc1.1kb.taj.cv <- gc1.1kb.taj.cv[which(gc1.1kb.taj.cv$CHROM=='scaffold-ma1'),]

gc1.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
gc1.1kb.taj.co <- gc1.1kb.taj.co[which(gc1.1kb.taj.co$CHROM=='scaffold-ma1'),]

gc1.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma1.cv1.co1.txt',header=T)

gc1.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
gc1.1kb.ihs.cv <- gc1.1kb.ihs.cv[which(gc1.1kb.ihs.cv$chrom=='scaffold-ma1'),]

gc1.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
gc1.1kb.ihs.co <- gc1.1kb.ihs.co[which(gc1.1kb.ihs.co$chrom=='scaffold-ma1'),]

gc1.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
gc1.1kb.beta.cv <- gc1.1kb.beta.cv[which(gc1.1kb.beta.cv$chrom=='scaffold-ma1'),]

gc1.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
gc1.1kb.beta.co <- gc1.1kb.beta.co[which(gc1.1kb.beta.co$chrom=='scaffold-ma1'),]

# Read in data for GC2
gc2.1kb.pi.cv <- read.table('./pixy/pixy.scaffold-mi7.1kb_pi.txt',header=T)
gc2.1kb.pi.cv <- gc2.1kb.pi.cv[which(gc2.1kb.pi.cv$pop=='CV1'),]

gc2.1kb.pi.co <- read.table('./pixy/pixy.scaffold-mi7.1kb_pi.txt',header=T)
gc2.1kb.pi.co <- gc2.1kb.pi.co[which(gc2.1kb.pi.co$pop=='CO1'),]

gc2.1kb.dxy <- read.table('./pixy/pixy.scaffold-mi7.1kb_dxy.txt',header=T)
gc2.1kb.dxy <- gc2.1kb.dxy[which(gc2.1kb.dxy$pop1=='CV1' & gc2.1kb.dxy$pop2=='CO1'),]

gc2.1kb.fst <- read.table('./pixy/pixy.scaffold-mi7.1kb_fst.txt',header=T)
gc2.1kb.fst <- gc2.1kb.fst[which(gc2.1kb.fst$pop1=='CV1' & gc2.1kb.fst$pop2=='CO1'),]

gc2.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
gc2.1kb.taj.cv <- gc2.1kb.taj.cv[which(gc2.1kb.taj.cv$CHROM=='scaffold-mi7'),]

gc2.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
gc2.1kb.taj.co <- gc2.1kb.taj.co[which(gc2.1kb.taj.co$CHROM=='scaffold-mi7'),]

gc2.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-mi7.cv1.co1.txt',header=T)

gc2.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
gc2.1kb.ihs.cv <- gc2.1kb.ihs.cv[which(gc2.1kb.ihs.cv$chrom=='scaffold-mi7'),]

gc2.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
gc2.1kb.ihs.co <- gc2.1kb.ihs.co[which(gc2.1kb.ihs.co$chrom=='scaffold-mi7'),]

gc2.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
gc2.1kb.beta.cv <- gc2.1kb.beta.cv[which(gc2.1kb.beta.cv$chrom=='scaffold-mi7'),]

gc2.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
gc2.1kb.beta.co <- gc2.1kb.beta.co[which(gc2.1kb.beta.co$chrom=='scaffold-mi7'),]

# Plot scans for GC1

# scaffold-ma1	256551622	256564040	Glutaminyl cyclase	maker-scaffold-ma1-augustus-gene-855.11	crovir-transcript-9447	vQC_1

par(mfrow=c(7,1))
plot(gc1.10kb.pi.cv$window_pos_1,gc1.10kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(256400000,256800000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')
plot(gc1.10kb.pi.co$window_pos_1,gc1.10kb.pi.co$avg_pi,type='l',col='navy',xlim=c(256400000,256800000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')

plot(gc1.10kb.dxy$window_pos_1,gc1.10kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(256400000,256800000),ylim=c(-0.004,0.03),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')

plot(gc1.10kb.fst$window_pos_1,gc1.10kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(256400000,256800000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')

plot(gc1.1kb.taj.cv$BIN_START,gc1.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(256400000,256800000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')
plot(gc1.1kb.taj.co$BIN_START,gc1.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(256400000,256800000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')

plot(gc1.1kb.df$start,gc1.1kb.df$df,type='l',col='darkgreen',xlim=c(256400000,256800000),ylim=c(-0.05,0.4),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')

plot(gc1.1kb.ihs.cv$start,gc1.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(256400000,256800000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')
plot(gc1.1kb.ihs.co$start,gc1.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(256400000,256800000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')

plot(gc1.1kb.beta.cv$start,gc1.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(256400000,256800000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')
plot(gc1.1kb.beta.co$start,gc1.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(256400000,256800000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=256551622,y0=vert,x1=256564040,y1=vert,col='darkblue')

# Plot scans for GC2

# scaffold-mi7	5091107	5094268	Glutaminyl cyclase	maker-scaffold-mi7-augustus-gene-17.22	crovir-transcript-1322	vQC_2

par(mfrow=c(7,1))
plot(gc2.1kb.pi.cv$window_pos_1,gc2.1kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(5080000,5105000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')
plot(gc2.1kb.pi.co$window_pos_1,gc2.1kb.pi.co$avg_pi,type='l',col='navy',xlim=c(5080000,5105000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')

plot(gc2.1kb.dxy$window_pos_1,gc2.1kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(5080000,5105000),ylim=c(-0.004,0.03),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')

plot(gc2.1kb.fst$window_pos_1,gc2.1kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(5080000,5105000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')

plot(gc2.1kb.taj.cv$BIN_START,gc2.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(5080000,5105000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')
plot(gc2.1kb.taj.co$BIN_START,gc2.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(5080000,5105000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')

plot(gc2.1kb.df$start,gc2.1kb.df$df,type='l',col='darkgreen',xlim=c(5080000,5105000),ylim=c(-0.05,0.4),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')

plot(gc2.1kb.ihs.cv$start,gc2.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(5080000,5105000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')
plot(gc2.1kb.ihs.co$start,gc2.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(5080000,5105000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')

plot(gc2.1kb.beta.cv$start,gc2.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(5080000,5105000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')
plot(gc2.1kb.beta.co$start,gc2.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(5080000,5105000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=5091107,y0=vert,x1=5094268,y1=vert,col='darkblue')

### LAAOs-------------------------------------------------------------------

# Coordinates:
# scaffold-ma2	4654769	4658293	L-amino acid oxidase	maker-scaffold-ma2-augustus-gene-15.17	crovir-transcript-13249	LAAO_1
# scaffold-ma2	4658599	4661642	L-amino acid oxidase	maker-scaffold-ma2-augustus-gene-15.18	crovir-transcript-13246	LAAO_2
# scaffold-ma4	85461961	85468906	L-amino acid oxidase	maker-scaffold-ma4-augustus-gene-284.10	crovir-transcript-3893	LAAO_3

# Read in data for LAAO1 and LAAO2
laao12.10kb.pi.cv <- read.table('./pixy/pi.10kb.cv1.txt',header=T)
laao12.10kb.pi.cv <- laao12.10kb.pi.cv[which(laao12.10kb.pi.cv$chromosome=='scaffold-ma2'),]

laao12.10kb.pi.co <- read.table('./pixy/pi.10kb.co1.txt',header=T)
laao12.10kb.pi.co <- laao12.10kb.pi.co[which(laao12.10kb.pi.co$chromosome=='scaffold-ma2'),]

laao12.10kb.dxy <- read.table('./pixy/dxy.10kb.cv1co1.txt',header=T)
laao12.10kb.dxy <- laao12.10kb.dxy[which(laao12.10kb.dxy$chromosome=='scaffold-ma2'),]

laao12.10kb.fst <- read.table('./pixy/fst.10kb.cv1co1.txt',header=T)
laao12.10kb.fst <- laao12.10kb.fst[which(laao12.10kb.fst$chromosome=='scaffold-ma2'),]

laao12.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
laao12.1kb.taj.cv <- laao12.1kb.taj.cv[which(laao12.1kb.taj.cv$CHROM=='scaffold-ma2'),]

laao12.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
laao12.1kb.taj.co <- laao12.1kb.taj.co[which(laao12.1kb.taj.co$CHROM=='scaffold-ma2'),]

laao12.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma2.cv1.co1.txt',header=T)

laao12.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
laao12.1kb.ihs.cv <- laao12.1kb.ihs.cv[which(laao12.1kb.ihs.cv$chrom=='scaffold-ma2'),]

laao12.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
laao12.1kb.ihs.co <- laao12.1kb.ihs.co[which(laao12.1kb.ihs.co$chrom=='scaffold-ma2'),]

laao12.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
laao12.1kb.beta.cv <- laao12.1kb.beta.cv[which(laao12.1kb.beta.cv$chrom=='scaffold-ma2'),]

laao12.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
laao12.1kb.beta.co <- laao12.1kb.beta.co[which(laao12.1kb.beta.co$chrom=='scaffold-ma2'),]

# Read in data for LAAO3
laao3.1kb.pi.cv <- read.table('./pixy/pixy.scaffold-ma4.1kb_pi.txt',header=T)
laao3.1kb.pi.cv <- laao3.1kb.pi.cv[which(laao3.1kb.pi.cv$pop=='CV1'),]

laao3.1kb.pi.co <- read.table('./pixy/pixy.scaffold-ma4.1kb_pi.txt',header=T)
laao3.1kb.pi.co <- laao3.1kb.pi.co[which(laao3.1kb.pi.co$pop=='CO1'),]

laao3.1kb.dxy <- read.table('./pixy/pixy.scaffold-ma4.1kb_dxy.txt',header=T)
laao3.1kb.dxy <- laao3.1kb.dxy[which(laao3.1kb.dxy$pop1=='CV1' & laao3.1kb.dxy$pop2=='CO1'),]

laao3.1kb.fst <- read.table('./pixy/pixy.scaffold-ma4.1kb_fst.txt',header=T)
laao3.1kb.fst <- laao3.1kb.fst[which(laao3.1kb.fst$pop1=='CV1' & laao3.1kb.fst$pop2=='CO1'),]

laao3.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
laao3.1kb.taj.cv <- laao3.1kb.taj.cv[which(laao3.1kb.taj.cv$CHROM=='scaffold-ma4'),]

laao3.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
laao3.1kb.taj.co <- laao3.1kb.taj.co[which(laao3.1kb.taj.co$CHROM=='scaffold-ma4'),]

laao3.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma4.cv1.co1.txt',header=T)

laao3.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
laao3.1kb.ihs.cv <- laao3.1kb.ihs.cv[which(laao3.1kb.ihs.cv$chrom=='scaffold-ma4'),]

laao3.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
laao3.1kb.ihs.co <- laao3.1kb.ihs.co[which(laao3.1kb.ihs.co$chrom=='scaffold-ma4'),]

laao3.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
laao3.1kb.beta.cv <- laao3.1kb.beta.cv[which(laao3.1kb.beta.cv$chrom=='scaffold-ma4'),]

laao3.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
laao3.1kb.beta.co <- laao3.1kb.beta.co[which(laao3.1kb.beta.co$chrom=='scaffold-ma4'),]

# Plot scans for LAAO 1 and 2

# scaffold-ma2	4654769	4658293	L-amino acid oxidase	maker-scaffold-ma2-augustus-gene-15.17	crovir-transcript-13249	LAAO_1
# scaffold-ma2	4658599	4661642	L-amino acid oxidase	maker-scaffold-ma2-augustus-gene-15.18	crovir-transcript-13246	LAAO_2

par(mfrow=c(7,1))
plot(laao12.10kb.pi.cv$window_pos_1,laao12.10kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(4620000,4690000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')
plot(laao12.10kb.pi.co$window_pos_1,laao12.10kb.pi.co$avg_pi,type='l',col='navy',xlim=c(4620000,4690000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')

plot(laao12.10kb.dxy$window_pos_1,laao12.10kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(4620000,4690000),ylim=c(-0.004,0.03),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')

plot(laao12.10kb.fst$window_pos_1,laao12.10kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(4620000,4690000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')

plot(laao12.1kb.taj.cv$BIN_START,laao12.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(4620000,4690000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')
plot(laao12.1kb.taj.co$BIN_START,laao12.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(4620000,4690000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')

plot(laao12.1kb.df$start,laao12.1kb.df$df,type='l',col='darkgreen',xlim=c(4620000,4690000),ylim=c(-0.05,0.4),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')

plot(laao12.1kb.ihs.cv$start,laao12.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(4620000,4690000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')
plot(laao12.1kb.ihs.co$start,laao12.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(4620000,4690000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')

plot(laao12.1kb.beta.cv$start,laao12.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(4620000,4690000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')
plot(laao12.1kb.beta.co$start,laao12.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(4620000,4690000),ylim=c(-3,20),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=4654769,y0=vert,x1=4658293,y1=vert,col='darkblue')
segments(x0=4658599,y0=vert,x1=4661642,y1=vert,col='darkblue')

# Plot scans for LAAO3

# scaffold-ma4	85461961	85468906	L-amino acid oxidase	maker-scaffold-ma4-augustus-gene-284.10	crovir-transcript-3893	LAAO_3

par(mfrow=c(7,1))
plot(laao3.1kb.pi.cv$window_pos_1,laao3.1kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(85420000,85500000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')
plot(laao3.1kb.pi.co$window_pos_1,laao3.1kb.pi.co$avg_pi,type='l',col='navy',xlim=c(85420000,85500000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')

plot(laao3.1kb.dxy$window_pos_1,laao3.1kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(85420000,85500000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')

plot(laao3.1kb.fst$window_pos_1,laao3.1kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(85420000,85500000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')

plot(laao3.1kb.taj.cv$BIN_START,laao3.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(85420000,85500000),ylim=c(-2,2),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')
plot(laao3.1kb.taj.co$BIN_START,laao3.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(85420000,85500000),ylim=c(-2,2),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -1.5
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')

plot(laao3.1kb.df$start,laao3.1kb.df$df,type='l',col='darkgreen',xlim=c(85420000,85500000),ylim=c(-0.05,0.2),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')

plot(laao3.1kb.ihs.cv$start,laao3.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(85420000,85500000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')
plot(laao3.1kb.ihs.co$start,laao3.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(85420000,85500000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')

plot(laao3.1kb.beta.cv$start,laao3.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(85420000,85500000),ylim=c(-4,28),xlab='Chromosome Position',ylab="ß")
vert <- -3
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')
plot(laao3.1kb.beta.co$start,laao3.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(85420000,85500000),ylim=c(-4,28),xlab='Chromosome Position',ylab="ß")
vert <- -3
segments(x0=85461961,y0=vert,x1=85468906,y1=vert,col='darkblue')


### VEGFs-------------------------------------------------------------------

# Coordinates:
# scaffold-ma1	260248287	260272500	Vascular endothelial growth factor	maker-scaffold-ma1-augustus-gene-867.2	crovir-transcript-9439	VEGF_1
# scaffold-ma7	40288572	40327884	Vascular endothelial growth factor	maker-scaffold-ma7-augustus-gene-134.2	crovir-transcript-4609	VEGF_2

# Read in data for VEGF1
vgf1.10kb.pi.cv <- read.table('./pixy/pi.10kb.cv1.txt',header=T)
vgf1.10kb.pi.cv <- vgf1.10kb.pi.cv[which(vgf1.10kb.pi.cv$chromosome=='scaffold-ma1'),]

vgf1.10kb.pi.co <- read.table('./pixy/pi.10kb.co1.txt',header=T)
vgf1.10kb.pi.co <- vgf1.10kb.pi.co[which(vgf1.10kb.pi.co$chromosome=='scaffold-ma1'),]

vgf1.10kb.dxy <- read.table('./pixy/dxy.10kb.cv1co1.txt',header=T)
vgf1.10kb.dxy <- vgf1.10kb.dxy[which(vgf1.10kb.dxy$chromosome=='scaffold-ma1'),]

vgf1.10kb.fst <- read.table('./pixy/fst.10kb.cv1co1.txt',header=T)
vgf1.10kb.fst <- vgf1.10kb.fst[which(vgf1.10kb.fst$chromosome=='scaffold-ma1'),]

vgf1.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
vgf1.1kb.taj.cv <- vgf1.1kb.taj.cv[which(vgf1.1kb.taj.cv$CHROM=='scaffold-ma1'),]

vgf1.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
vgf1.1kb.taj.co <- vgf1.1kb.taj.co[which(vgf1.1kb.taj.co$CHROM=='scaffold-ma1'),]

vgf1.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma1.cv1.co1.txt',header=T)

vgf1.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
vgf1.1kb.ihs.cv <- vgf1.1kb.ihs.cv[which(vgf1.1kb.ihs.cv$chrom=='scaffold-ma1'),]

vgf1.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
vgf1.1kb.ihs.co <- vgf1.1kb.ihs.co[which(vgf1.1kb.ihs.co$chrom=='scaffold-ma1'),]

vgf1.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
vgf1.1kb.beta.cv <- vgf1.1kb.beta.cv[which(vgf1.1kb.beta.cv$chrom=='scaffold-ma1'),]

vgf1.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
vgf1.1kb.beta.co <- vgf1.1kb.beta.co[which(vgf1.1kb.beta.co$chrom=='scaffold-ma1'),]

# Read in data for VEGF2
vgf2.1kb.pi.cv <- read.table('./pixy/pixy.scaffold-ma7.1kb_pi.txt',header=T)
vgf2.1kb.pi.cv <- vgf2.1kb.pi.cv[which(vgf2.1kb.pi.cv$pop=='CV1'),]

vgf2.1kb.pi.co <- read.table('./pixy/pixy.scaffold-ma7.1kb_pi.txt',header=T)
vgf2.1kb.pi.co <- vgf2.1kb.pi.co[which(vgf2.1kb.pi.co$pop=='CO1'),]

vgf2.1kb.dxy <- read.table('./pixy/pixy.scaffold-ma7.1kb_dxy.txt',header=T)
vgf2.1kb.dxy <- vgf2.1kb.dxy[which(vgf2.1kb.dxy$pop1=='CV1' & vgf2.1kb.dxy$pop2=='CO1'),]

vgf2.1kb.fst <- read.table('./pixy/pixy.scaffold-ma7.1kb_fst.txt',header=T)
vgf2.1kb.fst <- vgf2.1kb.fst[which(vgf2.1kb.fst$pop1=='CV1' & vgf2.1kb.fst$pop2=='CO1'),]

vgf2.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
vgf2.1kb.taj.cv <- vgf2.1kb.taj.cv[which(vgf2.1kb.taj.cv$CHROM=='scaffold-ma7'),]

vgf2.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
vgf2.1kb.taj.co <- vgf2.1kb.taj.co[which(vgf2.1kb.taj.co$CHROM=='scaffold-ma7'),]

vgf2.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma7.cv1.co1.txt',header=T)

vgf2.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
vgf2.1kb.ihs.cv <- vgf2.1kb.ihs.cv[which(vgf2.1kb.ihs.cv$chrom=='scaffold-ma7'),]

vgf2.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
vgf2.1kb.ihs.co <- vgf2.1kb.ihs.co[which(vgf2.1kb.ihs.co$chrom=='scaffold-ma7'),]

vgf2.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
vgf2.1kb.beta.cv <- vgf2.1kb.beta.cv[which(vgf2.1kb.beta.cv$chrom=='scaffold-ma7'),]

vgf2.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
vgf2.1kb.beta.co <- vgf2.1kb.beta.co[which(vgf2.1kb.beta.co$chrom=='scaffold-ma7'),]

# Plot scans for VEGF1

# scaffold-ma1	260248287	260272500	Vascular endothelial growth factor	maker-scaffold-ma1-augustus-gene-867.2	crovir-transcript-9439	VEGF_1

par(mfrow=c(7,1))
plot(vgf1.10kb.pi.cv$window_pos_1,vgf1.10kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(260220000,260290000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')
plot(vgf1.10kb.pi.co$window_pos_1,vgf1.10kb.pi.co$avg_pi,type='l',col='navy',xlim=c(260220000,260290000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')

plot(vgf1.10kb.dxy$window_pos_1,vgf1.10kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(260220000,260290000),ylim=c(-0.004,0.03),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')

plot(vgf1.10kb.fst$window_pos_1,vgf1.10kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(260220000,260290000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')

plot(vgf1.1kb.taj.cv$BIN_START,vgf1.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(260220000,260290000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')
plot(vgf1.1kb.taj.co$BIN_START,vgf1.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(260220000,260290000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')

plot(vgf1.1kb.df$start,vgf1.1kb.df$df,type='l',col='darkgreen',xlim=c(260220000,260290000),ylim=c(-0.05,0.4),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')

plot(vgf1.1kb.ihs.cv$start,vgf1.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(260220000,260290000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')
plot(vgf1.1kb.ihs.co$start,vgf1.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(260220000,260290000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')

plot(vgf1.1kb.beta.cv$start,vgf1.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(260220000,260290000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')
plot(vgf1.1kb.beta.co$start,vgf1.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(260220000,260290000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=260248287,y0=vert,x1=260272500,y1=vert,col='darkblue')

# Plot scans for VEGF2

# scaffold-ma7	40288572	40327884	Vascular endothelial growth factor	maker-scaffold-ma7-augustus-gene-134.2	crovir-transcript-4609	VEGF_2

par(mfrow=c(7,1))
plot(vgf2.1kb.pi.cv$window_pos_1,vgf2.1kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(40250000,40400000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')
plot(vgf2.1kb.pi.co$window_pos_1,vgf2.1kb.pi.co$avg_pi,type='l',col='navy',xlim=c(40250000,40400000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')

plot(vgf2.1kb.dxy$window_pos_1,vgf2.1kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(40250000,40400000),ylim=c(-0.004,0.03),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')

plot(vgf2.1kb.fst$window_pos_1,vgf2.1kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(40250000,40400000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')

plot(vgf2.1kb.taj.cv$BIN_START,vgf2.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(40250000,40400000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')
plot(vgf2.1kb.taj.co$BIN_START,vgf2.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(40250000,40400000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')

plot(vgf2.1kb.df$start,vgf2.1kb.df$df,type='l',col='darkgreen',xlim=c(40250000,40400000),ylim=c(-0.05,0.4),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')

plot(vgf2.1kb.ihs.cv$start,vgf2.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(40250000,40400000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')
plot(vgf2.1kb.ihs.co$start,vgf2.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(40250000,40400000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')

plot(vgf2.1kb.beta.cv$start,vgf2.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(40250000,40400000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')
plot(vgf2.1kb.beta.co$start,vgf2.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(40250000,40400000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=40288572,y0=vert,x1=40327884,y1=vert,col='darkblue')

### Vespryn-----------------------------------------------------------------

# Coordinates:
# scaffold-ma2	4377779	4385668	Vespryn/Ohanin	maker-scaffold-ma2-augustus-gene-14.15	crovir-transcript-13215	Vespryn

# Read in data for vespryn
vesp.10kb.pi.cv <- read.table('./pixy/pi.10kb.cv1.txt',header=T)
vesp.10kb.pi.cv <- vesp.10kb.pi.cv[which(vesp.10kb.pi.cv$chromosome=='scaffold-ma2'),]

vesp.10kb.pi.co <- read.table('./pixy/pi.10kb.co1.txt',header=T)
vesp.10kb.pi.co <- vesp.10kb.pi.co[which(vesp.10kb.pi.co$chromosome=='scaffold-ma2'),]

vesp.10kb.dxy <- read.table('./pixy/dxy.10kb.cv1co1.txt',header=T)
vesp.10kb.dxy <- vesp.10kb.dxy[which(vesp.10kb.dxy$chromosome=='scaffold-ma2'),]

vesp.10kb.fst <- read.table('./pixy/fst.10kb.cv1co1.txt',header=T)
vesp.10kb.fst <- vesp.10kb.fst[which(vesp.10kb.fst$chromosome=='scaffold-ma2'),]

vesp.1kb.taj.cv <- read.table('./tajima_d/cv.colorado.all.1kb.Tajima.D',header=T)
vesp.1kb.taj.cv <- vesp.1kb.taj.cv[which(vesp.1kb.taj.cv$CHROM=='scaffold-ma2'),]

vesp.1kb.taj.co <- read.table('./tajima_d/co.california.all.1kb.Tajima.D',header=T)
vesp.1kb.taj.co <- vesp.1kb.taj.co[which(vesp.1kb.taj.co$CHROM=='scaffold-ma2'),]

vesp.1kb.df <- read.table('./df/window.1kb.df_prop.scaffold-ma2.cv1.co1.txt',header=T)

vesp.1kb.ihs.cv <- read.table('./rehh/cv.all_ihs.1kb.txt',header=T)
vesp.1kb.ihs.cv <- vesp.1kb.ihs.cv[which(vesp.1kb.ihs.cv$chrom=='scaffold-ma2'),]

vesp.1kb.ihs.co <- read.table('./rehh/co.all_ihs.1kb.txt',header=T)
vesp.1kb.ihs.co <- vesp.1kb.ihs.co[which(vesp.1kb.ihs.co$chrom=='scaffold-ma2'),]

vesp.1kb.beta.cv <- read.table('./beta/results/cv1.phased.all.betascores.1kb.txt',header=T)
vesp.1kb.beta.cv <- vesp.1kb.beta.cv[which(vesp.1kb.beta.cv$chrom=='scaffold-ma2'),]

vesp.1kb.beta.co <- read.table('./beta/results/co1.phased.all.betascores.1kb.txt',header=T)
vesp.1kb.beta.co <- vesp.1kb.beta.co[which(vesp.1kb.beta.co$chrom=='scaffold-ma2'),]

# Plot scans for vespryn

# scaffold-ma2	4377779	4385668	Vespryn/Ohanin	maker-scaffold-ma2-augustus-gene-14.15	crovir-transcript-13215	Vespryn

par(mfrow=c(7,1))
plot(vesp.10kb.pi.cv$window_pos_1,vesp.10kb.pi.cv$avg_pi,type='l',col='darkorange',xlim=c(4300000,4450000),ylim=c(-0.004,0.02),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')
plot(vesp.10kb.pi.co$window_pos_1,vesp.10kb.pi.co$avg_pi,type='l',col='navy',xlim=c(4300000,4450000),ylim=c(-0.004,0.04),xlab='Chromosome Position',ylab='pi')
vert <- -0.002
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')

plot(vesp.10kb.dxy$window_pos_1,vesp.10kb.dxy$avg_dxy,type='l',col='darkgreen',xlim=c(4300000,4450000),ylim=c(-0.004,0.03),xlab='Chromosome Position',ylab='dxy')
vert <- -0.002
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')

plot(vesp.10kb.fst$window_pos_1,vesp.10kb.fst$avg_wc_fst,type='l',col='darkgreen',xlim=c(4300000,4450000),ylim=c(-0.1,1),xlab='Chromosome Position',ylab='fst')
vert <- -0.002
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')

plot(vesp.1kb.taj.cv$BIN_START,vesp.1kb.taj.cv$TajimaD,type='l',col='darkorange',xlim=c(4300000,4450000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')
plot(vesp.1kb.taj.co$BIN_START,vesp.1kb.taj.co$TajimaD,type='l',col='navy',xlim=c(4300000,4450000),ylim=c(-3,3),xlab='Chromosome Position',ylab="Tajima's D")
vert <- -2.5
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')

plot(vesp.1kb.df$start,vesp.1kb.df$df,type='l',col='darkgreen',xlim=c(4300000,4450000),ylim=c(-0.05,0.4),xlab='Chromosome Position',ylab='df')
vert <- -0.025
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')

plot(vesp.1kb.ihs.cv$start,vesp.1kb.ihs.cv$iHS,type='l',col='darkorange',xlim=c(4300000,4450000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')
plot(vesp.1kb.ihs.co$start,vesp.1kb.ihs.co$iHS,type='l',col='navy',xlim=c(4300000,4450000),ylim=c(-2,2),xlab='Chromosome Position',ylab="iHS")
vert <- -1.5
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')

plot(vesp.1kb.beta.cv$start,vesp.1kb.beta.cv$Beta1.,type='l',col='darkorange',xlim=c(4300000,4450000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')
plot(vesp.1kb.beta.co$start,vesp.1kb.beta.co$Beta1.,type='l',col='navy',xlim=c(4300000,4450000),ylim=c(-3,12),xlab='Chromosome Position',ylab="ß")
vert <- -2
segments(x0=4377779,y0=vert,x1=4385668,y1=vert,col='darkblue')

