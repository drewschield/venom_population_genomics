############################################################################
# Plotting ADMIXTURE results from C. viridis & C. oreganus populations
############################################################################

### Goal: plot Q tables from ADMIXTURE analyses of C. viridis and C. oreganus
### populations from Colorado/Montana and California/Idaho, respectively.

### Set working directory---------------------------------------------------

setwd('./population_structure/')

### Install and load libraries/dependencies---------------------------------

devtools::install_github('royfrancis/pophelper')
install.packages('viridis')
library(pophelper)
library(viridis)

### Read in data------------------------------------------------------------

# CV error
CVerror <- read.table('./CV_error.txt',header=T)

# Q tables
k2 <- readQ('./cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04.2.Q',filetype='auto')
k3 <- readQ('./cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04.3.Q',filetype='auto')
k4 <- readQ('./cvco.ingroup.mask.HardFilter.chrom.snp.maf05.thin1kb.miss04.4.Q',filetype='auto')

# Sample metadata
inds_grps <- read.table('sampleID_pop.txt',header=F,stringsAsFactors = F)
grps <- as.data.frame(inds_grps[,2])
grps[,1] <- as.character(grps[,1])
colnames(grps) <- c('subs')

# Add sample IDs as rownames to qlist
rownames(k2[[1]]) <- inds_grps$V1
rownames(k3[[1]]) <- inds_grps$V1
rownames(k4[[1]]) <- inds_grps$V1
rownames(k5[[1]]) <- inds_grps$V1

### Plot CV error------------------------------------------------------------

par(mfrow=c(1,1))
plot(CVerror$K,CVerror$CV_error,pch=20,xlab='K',ylab='CV error')
low <- min(CVerror$CV_error)
abline(h=low,col='navy',lty=2,lwd=0.5)

### Plot ADMIXTURE plots-----------------------------------------------------

# set population order
pop_order <- c('CV1','CV2','CO1','CO2')

# Groups individuals by subspecies abbreviation

plotQ(k2,imgtype = 'pdf',basesize=11,
      height=5,width=30,
      clustercol = viridis(2),
      showlegend = F,legendtextsize = 3, legendmargin = c(2,4,4,0),
      showyaxis = T,
      sortind = 'all',
      showindlab = T, useindlab = T, indlabsize = 5,sharedindlab = F,
      linesize=0.8,pointsize = 4,
      grplab = grps,grplabsize = 3,grplabangle = 45, grplabheight = 0.1,ordergrp = T,selgrp = 'subs',
      returndata = T,
      outputfilename = '/k2_AdmixtureROUGH_03.08.21',
      exportpath=getwd()
)


plotQ(k3,imgtype = 'pdf',basesize=11,
      height=5,width=30,
      clustercol = viridis(3),
      showlegend = F,legendtextsize = 3, legendmargin = c(2,4,4,0),
      showyaxis = T,
      sortind = 'all',
      showindlab = T, useindlab = T, indlabsize = 5,sharedindlab = F,
      linesize=0.8,pointsize = 4,
      grplab = grps,subsetgrp = pop_order,grplabsize = 3,grplabangle = 45, grplabheight = 0.1,ordergrp = T,selgrp = 'subs',
      returndata = T,
      outputfilename = '/k3_AdmixtureROUGH_03.08.21',
      exportpath=getwd()
)

plotQ(k4,imgtype = 'pdf',basesize=11,
      height=5,width=30,
      clustercol = viridis(4),
      showlegend = F,legendtextsize = 3, legendmargin = c(2,4,4,0),
      showyaxis = T,
      sortind = 'all',
      showindlab = T, useindlab = T, indlabsize = 5,sharedindlab = F,
      linesize=0.8,pointsize = 4,
      grplab = grps,grplabsize = 3,grplabangle = 45, grplabheight = 0.1,ordergrp = T,selgrp = 'subs',
      returndata = T,
      outputfilename = '/k4_AdmixtureROUGH_03.08.21',
      exportpath=getwd()
)
