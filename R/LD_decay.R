############################################################################
# Linkage disequilibrium (LD) decay in venom and flanking regions
############################################################################

### Goal: examine rate of LD decay with increasing physical distance between
### SNPs in venom and flanking regions.

### Steps in this script were adapted from rmf's post here: https://www.biostars.org/p/300381/

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory and load dependencies-----------------------------

setwd('./recombination/ld_decay/')

library(dplyr)
library(stringr)
library(ggplot2)

### Read in data-------------------------------------------------------------

cv.svmp.dfr <- read.delim("viridis.phased.svmp.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
cv.svmp.ffr <- read.delim("viridis.phased.svmp-flank.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)

cv.svsp.dfr <- read.delim("viridis.phased.svsp.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
cv.svsp.ffr <- read.delim("viridis.phased.svsp-flank.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)

cv.pla2.dfr <- read.delim("viridis.phased.pla2.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
cv.pla2.ffr <- read.delim("viridis.phased.pla2-flank.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)

co.svmp.dfr <- read.delim("oreganus.phased.svmp.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
co.svmp.ffr <- read.delim("oreganus.phased.svmp-flank.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)

co.svsp.dfr <- read.delim("oreganus.phased.svsp.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
co.svsp.ffr <- read.delim("oreganus.phased.svsp-flank.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)

co.pla2.dfr <- read.delim("oreganus.phased.pla2.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
co.pla2.ffr <- read.delim("oreganus.phased.pla2-flank.hap.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)

### Calculate mean and interquartile values with physical distance-----------

## C. viridis
# SVMP and flanking region
colnames(cv.svmp.dfr) <- c("dist","rsq")
cv.svmp.dfr$distc <- cut(cv.svmp.dfr$dist,breaks=seq(from=min(cv.svmp.dfr$dist)-1,to=max(cv.svmp.dfr$dist)+1,by=100))
cv.svmp.dfr1 <- cv.svmp.dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
cv.svmp.dfrh <- cv.svmp.dfr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
cv.svmp.dfrl <- cv.svmp.dfr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
cv.svmp.dfr1 <- cv.svmp.dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svmp.dfrh <- cv.svmp.dfrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svmp.dfrl <- cv.svmp.dfrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

colnames(cv.svmp.ffr) <- c("dist","rsq")
cv.svmp.ffr$distc <- cut(cv.svmp.ffr$dist,breaks=seq(from=min(cv.svmp.ffr$dist)-1,to=max(cv.svmp.ffr$dist)+1,by=100))
cv.svmp.ffr1 <- cv.svmp.ffr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
cv.svmp.ffrh <- cv.svmp.ffr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
cv.svmp.ffrl <- cv.svmp.ffr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
cv.svmp.ffr1 <- cv.svmp.ffr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svmp.ffrh <- cv.svmp.ffrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svmp.ffrl <- cv.svmp.ffrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

# SVSP and flanking region
colnames(cv.svsp.dfr) <- c("dist","rsq")
cv.svsp.dfr$distc <- cut(cv.svsp.dfr$dist,breaks=seq(from=min(cv.svsp.dfr$dist)-1,to=max(cv.svsp.dfr$dist)+1,by=100))
cv.svsp.dfr1 <- cv.svsp.dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
cv.svsp.dfrh <- cv.svsp.dfr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
cv.svsp.dfrl <- cv.svsp.dfr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
cv.svsp.dfr1 <- cv.svsp.dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svsp.dfrh <- cv.svsp.dfrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svsp.dfrl <- cv.svsp.dfrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

colnames(cv.svsp.ffr) <- c("dist","rsq")
cv.svsp.ffr$distc <- cut(cv.svsp.ffr$dist,breaks=seq(from=min(cv.svsp.ffr$dist)-1,to=max(cv.svsp.ffr$dist)+1,by=100))
cv.svsp.ffr1 <- cv.svsp.ffr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
cv.svsp.ffrh <- cv.svsp.ffr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
cv.svsp.ffrl <- cv.svsp.ffr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
cv.svsp.ffr1 <- cv.svsp.ffr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svsp.ffrh <- cv.svsp.ffrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.svsp.ffrl <- cv.svsp.ffrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

# PLA2 and flanking region
colnames(cv.pla2.dfr) <- c("dist","rsq")
cv.pla2.dfr$distc <- cut(cv.pla2.dfr$dist,breaks=seq(from=min(cv.pla2.dfr$dist)-1,to=max(cv.pla2.dfr$dist)+1,by=100))
cv.pla2.dfr1 <- cv.pla2.dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
cv.pla2.dfrh <- cv.pla2.dfr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
cv.pla2.dfrl <- cv.pla2.dfr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
cv.pla2.dfr1 <- cv.pla2.dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.pla2.dfrh <- cv.pla2.dfrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.pla2.dfrl <- cv.pla2.dfrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

colnames(cv.pla2.ffr) <- c("dist","rsq")
cv.pla2.ffr$distc <- cut(cv.pla2.ffr$dist,breaks=seq(from=min(cv.pla2.ffr$dist)-1,to=max(cv.pla2.ffr$dist)+1,by=100))
cv.pla2.ffr1 <- cv.pla2.ffr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
cv.pla2.ffrh <- cv.pla2.ffr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
cv.pla2.ffrl <- cv.pla2.ffr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
cv.pla2.ffr1 <- cv.pla2.ffr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.pla2.ffrh <- cv.pla2.ffrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
cv.pla2.ffrl <- cv.pla2.ffrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

## C. oreganus
# SVMP and flanking region
colnames(co.svmp.dfr) <- c("dist","rsq")
co.svmp.dfr$distc <- cut(co.svmp.dfr$dist,breaks=seq(from=min(co.svmp.dfr$dist)-1,to=max(co.svmp.dfr$dist)+1,by=100))
co.svmp.dfr1 <- co.svmp.dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
co.svmp.dfrh <- co.svmp.dfr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
co.svmp.dfrl <- co.svmp.dfr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
co.svmp.dfr1 <- co.svmp.dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svmp.dfrh <- co.svmp.dfrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svmp.dfrl <- co.svmp.dfrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))

colnames(co.svmp.ffr) <- c("dist","rsq")
co.svmp.ffr$distc <- cut(co.svmp.ffr$dist,breaks=seq(from=min(co.svmp.ffr$dist)-1,to=max(co.svmp.ffr$dist)+1,by=100))
co.svmp.ffr1 <- co.svmp.ffr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
co.svmp.ffrh <- co.svmp.ffr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
co.svmp.ffrl <- co.svmp.ffr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
co.svmp.ffr1 <- co.svmp.ffr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svmp.ffrh <- co.svmp.ffrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svmp.ffrl <- co.svmp.ffrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))

# SVSP and flanking region
colnames(co.svsp.dfr) <- c("dist","rsq")
co.svsp.dfr$distc <- cut(co.svsp.dfr$dist,breaks=seq(from=min(co.svsp.dfr$dist)-1,to=max(co.svsp.dfr$dist)+1,by=100))
co.svsp.dfr1 <- co.svsp.dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
co.svsp.dfrh <- co.svsp.dfr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
co.svsp.dfrl <- co.svsp.dfr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
co.svsp.dfr1 <- co.svsp.dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svsp.dfrh <- co.svsp.dfrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svsp.dfrl <- co.svsp.dfrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))

colnames(co.svsp.ffr) <- c("dist","rsq")
co.svsp.ffr$distc <- cut(co.svsp.ffr$dist,breaks=seq(from=min(co.svsp.ffr$dist)-1,to=max(co.svsp.ffr$dist)+1,by=100))
co.svsp.ffr1 <- co.svsp.ffr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
co.svsp.ffrh <- co.svsp.ffr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
co.svsp.ffrl <- co.svsp.ffr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
co.svsp.ffr1 <- co.svsp.ffr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svsp.ffrh <- co.svsp.ffrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.svsp.ffrl <- co.svsp.ffrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))

# PLA2 and flanking region
colnames(co.pla2.dfr) <- c("dist","rsq")
co.pla2.dfr$distc <- cut(co.pla2.dfr$dist,breaks=seq(from=min(co.pla2.dfr$dist)-1,to=max(co.pla2.dfr$dist)+1,by=100))
co.pla2.dfr1 <- co.pla2.dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
co.pla2.dfrh <- co.pla2.dfr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
co.pla2.dfrl <- co.pla2.dfr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
co.pla2.dfr1 <- co.pla2.dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.pla2.dfrh <- co.pla2.dfrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.pla2.dfrl <- co.pla2.dfrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))

colnames(co.pla2.ffr) <- c("dist","rsq")
co.pla2.ffr$distc <- cut(co.pla2.ffr$dist,breaks=seq(from=min(co.pla2.ffr$dist)-1,to=max(co.pla2.ffr$dist)+1,by=100))
co.pla2.ffr1 <- co.pla2.ffr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
co.pla2.ffrh <- co.pla2.ffr %>% group_by(distc) %>% summarise(up=quantile(rsq,0.75))
co.pla2.ffrl <- co.pla2.ffr %>% group_by(distc) %>% summarise(do=quantile(rsq,0.25))
co.pla2.ffr1 <- co.pla2.ffr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.pla2.ffrh <- co.pla2.ffrh %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))
co.pla2.ffrl <- co.pla2.ffrl %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                                        mid=start+((end-start)/2))

### Plot LD decay------------------------------------------------------------

par(mfrow=c(3,4))
plot(cv.svmp.dfr1$start,cv.svmp.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(cv.svmp.dfrh$start,cv.svmp.dfrh$up,col=alpha('lightblue3',0.4),pch=20)
lines(cv.svmp.dfrl$start,cv.svmp.dfrl$do,col=alpha('lightblue3',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(cv.svmp.ffr1$start,cv.svmp.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(cv.svmp.ffrh$start,cv.svmp.ffrh$up,col=alpha('grey',0.4),pch=20)
lines(cv.svmp.ffrl$start,cv.svmp.ffrl$do,col=alpha('grey',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(co.svmp.dfr1$start,co.svmp.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(co.svmp.dfrh$start,co.svmp.dfrh$up,col=alpha('lightblue3',0.4),pch=20)
lines(co.svmp.dfrl$start,co.svmp.dfrl$do,col=alpha('lightblue3',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(co.svmp.ffr1$start,co.svmp.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(co.svmp.ffrh$start,co.svmp.ffrh$up,col=alpha('grey',0.4),pch=20)
lines(co.svmp.ffrl$start,co.svmp.ffrl$do,col=alpha('grey',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(cv.svsp.dfr1$start,cv.svsp.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(cv.svsp.dfrh$start,cv.svsp.dfrh$up,col=alpha('aquamarine3',0.4),pch=20)
lines(cv.svsp.dfrl$start,cv.svsp.dfrl$do,col=alpha('aquamarine3',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(cv.svsp.ffr1$start,cv.svsp.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(cv.svsp.ffrh$start,cv.svsp.ffrh$up,col=alpha('grey',0.4),pch=20)
lines(cv.svsp.ffrl$start,cv.svsp.ffrl$do,col=alpha('grey',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(co.svsp.dfr1$start,co.svsp.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(co.svsp.dfrh$start,co.svsp.dfrh$up,col=alpha('aquamarine3',0.4),pch=20)
lines(co.svsp.dfrl$start,co.svsp.dfrl$do,col=alpha('aquamarine3',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(co.svsp.ffr1$start,co.svsp.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8),ylab='r2',xlab='Physical Distance (bp)')
lines(co.svsp.ffrh$start,co.svsp.ffrh$up,col=alpha('grey',0.4),pch=20)
lines(co.svsp.ffrl$start,co.svsp.ffrl$do,col=alpha('grey',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(cv.pla2.dfr1$start,cv.pla2.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,1),ylab='r2',xlab='Physical Distance (bp)')
lines(cv.pla2.dfrh$start,cv.pla2.dfrh$up,col=alpha('maroon',0.4),pch=20)
lines(cv.pla2.dfrl$start,cv.pla2.dfrl$do,col=alpha('maroon',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(cv.pla2.ffr1$start,cv.pla2.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,1),ylab='r2',xlab='Physical Distance (bp)')
lines(cv.pla2.ffrh$start,cv.pla2.ffrh$up,col=alpha('grey',0.4),pch=20)
lines(cv.pla2.ffrl$start,cv.pla2.ffrl$do,col=alpha('grey',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(co.pla2.dfr1$start,co.pla2.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,1),ylab='r2',xlab='Physical Distance (bp)')
lines(co.pla2.dfrh$start,co.pla2.dfrh$up,col=alpha('maroon',0.4),pch=20)
lines(co.pla2.dfrl$start,co.pla2.dfrl$do,col=alpha('maroon',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))

plot(co.pla2.ffr1$start,co.pla2.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,1),ylab='r2',xlab='Physical Distance (bp)')
lines(co.pla2.ffrh$start,co.pla2.ffrh$up,col=alpha('grey',0.4),pch=20)
lines(co.pla2.ffrl$start,co.pla2.ffrl$do,col=alpha('grey',0.4),pch=20)
abline(h=0.2,lty=3,col=alpha('black',0.5))


plot(cv.svsp.dfr1$start,cv.svsp.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(cv.svsp.dfrh$start,cv.svsp.dfrh$up,col=alpha('red',0.15),pch=20)
lines(cv.svsp.dfrl$start,cv.svsp.dfrl$do,col=alpha('red',0.15),pch=20)

plot(cv.svsp.ffr1$start,cv.svsp.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(cv.svsp.ffrh$start,cv.svsp.ffrh$up,col=alpha('blue',0.15),pch=20)
lines(cv.svsp.ffrl$start,cv.svsp.ffrl$do,col=alpha('blue',0.15),pch=20)

plot(co.svsp.dfr1$start,co.svsp.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(co.svsp.dfrh$start,co.svsp.dfrh$up,col=alpha('red',0.15),pch=20)
lines(co.svsp.dfrl$start,co.svsp.dfrl$do,col=alpha('red',0.15),pch=20)

plot(co.svsp.ffr1$start,co.svsp.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(co.svsp.ffrh$start,co.svsp.ffrh$up,col=alpha('blue',0.15),pch=20)
lines(co.svsp.ffrl$start,co.svsp.ffrl$do,col=alpha('blue',0.15),pch=20)

plot(cv.pla2.dfr1$start,cv.pla2.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(cv.pla2.dfrh$start,cv.pla2.dfrh$up,col=alpha('red',0.15),pch=20)
lines(cv.pla2.dfrl$start,cv.pla2.dfrl$do,col=alpha('red',0.15),pch=20)

plot(cv.pla2.ffr1$start,cv.pla2.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(cv.pla2.ffrh$start,cv.pla2.ffrh$up,col=alpha('blue',0.15),pch=20)
lines(cv.pla2.ffrl$start,cv.pla2.ffrl$do,col=alpha('blue',0.15),pch=20)

plot(co.pla2.dfr1$start,co.pla2.dfr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(co.pla2.dfrh$start,co.pla2.dfrh$up,col=alpha('red',0.15),pch=20)
lines(co.pla2.dfrl$start,co.pla2.dfrl$do,col=alpha('red',0.15),pch=20)

plot(co.pla2.ffr1$start,co.pla2.ffr1$mean,type='l',xlim=c(0,10000),ylim=c(0,0.8))
lines(co.pla2.ffrh$start,co.pla2.ffrh$up,col=alpha('blue',0.15),pch=20)
lines(co.pla2.ffrl$start,co.pla2.ffrl$do,col=alpha('blue',0.15),pch=20)
