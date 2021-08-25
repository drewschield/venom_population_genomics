############################################################################
# Compare ß estimates from series of parameter settings in BetaScan
############################################################################

### Goal: evaluate whether ß estimates are robust to various -w and -p settings
### in BetaScan analysis.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & load dependencies-------------------------------

setwd('./beta/permutations/')

install.packages("rstatix")
library(rstatix)
install.packages("ggpubr")
library(ggpubr)

### Read in data------------------------------------------------------------

beta.perm <- read.table('combined.cv1.phased.permutations.betascores.txt',header=T)

### Compute correlation matrix for all variables----------------------------

beta.perm.prune <- beta.perm[, c(2,3,4,5,6,7,8,9,10,11,12,13)]
head(beta.perm.prune, 5)

# Set correlation matrix
beta.cor_mat <- beta.perm.prune %>% cor_mat()
# Get p-values
beta.cor_mat %>% cor_get_pval()

write.table(beta.cor_mat, "correlation_table.beta.txt",quote=FALSE,row.names=FALSE,sep='\t')
write.table(beta.cor_mat %>% cor_get_pval(), "correlation_table_p-val.beta.txt",quote=FALSE,row.names=FALSE,sep='\t')


cor.test(beta.perm$beta_p2_w500,beta.perm$beta_p2_w1000,method='spearman')
plot(beta.perm$beta_p2_w500,beta.perm$beta_p2_w1000,pch=20)
hist(beta.perm$beta_p2_w500,breaks=1000)
