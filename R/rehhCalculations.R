############################################################################
# Genomic scans haplotype diversity statistics
############################################################################

### Goal: look for evidence of selection based on extended haplotype lengths
### using the 'rehh' package and phased haplotypes for C. viridis and C.
### oreganus from the recombination study.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Set working directory & dependencies------------------------------------

setwd('/Volumes/GoogleDrive/My Drive/projects/venom_population_genomics/3_selection/rehh/')

install.packages('rehh')
install.packages('tidyselect')
install.packages('tidyverse')
install.packages('vcfR')
library(rehh)
library(vcfR)
library(tidyverse)
library(data.table)

### Read in data------------------------------------------------------------

cvhh.ma1 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-ma1.vcf", polarize_vcf = FALSE)
cvhh.ma2 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-ma2.vcf", polarize_vcf = FALSE)
cvhh.ma3 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-ma3.vcf", polarize_vcf = FALSE)
cvhh.ma4 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-ma4.vcf", polarize_vcf = FALSE)
cvhh.ma5 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-ma5.vcf", polarize_vcf = FALSE)
cvhh.ma6 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-ma6.vcf", polarize_vcf = FALSE)
cvhh.ma7 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-ma7.vcf", polarize_vcf = FALSE)
cvhh.Z <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-Z.vcf", polarize_vcf = FALSE)
cvhh.mi1 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi1.vcf", polarize_vcf = FALSE)
cvhh.mi2 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi2.vcf", polarize_vcf = FALSE)
cvhh.mi3 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi3.vcf", polarize_vcf = FALSE)
cvhh.mi4 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi4.vcf", polarize_vcf = FALSE)
cvhh.mi5 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi5.vcf", polarize_vcf = FALSE)
cvhh.mi6 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi6.vcf", polarize_vcf = FALSE)
cvhh.mi7 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi7.vcf", polarize_vcf = FALSE)
cvhh.mi8 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi8.vcf", polarize_vcf = FALSE)
cvhh.mi9 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi9.vcf", polarize_vcf = FALSE)
cvhh.mi10 <- data2haplohh(hap_file = "./data/viridis.phased.scaffold-mi10.vcf", polarize_vcf = FALSE)

cohh.ma1 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-ma1.vcf", polarize_vcf = FALSE)
cohh.ma2 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-ma2.vcf", polarize_vcf = FALSE)
cohh.ma3 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-ma3.vcf", polarize_vcf = FALSE)
cohh.ma4 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-ma4.vcf", polarize_vcf = FALSE)
cohh.ma5 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-ma5.vcf", polarize_vcf = FALSE)
cohh.ma6 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-ma6.vcf", polarize_vcf = FALSE)
cohh.ma7 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-ma7.vcf", polarize_vcf = FALSE)
cohh.Z <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-Z.vcf", polarize_vcf = FALSE)
cohh.mi1 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi1.vcf", polarize_vcf = FALSE)
cohh.mi2 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi2.vcf", polarize_vcf = FALSE)
cohh.mi3 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi3.vcf", polarize_vcf = FALSE)
cohh.mi4 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi4.vcf", polarize_vcf = FALSE)
cohh.mi5 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi5.vcf", polarize_vcf = FALSE)
cohh.mi6 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi6.vcf", polarize_vcf = FALSE)
cohh.mi7 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi7.vcf", polarize_vcf = FALSE)
cohh.mi8 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi8.vcf", polarize_vcf = FALSE)
cohh.mi9 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi9.vcf", polarize_vcf = FALSE)
cohh.mi10 <- data2haplohh(hap_file = "./data/oreganus.phased.scaffold-mi10.vcf", polarize_vcf = FALSE)

### Filter on MAF-----------------------------------------------------------

cvhh.ma1_f <- subset(cvhh.ma1, min_maf = 0.05)
cvhh.ma2_f <- subset(cvhh.ma2, min_maf = 0.05)
cvhh.ma3_f <- subset(cvhh.ma3, min_maf = 0.05)
cvhh.ma4_f <- subset(cvhh.ma4, min_maf = 0.05)
cvhh.ma5_f <- subset(cvhh.ma5, min_maf = 0.05)
cvhh.ma6_f <- subset(cvhh.ma6, min_maf = 0.05)
cvhh.ma7_f <- subset(cvhh.ma7, min_maf = 0.05)
cvhh.Z_f <- subset(cvhh.Z, min_maf = 0.05)
cvhh.mi1_f <- subset(cvhh.mi1, min_maf = 0.05)
cvhh.mi2_f <- subset(cvhh.mi2, min_maf = 0.05)
cvhh.mi3_f <- subset(cvhh.mi3, min_maf = 0.05)
cvhh.mi4_f <- subset(cvhh.mi4, min_maf = 0.05)
cvhh.mi5_f <- subset(cvhh.mi5, min_maf = 0.05)
cvhh.mi6_f <- subset(cvhh.mi6, min_maf = 0.05)
cvhh.mi7_f <- subset(cvhh.mi7, min_maf = 0.05)
cvhh.mi8_f <- subset(cvhh.mi8, min_maf = 0.05)
cvhh.mi9_f <- subset(cvhh.mi9, min_maf = 0.05)
cvhh.mi10_f <- subset(cvhh.mi10, min_maf = 0.05)

cohh.ma1_f <- subset(cohh.ma1, min_maf = 0.05)
cohh.ma2_f <- subset(cohh.ma2, min_maf = 0.05)
cohh.ma3_f <- subset(cohh.ma3, min_maf = 0.05)
cohh.ma4_f <- subset(cohh.ma4, min_maf = 0.05)
cohh.ma5_f <- subset(cohh.ma5, min_maf = 0.05)
cohh.ma6_f <- subset(cohh.ma6, min_maf = 0.05)
cohh.ma7_f <- subset(cohh.ma7, min_maf = 0.05)
cohh.Z_f <- subset(cohh.Z, min_maf = 0.05)
cohh.mi1_f <- subset(cohh.mi1, min_maf = 0.05)
cohh.mi2_f <- subset(cohh.mi2, min_maf = 0.05)
cohh.mi3_f <- subset(cohh.mi3, min_maf = 0.05)
cohh.mi4_f <- subset(cohh.mi4, min_maf = 0.05)
cohh.mi5_f <- subset(cohh.mi5, min_maf = 0.05)
cohh.mi6_f <- subset(cohh.mi6, min_maf = 0.05)
cohh.mi7_f <- subset(cohh.mi7, min_maf = 0.05)
cohh.mi8_f <- subset(cohh.mi8, min_maf = 0.05)
cohh.mi9_f <- subset(cohh.mi9, min_maf = 0.05)
cohh.mi10_f <- subset(cohh.mi10, min_maf = 0.05)

### Perform scans-----------------------------------------------------------

# Fix: replace inputs with '_f' appended to designated being filtered
cv.ma1_scan <- scan_hh(cvhh.ma1_f, polarized = FALSE)
cv.ma2_scan <- scan_hh(cvhh.ma2_f, polarized = FALSE)
cv.ma3_scan <- scan_hh(cvhh.ma3_f, polarized = FALSE)
cv.ma4_scan <- scan_hh(cvhh.ma4_f, polarized = FALSE)
cv.ma5_scan <- scan_hh(cvhh.ma5_f, polarized = FALSE)
cv.ma6_scan <- scan_hh(cvhh.ma6_f, polarized = FALSE)
cv.ma7_scan <- scan_hh(cvhh.ma7_f, polarized = FALSE)
cv.Z_scan <- scan_hh(cvhh.Z_f, polarized = FALSE)

co.ma1_scan <- scan_hh(cohh.ma1_f, polarized = FALSE)
co.ma2_scan <- scan_hh(cohh.ma2_f, polarized = FALSE)
co.ma3_scan <- scan_hh(cohh.ma3_f, polarized = FALSE)
co.ma4_scan <- scan_hh(cohh.ma4_f, polarized = FALSE)
co.ma5_scan <- scan_hh(cohh.ma5_f, polarized = FALSE)
co.ma6_scan <- scan_hh(cohh.ma6_f, polarized = FALSE)
co.ma7_scan <- scan_hh(cohh.ma7_f, polarized = FALSE)
co.Z_scan <- scan_hh(cohh.Z_f, polarized = FALSE)

cv.mi1_scan <- scan_hh(cvhh.mi1_f, polarized = FALSE)
cv.mi2_scan <- scan_hh(cvhh.mi2_f, polarized = FALSE)
cv.mi3_scan <- scan_hh(cvhh.mi3_f, polarized = FALSE)
cv.mi4_scan <- scan_hh(cvhh.mi4_f, polarized = FALSE)
cv.mi5_scan <- scan_hh(cvhh.mi5_f, polarized = FALSE)
cv.mi6_scan <- scan_hh(cvhh.mi6_f, polarized = FALSE)
cv.mi7_scan <- scan_hh(cvhh.mi7_f, polarized = FALSE)
cv.mi8_scan <- scan_hh(cvhh.mi8_f, polarized = FALSE)
cv.mi9_scan <- scan_hh(cvhh.mi9_f, polarized = FALSE)
cv.mi10_scan <- scan_hh(cvhh.mi10_f, polarized = FALSE)

co.mi1_scan <- scan_hh(cohh.mi1_f, polarized = FALSE)
co.mi2_scan <- scan_hh(cohh.mi2_f, polarized = FALSE)
co.mi3_scan <- scan_hh(cohh.mi3_f, polarized = FALSE)
co.mi4_scan <- scan_hh(cohh.mi4_f, polarized = FALSE)
co.mi5_scan <- scan_hh(cohh.mi5_f, polarized = FALSE)
co.mi6_scan <- scan_hh(cohh.mi6_f, polarized = FALSE)
co.mi7_scan <- scan_hh(cohh.mi7_f, polarized = FALSE)
co.mi8_scan <- scan_hh(cohh.mi8_f, polarized = FALSE)
co.mi9_scan <- scan_hh(cohh.mi9_f, polarized = FALSE)
co.mi10_scan <- scan_hh(cohh.mi10_f, polarized = FALSE)

# Older
cv.mi1_scan <- scan_hh(cvhh.mi1_f, polarized = FALSE)
co.mi1_scan <- scan_hh(cohh.mi1_f, polarized = FALSE)

### Calculate iHS-----------------------------------------------------------

cv.ma1_ihs <- ihh2ihs(cv.ma1_scan, freqbin = 1)
cv.ma2_ihs <- ihh2ihs(cv.ma2_scan, freqbin = 1)
cv.ma3_ihs <- ihh2ihs(cv.ma3_scan, freqbin = 1)
cv.ma4_ihs <- ihh2ihs(cv.ma4_scan, freqbin = 1)
cv.ma5_ihs <- ihh2ihs(cv.ma5_scan, freqbin = 1)
cv.ma6_ihs <- ihh2ihs(cv.ma6_scan, freqbin = 1)
cv.ma7_ihs <- ihh2ihs(cv.ma7_scan, freqbin = 1)
cv.Z_ihs <- ihh2ihs(cv.Z_scan, freqbin = 1)
cv.mi1_ihs <- ihh2ihs(cv.mi1_scan, freqbin = 1)
cv.mi2_ihs <- ihh2ihs(cv.mi2_scan, freqbin = 1)
cv.mi3_ihs <- ihh2ihs(cv.mi3_scan, freqbin = 1)
cv.mi4_ihs <- ihh2ihs(cv.mi4_scan, freqbin = 1)
cv.mi5_ihs <- ihh2ihs(cv.mi5_scan, freqbin = 1)
cv.mi6_ihs <- ihh2ihs(cv.mi6_scan, freqbin = 1)
cv.mi7_ihs <- ihh2ihs(cv.mi7_scan, freqbin = 1)
cv.mi8_ihs <- ihh2ihs(cv.mi8_scan, freqbin = 1)
cv.mi9_ihs <- ihh2ihs(cv.mi9_scan, freqbin = 1)
cv.mi10_ihs <- ihh2ihs(cv.mi10_scan, freqbin = 1)

co.ma1_ihs <- ihh2ihs(co.ma1_scan, freqbin = 1)
co.ma2_ihs <- ihh2ihs(co.ma2_scan, freqbin = 1)
co.ma3_ihs <- ihh2ihs(co.ma3_scan, freqbin = 1)
co.ma4_ihs <- ihh2ihs(co.ma4_scan, freqbin = 1)
co.ma5_ihs <- ihh2ihs(co.ma5_scan, freqbin = 1)
co.ma6_ihs <- ihh2ihs(co.ma6_scan, freqbin = 1)
co.ma7_ihs <- ihh2ihs(co.ma7_scan, freqbin = 1)
co.Z_ihs <- ihh2ihs(co.Z_scan, freqbin = 1)
co.mi1_ihs <- ihh2ihs(co.mi1_scan, freqbin = 1)
co.mi2_ihs <- ihh2ihs(co.mi2_scan, freqbin = 1)
co.mi3_ihs <- ihh2ihs(co.mi3_scan, freqbin = 1)
co.mi4_ihs <- ihh2ihs(co.mi4_scan, freqbin = 1)
co.mi5_ihs <- ihh2ihs(co.mi5_scan, freqbin = 1)
co.mi6_ihs <- ihh2ihs(co.mi6_scan, freqbin = 1)
co.mi7_ihs <- ihh2ihs(co.mi7_scan, freqbin = 1)
co.mi8_ihs <- ihh2ihs(co.mi8_scan, freqbin = 1)
co.mi9_ihs <- ihh2ihs(co.mi9_scan, freqbin = 1)
co.mi10_ihs <- ihh2ihs(co.mi10_scan, freqbin = 1)

### Plot iHS statistics-----------------------------------------------------

# Statistic
ggplot(cv.mi1_ihs$ihs, aes(POSITION, IHS)) + geom_point()
ggplot(co.mi1_ihs$ihs, aes(POSITION, IHS)) + geom_point()

# P-value
ggplot(cv.mi1_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()
ggplot(co.mi1_ihs$ihs, aes(POSITION, LOGPVALUE)) + geom_point()

### Write results-----------------------------------------------------------

cv.ma1_ihs <- tbl_df(cv.ma1_ihs$ihs)
colnames(cv.ma1_ihs) <- tolower(colnames(cv.ma1_ihs))
write_tsv(cv.ma1_ihs, "./cv.ma1_ihs.txt")

cv.ma2_ihs <- tbl_df(cv.ma2_ihs$ihs)
colnames(cv.ma2_ihs) <- tolower(colnames(cv.ma2_ihs))
write_tsv(cv.ma2_ihs, "./cv.ma2_ihs.txt")

cv.ma3_ihs <- tbl_df(cv.ma3_ihs$ihs)
colnames(cv.ma3_ihs) <- tolower(colnames(cv.ma3_ihs))
write_tsv(cv.ma3_ihs, "./cv.ma3_ihs.txt")

cv.ma4_ihs <- tbl_df(cv.ma4_ihs$ihs)
colnames(cv.ma4_ihs) <- tolower(colnames(cv.ma4_ihs))
write_tsv(cv.ma4_ihs, "./cv.ma4_ihs.txt")

cv.ma5_ihs <- tbl_df(cv.ma5_ihs$ihs)
colnames(cv.ma5_ihs) <- tolower(colnames(cv.ma5_ihs))
write_tsv(cv.ma5_ihs, "./cv.ma5_ihs.txt")

cv.ma6_ihs <- tbl_df(cv.ma6_ihs$ihs)
colnames(cv.ma6_ihs) <- tolower(colnames(cv.ma6_ihs))
write_tsv(cv.ma6_ihs, "./cv.ma6_ihs.txt")

cv.ma7_ihs <- tbl_df(cv.ma7_ihs$ihs)
colnames(cv.ma7_ihs) <- tolower(colnames(cv.ma7_ihs))
write_tsv(cv.ma7_ihs, "./cv.ma7_ihs.txt")

cv.Z_ihs <- tbl_df(cv.Z_ihs$ihs)
colnames(cv.Z_ihs) <- tolower(colnames(cv.Z_ihs))
write_tsv(cv.Z_ihs, "./cv.Z_ihs.txt")

cv.mi1_ihs <- tbl_df(cv.mi1_ihs$ihs)
colnames(cv.mi1_ihs) <- tolower(colnames(cv.mi1_ihs))
write_tsv(cv.mi1_ihs, "./cv.mi1_ihs.txt")

cv.mi2_ihs <- tbl_df(cv.mi2_ihs$ihs)
colnames(cv.mi2_ihs) <- tolower(colnames(cv.mi2_ihs))
write_tsv(cv.mi2_ihs, "./cv.mi2_ihs.txt")

cv.mi3_ihs <- tbl_df(cv.mi3_ihs$ihs)
colnames(cv.mi3_ihs) <- tolower(colnames(cv.mi3_ihs))
write_tsv(cv.mi3_ihs, "./cv.mi3_ihs.txt")

cv.mi4_ihs <- tbl_df(cv.mi4_ihs$ihs)
colnames(cv.mi4_ihs) <- tolower(colnames(cv.mi4_ihs))
write_tsv(cv.mi4_ihs, "./cv.mi4_ihs.txt")

cv.mi5_ihs <- tbl_df(cv.mi5_ihs$ihs)
colnames(cv.mi5_ihs) <- tolower(colnames(cv.mi5_ihs))
write_tsv(cv.mi5_ihs, "./cv.mi5_ihs.txt")

cv.mi6_ihs <- tbl_df(cv.mi6_ihs$ihs)
colnames(cv.mi6_ihs) <- tolower(colnames(cv.mi6_ihs))
write_tsv(cv.mi6_ihs, "./cv.mi6_ihs.txt")

cv.mi7_ihs <- tbl_df(cv.mi7_ihs$ihs)
colnames(cv.mi7_ihs) <- tolower(colnames(cv.mi7_ihs))
write_tsv(cv.mi7_ihs, "./cv.mi7_ihs.txt")

cv.mi8_ihs <- tbl_df(cv.mi8_ihs$ihs)
colnames(cv.mi8_ihs) <- tolower(colnames(cv.mi8_ihs))
write_tsv(cv.mi8_ihs, "./cv.mi8_ihs.txt")

cv.mi9_ihs <- tbl_df(cv.mi9_ihs$ihs)
colnames(cv.mi9_ihs) <- tolower(colnames(cv.mi9_ihs))
write_tsv(cv.mi9_ihs, "./cv.mi9_ihs.txt")

cv.mi10_ihs <- tbl_df(cv.mi10_ihs$ihs)
colnames(cv.mi10_ihs) <- tolower(colnames(cv.mi10_ihs))
write_tsv(cv.mi10_ihs, "./cv.mi10_ihs.txt")


co.ma1_ihs <- tbl_df(co.ma1_ihs$ihs)
colnames(co.ma1_ihs) <- tolower(colnames(co.ma1_ihs))
write_tsv(co.ma1_ihs, "./co.ma1_ihs.txt")

co.ma2_ihs <- tbl_df(co.ma2_ihs$ihs)
colnames(co.ma2_ihs) <- tolower(colnames(co.ma2_ihs))
write_tsv(co.ma2_ihs, "./co.ma2_ihs.txt")

co.ma3_ihs <- tbl_df(co.ma3_ihs$ihs)
colnames(co.ma3_ihs) <- tolower(colnames(co.ma3_ihs))
write_tsv(co.ma3_ihs, "./co.ma3_ihs.txt")

co.ma4_ihs <- tbl_df(co.ma4_ihs$ihs)
colnames(co.ma4_ihs) <- tolower(colnames(co.ma4_ihs))
write_tsv(co.ma4_ihs, "./co.ma4_ihs.txt")

co.ma5_ihs <- tbl_df(co.ma5_ihs$ihs)
colnames(co.ma5_ihs) <- tolower(colnames(co.ma5_ihs))
write_tsv(co.ma5_ihs, "./co.ma5_ihs.txt")

co.ma6_ihs <- tbl_df(co.ma6_ihs$ihs)
colnames(co.ma6_ihs) <- tolower(colnames(co.ma6_ihs))
write_tsv(co.ma6_ihs, "./co.ma6_ihs.txt")

co.ma7_ihs <- tbl_df(co.ma7_ihs$ihs)
colnames(co.ma7_ihs) <- tolower(colnames(co.ma7_ihs))
write_tsv(co.ma7_ihs, "./co.ma7_ihs.txt")

co.Z_ihs <- tbl_df(co.Z_ihs$ihs)
colnames(co.Z_ihs) <- tolower(colnames(co.Z_ihs))
write_tsv(co.Z_ihs, "./co.Z_ihs.txt")

co.mi1_ihs <- tbl_df(co.mi1_ihs$ihs)
colnames(co.mi1_ihs) <- tolower(colnames(co.mi1_ihs))
write_tsv(co.mi1_ihs, "./co.mi1_ihs.txt")

co.mi2_ihs <- tbl_df(co.mi2_ihs$ihs)
colnames(co.mi2_ihs) <- tolower(colnames(co.mi2_ihs))
write_tsv(co.mi2_ihs, "./co.mi2_ihs.txt")

co.mi3_ihs <- tbl_df(co.mi3_ihs$ihs)
colnames(co.mi3_ihs) <- tolower(colnames(co.mi3_ihs))
write_tsv(co.mi3_ihs, "./co.mi3_ihs.txt")

co.mi4_ihs <- tbl_df(co.mi4_ihs$ihs)
colnames(co.mi4_ihs) <- tolower(colnames(co.mi4_ihs))
write_tsv(co.mi4_ihs, "./co.mi4_ihs.txt")

co.mi5_ihs <- tbl_df(co.mi5_ihs$ihs)
colnames(co.mi5_ihs) <- tolower(colnames(co.mi5_ihs))
write_tsv(co.mi5_ihs, "./co.mi5_ihs.txt")

co.mi6_ihs <- tbl_df(co.mi6_ihs$ihs)
colnames(co.mi6_ihs) <- tolower(colnames(co.mi6_ihs))
write_tsv(co.mi6_ihs, "./co.mi6_ihs.txt")

co.mi7_ihs <- tbl_df(co.mi7_ihs$ihs)
colnames(co.mi7_ihs) <- tolower(colnames(co.mi7_ihs))
write_tsv(co.mi7_ihs, "./co.mi7_ihs.txt")

co.mi8_ihs <- tbl_df(co.mi8_ihs$ihs)
colnames(co.mi8_ihs) <- tolower(colnames(co.mi8_ihs))
write_tsv(co.mi8_ihs, "./co.mi8_ihs.txt")

co.mi9_ihs <- tbl_df(co.mi9_ihs$ihs)
colnames(co.mi9_ihs) <- tolower(colnames(co.mi9_ihs))
write_tsv(co.mi9_ihs, "./co.mi9_ihs.txt")

co.mi10_ihs <- tbl_df(co.mi10_ihs$ihs)
colnames(co.mi10_ihs) <- tolower(colnames(co.mi10_ihs))
write_tsv(co.mi10_ihs, "./co.mi10_ihs.txt")
