#!/usr/bin/env Rscript
# usage: Rscript --vanilla MR.R $uid $exposure_file $outcome_file
## uid='30880_raw/30880_raw.male'
## exposure_file='30880_raw/30880_raw.male.5e8'
## outcome_file='common/male_pd_gsmr.raw'

# arguments
args = commandArgs(trailingOnly=TRUE)
uid=args[1]
exposure_file = args[2]
outcome_file = args[3]

# set up
library("TwoSampleMR")
library("data.table")

# load files
exposure_dat=read_exposure_data(filename=exposure_file, sep='\t')
exposure_dat <- clump_data(exposure_dat, clump_kb = 10000, clump_r2 = 0.001,clump_p1 = 1,clump_p2 = 1,pop = "EUR")
outcome_dat <- read_outcome_data(
    snps = exposure_dat$SNP,
    filename = outcome_file,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "b",
    se_col = "se",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    eaf_col = "freq",
    pval_col = "p",
    samplesize_col = "N"
    )

# harmonize
dat <- harmonise_data(exposure_dat, outcome_dat)

# do MR
res<-mr(dat)
het<-mr_heterogeneity(dat)
plt<-mr_pleiotropy_test(dat)
sin<-mr_singlesnp(dat)

# write data and results
fwrite(dat, paste0(uid, '.mr_dat.csv'))
fwrite(res, paste0(uid, '.mr_res.csv'))
fwrite(sin, paste0(uid, '.mr_sin.csv'))
fwrite(het, paste0(uid, '.mr_het.csv'))
fwrite(plt, paste0(uid, '.mr_plt.csv'))

# To retrieve
# dat = fread(paste0(uid, '.mr_dat.csv'))
# res = fread(paste0(uid, '.mr_res.csv'))
# sin = fread(paste0(uid, '.mr_sin.csv'))
# het = fread(paste0(uid, '.mr_het.csv'))
# plt = fread(paste0(uid, '.mr_plt.csv'))
