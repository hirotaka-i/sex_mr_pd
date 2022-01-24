#!/usr/bin/env Rscript
# usage: Rscript --vanilla MRpower.R ${iid} ${harmonized_data}
args = commandArgs(trailingOnly=TRUE)
uid=args[1]
harmonized_data=args[2]

# Other parameter
OR = 1.2 # the OR for 1SD increase of exposure.
alpha = 0.05 # Type-I error rate


# set up
library("TwoSampleMR")

# functions
power_brion = function(R2xz, N, alpha=.05, K, OR, k = 1){
    # Brion's method
    ## online app://shiny. https:cnsgenomics.com/mRnd/
    ## paper: https://doi.org/10.1093/ije/dyt179

    ### R2xz: exposure variance explained by a snp
    ### N: total sample size
    ### alpha: alpha error
    ### K: proportion of cases in the study
    ### OR: Odds ratio for the outcome by 1 SD increase in exposure
    ### number of instrumental variants with the mean Rxy
    
    b01 = K * (OR / (1+ K * (OR - 1)) - 1)
    bMR = K * (OR / (1+ K * (OR - 1)) - 1)
    VbMR = (K * (1-K) - b01^2)/ (N * R2xz * k) # *k is added to adjust multi-instrumental MR
    
    NCP = bMR^2 / VbMR
    qthres = qchisq(1-alpha,df=1)
    power_brion = pchisq(qthres, df=1, ncp=NCP, lower.tail=F)
    return(power_brion)
}

power_burgess = function(R2xz, N, alpha=.05, K, OR, k=1){
    # Burgess method
    ## online app: https://sb452.shinyapps.io/power/
    ## paper: https://doi.org/10.1093/ije/dyu005
    b1 = log(OR)
    rxz = sqrt(R2xz * k) # *k is added to adjust multi-indtrumental MR
    zb = b1 * rxz * sqrt(N * K * (1-K))
    za = qnorm(1-alpha/2)
    power_burgess = pnorm(zb-za)
    return(power_burgess)
}

calculate_F = function(R2xz, n, k){ # k: number of instrument SNPs in the model
    F = R2xz * (n - 1 -k) / ((1 - R2xz) * k)
    return(F)
}



# load data
dat = read.csv(harmonized_data)
if(nrow(dat)==0){print('No Instrumental SNPs');q('y')}
n = mean(dat$samplesize.exposure)
k = nrow(dat)

# parameters for Male/Female GWAS
if(grepl('female',uid)){
    N = 7384 + 12389
    K = 7384 / N
}else{
    N = 12054 + 11999
    K = 12054 / N
}


# get F
# mean F
R2s = get_r_from_pn(dat$pval.exposure, dat$samplesize.exposure)^2
meanR2xz =  mean(R2s) # exposure variance explained by instruments (mean)
mean_F = calculate_F(meanR2xz, n, 1)
# max F
maxR2xz =  max(R2s) # exposure variance explained by instruments (mean)
max_F = calculate_F(maxR2xz, n, 1)

# all_F
all_F = calculate_F(sum(R2s), n, length(R2s)) 
#NOTE: F cannot be compared if the sample size and n_snps are different

df_power = data.frame(
    R2_used = c('mean', 'mean', 'max', 'max', 'all', 'all'),
    method = c('Brion', 'Burgess', 'Brion', 'Burgess','Brion', 'Burgess'),
    n_snp = c(1,1,1,1,k,k),
    power = c(
        power_brion(meanR2xz, N, alpha=.05, K, OR),
        power_burgess(meanR2xz, N, alpha=.05, K, OR),
        power_brion(maxR2xz, N, alpha=.05, K, OR),
        power_burgess(maxR2xz, N, alpha=.05, K, OR),
        power_brion(meanR2xz, N, alpha=.05, K, OR, k),
        power_burgess(meanR2xz, N, alpha=.05, K, OR, k)
    ),
    F = c(mean_F, mean_F, max_F, max_F, all_F, all_F)
)

write.csv(df_power, paste0(uid,'.mr_power.csv'),row.names=F)