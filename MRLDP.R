#!/usr/bin/env Rscript (4.1.0)
# usage: Rscript --vanilla MRLDP.R $uid $filescreen $fileexposure $fileoutcome $stringname3 $blockfile $nthread
# uid='30880_raw/30880_raw.male'
# filescreen='30880_raw/30880_raw.male.1e4'
# fileexposure='30880_raw/30880_raw.male.1e4'
# fileoutcome='common/male_pd_mrldp.raw'
# stringname3='common/MR-LDP/all_chr_1000G'
# blockfile='common/MR-LDP/fourier_ls-all.bed'
# nthread=2
# stringname3='/data/CARD/GENERAL/1000g_p3/euro_b/1000g_euro' # This can be also used?
#### example from the paper ########################################
# # filescreen='common/MR-LDP/heart attack_myocardial infarction.txt'
# # fileexposure='common/MR-LDP/c4d.txt'
# # fileoutcome='common/MR-LDP/cardiogram.txt'
######################################################################

# arguments
args = commandArgs(trailingOnly=TRUE)
uid=args[1]
filescreen = args[2]
fileexposure = args[3]
fileoutcome=args[4]
stringname3=args[5]
blockfile=args[6]
nthread=as.numeric(args[7])

# set up
library("MR.LDP")

# parameter
pva_cutoff = 1e-4
scrres = matchscreen(filescreen, fileexposure, fileoutcome, stringname3, pva_cutoff, matchExp = FALSE )
bh1 = as.numeric ( scrres$bh1 )
bh2 = as.numeric ( scrres$bh2 )
s12 = as.numeric ( scrres$s12 )
s22 = as.numeric ( scrres$s22 )
chr = as.numeric ( scrres$chr )
bp = scrres$bp
rsname = scrres$rsname
avbIndex = scrres$idxin
idx4panel = scrres$idx4panel


# remove MHC regions
QCindex = 1;
if( QCindex ) {
    QCresult = summaryQC ( mhcstart , mhcend , bh1 , bh2 , s12 , s22 , bp , chr , rsname , avbIndex , idx4panel , Inf , Inf )
    bh1new = QCresult$bh1new
    bh2new = QCresult$bh2new
    s12new = QCresult$s12new
    s22new = QCresult$s22new
    bpnew = QCresult$bpnew
    chrnew = QCresult$chrnew
    avbIndexnew = QCresult$avbIndexnew
    idx4panelnew = QCresult$idx4panel
    rsnamenew = QCresult$rsnamenew ;
    } else {
    bh1new = bh1
    bh2new = bh2
    s12new = s12
    s22new = s22
    bpnew = bp
    chrnew = chr
    rsnamenew = rsname
    idx4panelnew = idx4panel
    avbIndexnew = avbIndex
}
p = length ( avbIndexnew )





# fit MR-LDP with various lambda
# initialize
gamma = rep (0.01 , p )
alpha = rep (0.01 , p )
sgga2 = 0.01
sgal2 = 0.01
beta0 = 0
maxIter = 3000
coreNum = nthread  ### Number of cores!
epsStopLogLik = 1e-7
DD = data.frame(nsnp=NA, lam=NA, b=NA, se=NA, pval=NA, Tstat_MRLDP=NA, iter_Hb=NA, iter_H0=NA)

for (i in 1:7){
    print(paste('iter', i, '/ 7'))
    # lambda
    lam=0.5+0.05*(i)

    RealMRLDP_Hb = MRLDP_RealPXvb_block ( 
        bpnew , chrnew , avbIndexnew-1 , idx4panelnew , blockfile , stringname3 , 
        bh1new , bh2new , s12new , s22new , 
        gamma , alpha , beta0 , sgga2 , sgal2 , 
        coreNum , lam , 0 , epsStopLogLik , maxIter , model = 2)
    RealMRLDP_H0 = MRLDP_RealPXvb_block ( 
        bpnew , chrnew , avbIndex-1 , idx4panelnew , blockfile , stringname3 , 
        bh1new , bh2new , s12new , s22new , 
        gamma , alpha , beta0 , sgga2 , sgal2 , 
        coreNum , lam , 1 , epsStopLogLik , maxIter , model = 2)
    beta0_MRLDP = RealMRLDP_Hb$beta0
    Tstat_MRLDP <- 2*( RealMRLDP_Hb$tstat - RealMRLDP_H0$tstat ) 
    MRLDP_se = abs( RealMRLDP_Hb$beta0 / sqrt (Tstat_MRLDP))
    pval = pnorm(-abs(beta0_MRLDP/MRLDP_se))*2
    DD[i,] = c(p, lam, beta0_MRLDP, MRLDP_se, pval, Tstat_MRLDP, RealMRLDP_Hb$Iteration, RealMRLDP_H0$Iteration)
    print(DD[i,])
}

print(DD)

write.csv(DD, paste0(uid, '.mrldp.csv'), row.names = FALSE)