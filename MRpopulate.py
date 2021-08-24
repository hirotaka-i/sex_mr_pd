# populate results
import sys
import os
import pandas as pd
import numpy as np
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Arguments for preprocessing summarystatistics')
parser.add_argument('--pheno', type=str, default='nope', help='unique id for data eg.30880_raw/30880_raw.male')
args = parser.parse_args()

cols = ['pheno_ukbb', 'sex', 'n_snp', 'n_exp', 'n_out', 
        'p_IVWT', 'p_Egger','p_plt', 'het_flag', 'mean_F','power_max1', 'power_IVWT', 
        'b_IVWT', 'se_IVWT', 'b_Egger', 'se_Egger']

pheno=args.pheno

d = pd.DataFrame(data={'pheno_ukbb':pheno, 'sex':['male', 'female']}, columns=cols)

for sex in ['male', 'female']:
    file = f'{pheno}/{pheno}.{sex}.mr_dat.csv'
    if os.path.exists(file):
        dat = pd.read_csv(file)
        d.loc[d.sex==sex, ['n_exp']] = np.mean(dat['samplesize.exposure'])
        d.loc[d.sex==sex, ['n_out']] = np.mean(dat['samplesize.outcome'])

    file = f'{pheno}/{pheno}_MR_{sex}_res.csv'
    if os.path.exists(file):
        res = pd.read_csv(file)
        d.loc[d.sex==sex, ['n_snp', 'b_IVWT', 'se_IVWT','p_IVWT']] = res.loc[res.method=='Inverse variance weighted',['nsnp','b', 'se', 'pval']].values
        d.loc[d.sex==sex, ['b_Egger', 'se_Egger','p_Egger']] = res.loc[res.method=='MR Egger',['b', 'se', 'pval']].values
    
    file = f'{pheno}/{pheno}.{sex}.mr_het.csv'
    if os.path.exists(file):
        het = pd.read_csv(file)
        d.loc[d.sex==sex, ['het_flag']]= (',').join(het.method[het.Q_pval<0.05])

    file = f'{pheno}/{pheno}.{sex}.mr_plt.csv'
    if os.path.exists(file):
        plt = pd.read_csv(file)
        d.loc[d.sex==sex, ['p_plt']]= plt.pval[0]

    file = f'{pheno}/{pheno}.{sex}.mr_power.csv'
    if os.path.exists(file):
        power = pd.read_csv(file)
        d.loc[d.sex==sex, ['mean_F']] = power[(power.R2_used=='mean')&(power.method=='Burgess')&(power.n_snp==1)].F.values
        d.loc[d.sex==sex, ['power_max1']] = power[(power.R2_used=='max')&(power.method=='Burgess')&(power.n_snp==1)].power.values
        d.loc[d.sex==sex, ['power_IVWT']] = power[(power.R2_used=='mean')&(power.method=='Burgess')].drop_duplicates(subset='method',keep='last').power.values

d.to_csv(f'{pheno}/{pheno}_MR.csv', index=False)

print(f'{pheno}: Results were populated')