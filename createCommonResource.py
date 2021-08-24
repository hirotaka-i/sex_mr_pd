# This is to process the original PD summary stats for our analysis.
# Also MR-LDP folder were downloaded [from this repo](https://drive.google.com/drive/folders/1IAs3daG9TIvjnneR32j1pfV0niz8PHyu)
wkdir='/data/CARD/projects/sexMRforPD'

import pandas as pd
import numpy as np

for sex in ['male', 'female']:
    print(sex)
    SEX=sex.upper()

    d = pd.read_csv(f"/data/LNG/CORNELIS_TEMP/MALE_FEMALE_PD/RESULTS/to_meta_files/{SEX}_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt", sep='\t', engine='c') # 1min
    d['FreqDiff']=d.MaxFreq - d.MinFreq
    print('original', d.shape)
    d = d[pd.notna(d.ID)&(d.FreqDiff<0.15)&(d.HetISq<80)&(d.Effect<5)&(d.StdErr<5)].copy()
    print('filtered - major', d.shape)
    
    # process flipping
    d['flip'] = [1 if ((a1.upper()==ref)&(a2.upper()==alt)) 
             else 0 if ((a1.upper()==alt)&(a2.upper()==ref))
             else np.nan for a1,a2,ref,alt in zip(d.Allele1, d.Allele2, d.REF, d.ALT)]

    print('flipping (1)\n', d.flip.value_counts(dropna=False))

    # remove ambiguous flipping and then create variables. 
    d = d[pd.notna(d.flip)].copy()
    d['A1'] = d.ALT
    d['A2'] = d.REF
    d['freq'] = d.Freq1.where(d.flip==0, other=1-d.Freq1)
    d['b'] = d.Effect.where(d.flip==0, other=-d.Effect)
    d['se'] = d.StdErr
    d['p'] = d['P-value']
    if sex=='male':
        d['N'] = 24053
    elif sex=='female':
        d['N'] = 19773
    d['SNP'] =d.ID
    
    
    
    dt = d.loc[(d['P-value']<5e-8), 'ID'].copy()
    dt.to_csv(f'{wkdir}/common/{sex}_pd_iv_5e8.txt', index=False, header=False)
    print(dt.shape)
    dt = d.loc[(d['P-value']<1e-4), 'ID'].copy()
    dt.to_csv(f'{wkdir}/common/{sex}_pd_iv_1e4.txt', index=False, header=False)
    print(dt.shape)

    # gsmr
    d[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(f'{wkdir}/common/{sex}_pd_gsmr.raw', sep='\t', index=False)

    # For MR-LDP
    t = d['MarkerName'].str.split(':', expand=True)
    d['chr']=t[0]
    d['BP']=t[1]
    d['beta']=d.b
    d['pvalue']=d.p
    d[['SNP', 'chr', 'BP', 'A1', 'A2', 'beta', 'se', 'pvalue']].to_csv(f'{wkdir}/common/{sex}_pd_mrldp.raw', sep='\t', index=False)

    
## Output
# male
# original (7234141, 19)
# filtered - major (7106595, 19)
# flipping (1)
#  0    3610445
# 1    3496150
# Name: flip, dtype: int64
# (604,)
# (4007,)
# female
# original (7309889, 19)
# filtered - major (7182330, 19)
# flipping (1)
#  0    3625649
# 1    3556681
# Name: flip, dtype: int64
# (635,)
# (3628,)