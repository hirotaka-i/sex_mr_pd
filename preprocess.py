import sys
import os
import pandas as pd
import numpy as np
import argparse
import subprocess

"""
This fuction converts ukbb sumstats to ldsc/mtag format
1. reduce variants by maf (ldsc's default is maf>0.01)
2. use annovar to add rsid
3. add a header (-->for mtag)
4. munge the data (-->for ldsc)
Resource: mem 10G, 1 thread, 2hr would be enough
Need annovar/2019-10-24 and ldsc/3d0c4464
In addition, it also creates some files can be used for 
several MR softwares
"""

parser = argparse.ArgumentParser(description='Arguments for preprocessing summarystatistics')
parser.add_argument('--input', type=str, default='nope', help='Path to the input summary statistics eg. 30880_raw/30880_raw.gwas.imputed_v3.male.tsv.bgz')
parser.add_argument('--uid', type=str, default='nope', help='The unique id for data eg.30880_raw/30880_raw.male')
parser.add_argument('--nthread', type=str, default='1', help='Number of threads to use. default=1')
args = parser.parse_args()

# arguments
input_file=args.input
uid=args.uid
nthread=args.nthread

# parameter
variant_file = '/data/CARD/projects/sexMRforPD/common/variants.tsv.bgz'
w_hm3_snplist="/data/CARD/projects/sexMRforPD/common/mtag/ld_ref_panel/eur_w_ld_chr/w_hm3.snplist"


def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))


# add RSID
d = pd.read_csv(input_file, compression='gzip', sep='\t', engine='c', low_memory=False)
v = pd.read_csv(variant_file, compression='gzip', sep='\t', engine='c', low_memory=False)
print(d.shape, 'original')

if sum(d.variant!=v.variant)==0: # order check
    d['snpid']=v.rsid
    d['chr']=v.chr
    d['pos']=v.pos
    d['ref']=v.ref
    d['alt']=v.alt
    d['freq']=d.minor_AF.where(d.alt==d.minor_allele, other=1-d.minor_AF)
    d['N'] = d.n_complete_samples
else:
    print('UKBB varinats and sumstats orders are different')

# Filter Reduce maf > 0.01, low_confidence_variant
d = d[d.minor_AF>0.01].copy()
d = d[~d.low_confidence_variant].copy()
print(d.shape, 'filtered')


# mtag format for ldsc
t = d.rename(columns={'alt':'a1', 'ref':'a2', 'N':'n', 'pos':'bpos'}).copy()
t['z'] = t.beta/t.se
t[['snpid', 'chr', 'bpos', 'a1', 'a2', 'freq', 'beta', 'se', 'z', 'pval', 'n']].to_csv(f'{uid}.mtag', index=False, sep='\t')

# for TwoSampleMR
t = d.rename(columns={'snpid':'SNP', 'alt':'effect_allele', 'ref':'other_allele', 'freq':'eaf', 'pos':'position', 'N':'samplesize'}).copy()
t.loc[t.pval<5e-8, 
      ['SNP', 'beta', 'se', 'other_allele', 'effect_allele', 'pval', 'eaf', 'chr', 'position', 'samplesize'
      ]].to_csv(f'{uid}.5e8', sep='\t', index=False)

# for MR-LDP
t = d.rename(columns={'snpid':'SNP', 'pos':'BP', 'alt':'A1', 'ref':'A2','pval':'pvalue'}).copy()
t.loc[t.pvalue<1e-4, ['SNP', 'chr', 'BP', 'A2', 'A1', 'beta', 'se', 'pvalue']].to_csv(f'{uid}.1e4', sep='\t', index=False)


# for GSMR
t = d.rename(columns={'snpid':'SNP', 'alt':'A1', 'ref':'A2', 'beta':'b', 'pval':'p'}).copy()
t = t.drop_duplicates(subset='SNP',keep='first').copy()
print(t.shape,'GSMR')
t[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(f'{uid}.raw', index=False, sep='\t')

# munge for ldsc (ldsc/3d0c4464)
shell_do(f"munge_sumstats.py --sumstats {uid}.mtag --ignore beta --out {uid}.ldsc --merge-alleles {w_hm3_snplist} --chunksize 500000")
