import sys
import os
import argparse
import subprocess

parser = argparse.ArgumentParser(description='Arguments for preprocessing summarystatistics')
parser.add_argument('--uid', type=str, default='nope', help='unique id for data eg.30880_raw/30880_raw.male')
parser.add_argument('--bfile', type=str, default='nope', help='Path to the reference plink binary eg. /data/CARD/GENERAL/1000g_p3/euro_b/1000g_euro')
parser.add_argument('--exposure', type=str, default='nope', help='Path to the exposure summary statistics eg. 30880_raw/30880_raw.male.raw')
parser.add_argument('--outcome', type=str, default='nope', help='Path to the outcome summary statistics eg. common/male_pd_gsmr.raw')
args = parser.parse_args()

# arguments
uid=args.uid
bfile=args.bfile
exposure_file=args.exposure
outcome_file=args.outcome

def shell_do(command, log=False, return_log=False):
    print(f'Executing: {(" ").join(command.split())}', file=sys.stderr)

    res=subprocess.run(command.split(), stdout=subprocess.PIPE)

    if log:
        print(res.stdout.decode('utf-8'))
    if return_log:
        return(res.stdout.decode('utf-8'))


# create gsmr input files
gsmr_exposure=f'{uid}.gsmr_exposure.txt'
gsmr_outcome=f'{uid}.gsmr_outcome.txt'

with open(gsmr_exposure, 'w')as f:
    f.write(f'{uid} {exposure_file}')
with open(gsmr_outcome, 'w')as f:
    f.write(f'PD {outcome_file}')

# conduct gsmr
shell_do(f'gcta64 --bfile {bfile} --gsmr-file  {gsmr_exposure} {gsmr_outcome} --gsmr-direction 0 --out {uid}.gsmr.res')