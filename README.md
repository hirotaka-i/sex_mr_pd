# Sex specific MR analysis for PD

 - **Project:** sex specific biomarker MR analysis for PD
 - **Author(s):** Hirotaka Iwaki 
 - **Date Last Updated:** Oct 29 2021
  - **Status:** Incomplete
    - **Update Description:** Starting analysis



## Overview
In this analysis, we conducted sex-stratified MR analyses for biomarkers and PD


## Exposure GWAS
The following potentially sex-associated biomarkers were selected and the sex-specific summary-stats were downloaded from [Neale Lab GWAS of UK Biobank phenotypes](https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?usp=sharing).

### Phenotypes
    pheno_dic = {
        '30880_irnt':'Urate (quantile)', 
        '30880_raw':'Urate (umol/L)', 
        '30800_irnt':'Oestradiol (quantile)', 
        '30800_raw':'Oestradiol (pmol/L)', 
        '30850_irnt':'Testosterone (quantile)', 
        '30850_raw':'Testosterone (nmol/L)', 
        '30770_irnt':'IGF-1 (quantile)', 
        '20003_1140883066':'Treatment/medication code: insulin product',
        '100240':'Coffee consumed', 
        '20160':'Ever smoked', 
        '2644':'Light smokers, at least 100 smokes in lifetime',
        '30840_irnt':'Total bilirubin (quantile)', 
        '30830_irnt':'SHBG (quantile)',
        "30680_irnt":"Calcium (quantile)"
    }

## Outcome GWAS
The sex specific PD GWAS summary stats were downloaded from [IPDGC resource for Male and Female specific summary stats](https://pdgenetics.org/resources).    
\* [The previous research](https://onlinelibrary.wiley.com/doi/pdf/10.1002/ana.26090) showed that there are no autosomal riks structure differences in PD. 

## Overview
* LD score regression between men and women for the phenotype    
  Conduct [LDSC regression](https://github.com/bulik/ldsc/wiki) between male-biomarker-GWAS and female-biomarker_GWAS to check if the genetic architecture were different beetween men and women for that biomarker

* Used MR models: Software
  * IVW, Egger: [TwoSampleMR](https://mrcieu.github.io/TwoSampleMR/articles/gwas2020.html)
  * GSMR: [GCTA-GSMR](https://cnsgenomics.com/software/gcta/#GSMR)
  * MR-LDP: [MR-LDP](https://github.com/QingCheng0218/MR.LDP)


## Analysis script
Please refer to the Main.ipynb