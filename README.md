# Zero inflated poison factor analysis

<br>

**The method has a Matlab version and R version. The Matlab version codes are in this repository. The R version is on CRAN [(click)](https://cran.rstudio.com/web/packages/ZIPFA/index.html)**

**Here is the tutorial on both versions. [(click)](https://zjph602xtc.github.io/ZIPFA/)**

**This repository contains these files:**

Simulation folder:\
'uvax.mat': simulated ln lambda in Web Figure 1;\
'examplefitted.csv': Example data for Web Figure 2.\


Real Data folder:\
'alldataX.mat': ORIGINS data. (not included. Please email us if you need it.)\
'originsdataMar2017.csv': Original ORIGINS data. (not included. Please email us if you need it.)\
'finaluv.mat': ZIPFA results on ORIGINS data.\
'final_u.csv'/'final_v.csv': ZIPFA results on ORIGINS data.\
'lsvd_u.csv'/'lsvd_v.csv': log-PCA results on ORIGINS data.\
'gomms_u.csv'/'gomms_v.csv': gomms results on ORIGINS data.\
'psvdos_u.csv'/'psvdos_v.csv': psvdos results on ORIGINS data.\
'pcoa_u.csv': PCoA results on ORIGINS data.\
'nmds_u.csv': nMDS results on ORIGINS data.\
'taxa.csv': Taxa name list.\
'our_result.xlsx': Data for figure 5.\



'ZIPFA.m': Zero inflated Poisson factor analysis. \
'cv_ZIPFA.m': Cross validation on ZIPFA.\
'EMzeropoisson_mat.m': Zero inflated Poisson regression.\
'Examples codes.m': Example runs for 'EMzeropoisson_mat','ZIPFA','cv_ZIPFA'.\

'Table and Figure in paper.m': All tables and figures in the paper. \

'gomms.R': GOMMS model.\
'Simu.R': Some tables and figures in the paper. \

Accessory folder:\
'cluster_fit.m'/'heatmap.m'/'pnew.m'/'hatchfill2_r8'/'FTM': Accessory codes. \
   

