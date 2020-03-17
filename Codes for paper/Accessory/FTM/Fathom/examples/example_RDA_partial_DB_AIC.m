% Example of using AIC-based forward addition to select an optimal subset of
% explanatory variables when covariates are present - Partial db-RDA

% by David L. Jones, Apr-2016
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Forest Birds:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example uses data from 'Case Study 1' in Leps & Smilauer (2003)

% Clear workspace:
clz

% Import & parse biotic data:
raw     = f_importCSV('birds_bio.csv',2);
bio     = raw.dat;
bio_txt = raw.txt;
site    = raw.site;

% Import & parse environmental data:
raw     = f_importCSV('birds_env.csv',2);
env     = raw.dat;
env_txt = raw.txt;
clear raw;

% Apply a square root transform to the biotic data:
bio_2 = f_normal(bio,'2');

% Note: We need to remove redundant explanatory variables to prevent a 'rank
% deficient' error when performing the regression analysis. The 2 variables
% Rocks/NoRocks are actually dummy codes for a single categorical variable
% having 2 levels, so we need to use only one or the other in the analysis. The
% same applies to Warm/Cold.
% 
tar = [2:10 12]; % create index to target environmental variables
env_txt(tar)'    % show target environmental variables

% ans = 
%     'Forest'
%     'ForDens'
%     'BrLeaf'
%     'E2'
%     'E2Con'
%     'E1'
%     'E1Height'
%     'Slope'
%     'Rocks'
%     'Warm'

% Perform partial RDA to determine whether there are any significant effects of
% the remaining environmental variables after removing the variability in bird
% assemblages explained by ALTITUDE.

% First select an optimal subset of predictors:
X     = env(:,tar);   % predictor  variables
X_txt = env_txt(tar); % predictor labels
W     = env(:,1);     % set ALTIT as a covariate

best = f_rdaDB_AIC(f_dis(bio_2,'bc'),size(bio_2,2),X,W,0,1,2,X_txt);

% ==================================================
% AIC-based forward addition: partial db-RDA 
% --------------------------------------------------
% Conditional Tests: (each variable separately)
%     'RSS'       'R2'          'R2adj'        'AIC'        'delta'     'wts'          'ratio'     'var'     
%     [3.1597]    [  0.0995]    [ 0.077537]    [-107.96]    [     0]    [  0.24629]    [     1]    'BrLeaf'  
%     [ 3.241]    [ 0.07948]    [ 0.057029]    [-106.87]    [1.0924]    [  0.14264]    [1.7267]    'Warm'    
%     [3.2559]    [0.075815]    [ 0.053274]    [-106.67]    [1.2894]    [  0.12926]    [1.9054]    'E2Con'   
%     [3.2725]    [0.071735]    [ 0.049094]    [-106.45]    [1.5077]    [  0.11589]    [2.1252]    'E1'      
%     [3.3065]    [0.063365]    [ 0.040521]    [-106.01]    [ 1.952]    [ 0.092806]    [2.6538]    'Slope'   
%     [3.3415]    [0.054741]    [ 0.031686]    [-105.56]    [2.4051]    [ 0.073993]    [3.3285]    'Forest'  
%     [3.3721]    [0.047211]    [ 0.023973]    [-105.16]    [2.7968]    [ 0.060833]    [4.0487]    'E1Height'
%     [3.3826]    [0.044613]    [ 0.021311]    [-105.03]    [2.9311]    [  0.05688]    [  4.33]    'E2'      
%     [3.4319]    [0.032488]    [0.0088898]    [-104.41]    [3.5526]    [ 0.041689]    [5.9078]    'Rocks'   
%     [3.4538]    [0.027096]    [0.0033662]    [-104.13]    [ 3.826]    [ 0.036361]    [6.7735]    'ForDens' 
%     [ 4.061]    [     NaN]    [      NaN]    [-99.373]    [8.5878]    [0.0033623]    [ 73.25]    'none'    
% --------------------------------------------------
% 
% Marginal Tests: (sequential variable addition)
%     'RSS'       'R2'        'R2adj'       'AIC'        'wts'        'deltaN'     'var'       'idx'
%     [3.1597]    [0.0995]    [0.077537]    [-107.96]    [0.24629]    [ 8.5878]    'BrLeaf'    [  3]
%     [3.1597]    [   NaN]    [     NaN]    [-107.96]    [0.11058]    [0.99111]    'none'         []
% --------------------------------------------------
% 
% RSS    = residual sum-of-squares 
% R2     = fraction of total variance explained 
% R2adj  = fraction of adjusted total variance explained 
% AIC    = corrected AIC 
% deltaN = delta associated with NO variable addition 
% wts    = AIC weights 
% var    = variable labels 
% idx    = index to selected variables 
% 
% (Note: RSS, R2, and R2adj in Marginal tests are CUMULATIVE) 


% Perform final partial db-RDA using optimal subset:
rdaDB = f_rdaDB(f_dis(bio_2,'bc'),size(bio_2,2),X(:,best.idx),W,1000,1);
% 
% Permuting residuals under a reduced model 999 times...
% 
% ==================================================
% Partial db-REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 5.1152    p    =  0.0010 (iter=1000) 
% R2  = 0.0995   R2adj =  0.0775 
% R2w = 0.1224 (covariate) 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.4041
% Residual Eigenvalues:
%   0.6017  0.5346  0.4350  0.3272  0.3119  0.2423  0.2058  0.1716  0.1633  0.1586  0.1338  0.1227  0.0991  0.0931  0.0605  0.0555  0.0445  0.0356  0.0291  0.0180  0.0102  0.0068  0.0015  -0.0053  -0.0079  -0.0125  -0.0166  -0.0208  -0.0256  -0.0264  -0.0309  -0.0441  -0.0499  -0.0530
% 
% Species-Environment Correlations (r):
%   0.6808
% 
% Fraction of RESIDUAL variance explained (covariate removed):
% ------------------------------
% Canonical axes (total = 0.1017): 
%   0.1017
% Cumulative: 
%   0.1017
% 
% Residual axes  (total = 0.8983):
%   0.1514  0.1345  0.1095  0.0823  0.0785  0.0610  0.0518  0.0432  0.0411  0.0399  0.0337  0.0309  0.0249  0.0234  0.0152  0.0140  0.0112  0.0089  0.0073  0.0045  0.0026  0.0017  0.0004  -0.0013  -0.0020  -0.0032  -0.0042  -0.0052  -0.0064  -0.0066  -0.0078  -0.0111  -0.0125  -0.0133
% Cumulative: 
%   0.1514  0.2860  0.3954  0.4778  0.5563  0.6173  0.6691  0.7123  0.7534  0.7933  0.8270  0.8578  0.8828  0.9062  0.9214  0.9354  0.9466  0.9555  0.9629  0.9674  0.9699  0.9716  0.9720  0.9707  0.9687  0.9656  0.9614  0.9562  0.9497  0.9431  0.9353  0.9242  0.9116  0.8983
% ------------------------------
% 
% (X & W variables standardized)
% ==================================================


% -> 10% of the variablilily that remains after the covariate ALTITUDE has been
% accounted for is explained by BrLeaf
