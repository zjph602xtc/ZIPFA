% Examples of using partial RDA when covariates are present.
% 
% by David L. Jones, Feb-2016
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               CORAL REEF FISH:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/legendre_113.mat'
% 
% This example reproduces the partial RDA analysis of coral reef fish abundances
% from Legendre & Legendre (2012:p.653).

% Clear workspace:
clz

% Load data file:
load legendre_113.mat

% Create categorical variable 'substrate type' (1=coral, 2=sand, 3=other)
sub    = [2 2 2 3 1 3 1 3 1 3]';
x_txt  = {'depth' 'coral' 'sand' 'other'}';
 
% Dummy code the substrate type (trim last column to avoid a singular matrix):
subRx = f_dummy(sub,1);

% Determine the partial contribution of SUBSTRATE while controlling for the
% effect of DEPTH:
rda = f_rda(y,subRx,depth,1000,1);
% 
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 54.5597    p    =  0.00100 (iter=1000) 
% R2  = 0.7327   R2adj =  0.65635 
% R2w = 0.2270 (covariate) 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   64.1283  18.5868
% Residual Eigenvalues:
%   4.1888  0.3139  0.0370  0.0085
% 
% Species-Environment Correlations (r):
%   0.9428  0.9022
% 
% Fraction of RESIDUAL variance explained (covariate removed):
% ------------------------------
% Canonical axes (total = 0.9479): 
%   0.7349  0.2130
% Cumulative: 
%   0.7349  0.9479
% 
% Residual axes  (total = 0.0521):
%   0.0480  0.0036  0.0004  0.0001
% Cumulative: 
%   0.0480  0.0516  0.0520  0.0521
% ------------------------------
% 
% (X & W variables standardized)
% ==================================================
% 
% -> Since SUBSTRATE is correlated with DEPTH, we can see that SUBSTRATE
%    explains 73.3% of the variabiliity in fish abundance once the effect of
%    DEPTH has been removed (semipartial R^2 = 0.7327)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Forest Birds:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This example follows 'Case Study 1' in Leps & Smilauer (2003).

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

% Apply a natural log transform to the biotic data:
bio_ln = f_normal(bio,'ln');

% Perform RDA to determine how much of the variability in the bird assemblages
% can be explained by ALTITUDE:
rda_1 = f_rda(bio_ln,env(:,1),0,1000,1);
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 5.3233    p    =  0.00100 (iter=1000) 
% R2  = 0.1149   R2adj =  0.09333 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.1990
% Residual Eigenvalues:
%   0.2633  0.2046  0.1295  0.1163  0.1064  0.0880  0.0768  0.0680  0.0631  0.0554  0.0526  0.0484  0.0401  0.0319  0.0293  0.0213  0.0201  0.0176  0.0158  0.0150  0.0131  0.0112  0.0079  0.0074  0.0056  0.0049  0.0041  0.0037  0.0030  0.0025  0.0021  0.0015  0.0012  0.0006  0.0004  0.0003  0.0001
% 
% Species-Environment Correlations (r):
%   0.7918
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.1149): 
%   0.1149
% Cumulative: 
%   0.1149
% 
% Residual axes  (total = 0.8851):
%   0.1520  0.1181  0.0748  0.0672  0.0614  0.0508  0.0444  0.0392  0.0364  0.0320  0.0304  0.0280  0.0232  0.0184  0.0169  0.0123  0.0116  0.0101  0.0091  0.0087  0.0076  0.0064  0.0046  0.0043  0.0032  0.0028  0.0024  0.0021  0.0017  0.0015  0.0012  0.0008  0.0007  0.0003  0.0002  0.0002  0.0001
% Cumulative: 
%   0.1520  0.2701  0.3449  0.4121  0.4735  0.5243  0.5687  0.6079  0.6443  0.6763  0.7067  0.7347  0.7578  0.7762  0.7931  0.8054  0.8170  0.8271  0.8362  0.8449  0.8525  0.8589  0.8635  0.8678  0.8710  0.8738  0.8762  0.8783  0.8801  0.8815  0.8828  0.8836  0.8843  0.8846  0.8848  0.8850  0.8851
% ------------------------------
% 
% (X variables standardized)
% ==================================================
% 
% -> ALTITUDE explains a significant portion (11.5%) of the variability in bird
%   assemblages.



% -----Partial RDA:-----
% 
% Perform partial RDA to determine whether there are any significant effects of
% the remaining environmental variables after removing the variability in bird
% assemblages explained by ALTITUDE.

% Note: We need to remove redundant explanatory variables to prevent a 'rank
% deficient' error when performing the regression analysis. The 2 variables
% Rocks/NoRocks are actually dummy codes for a single categorical variable
% having 2 levels, so we need to use only one in the analysis. The same applies
% to Warm/Cold.
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

% Peform partial RDA:
rda_2 = f_rda(bio_ln,env(:,tar),env(:,1),1000,1);
% 
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 1.8214    p    =  0.00100 (iter=1000) 
% R2  = 0.3276   R2adj =  0.11744 
% R2w = 0.1149 (covariate) 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.1938  0.1167  0.0721  0.0590  0.0363  0.0270  0.0230  0.0168  0.0122  0.0104
% Residual Eigenvalues:
%   0.1471  0.1190  0.0917  0.0887  0.0687  0.0619  0.0536  0.0460  0.0435  0.0381  0.0301  0.0264  0.0221  0.0212  0.0173  0.0145  0.0133  0.0122  0.0105  0.0091  0.0064  0.0054  0.0049  0.0038  0.0035  0.0024  0.0015  0.0013  0.0007  0.0005  0.0002
% 
% Species-Environment Correlations (r):
%   0.8297  0.7564  0.8129  0.6533  0.7077  0.7076  0.6291  0.4869  0.5253  0.4942
% 
% Fraction of RESIDUAL variance explained (covariate removed):
% ------------------------------
% Canonical axes (total = 0.3701): 
%   0.1264  0.0761  0.0471  0.0385  0.0237  0.0176  0.0150  0.0110  0.0079  0.0068
% Cumulative: 
%   0.1264  0.2025  0.2496  0.2881  0.3118  0.3294  0.3444  0.3553  0.3633  0.3701
% 
% Residual axes  (total = 0.6299):
%   0.0959  0.0776  0.0598  0.0578  0.0448  0.0403  0.0350  0.0300  0.0284  0.0248  0.0196  0.0172  0.0144  0.0139  0.0113  0.0095  0.0087  0.0080  0.0068  0.0059  0.0042  0.0035  0.0032  0.0025  0.0023  0.0016  0.0010  0.0008  0.0004  0.0003  0.0001
% Cumulative: 
%   0.0959  0.1736  0.2334  0.2912  0.3361  0.3764  0.4114  0.4414  0.4698  0.4946  0.5142  0.5314  0.5459  0.5597  0.5710  0.5804  0.5891  0.5971  0.6039  0.6099  0.6140  0.6176  0.6208  0.6233  0.6256  0.6272  0.6281  0.6290  0.6294  0.6298  0.6299
% ------------------------------
% 
% (X & W variables standardized)
% ==================================================
% 
% 
% -> these 10 environmental variables account for 32.8% of the total variability
% in the bird assemblage data; that is 32% in addition to the 11.5% already
% explained by the covariate ALTITUDE.
% 
% -> 37% of the variablilily that remains after the covariate ALTITUDE has been
% accounted for is explained by these 10 environmental variables.
% 
% -> i.e., 0.3701 * (1 - 0.1149) = 0.32758
