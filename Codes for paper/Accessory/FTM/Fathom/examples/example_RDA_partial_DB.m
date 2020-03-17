% This example demonstrates that traditional RDA and Distance-based RDA produce
% equivalent results when based on the same underlying distance/dissimilarity
% metric.
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
% This example follows the partial RDA analysis of coral reef fish abundances
% from Legendre & Legendre (2012:p.653), but employs the Hellinger dissimilarity
% metric instead of the Euclidean distance.

% Clear workspace:
clz

% Load data file:
load legendre_113.mat

% Create categorical variable 'substrate type' (1=coral, 2=sand, 3=other)
sub    = [2 2 2 3 1 3 1 3 1 3]';
x_txt  = {'depth' 'coral' 'sand' 'other'}';
 
% Dummy code the substrate type (trim last column to avoid a singular matrix):
subRx = f_dummy(sub,1);

% Hellinger transform the abundance data:
y_H = f_hellinger(y);

% Determine the partial contribution of SUBSTRATE while controlling for the
% effect of DEPTH using traditional RDA on Hellinger-transformed data:
rda = f_rda(y_H,subRx,depth,1000,1);
% Permuting residuals under a reduced model 999 times...
% 
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 5.1962    p    =  0.00100 (iter=1000) 
% R2  = 0.5545   R2adj =  0.42725 
% R2w = 0.1253 (covariate) 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.1772  0.0499
% Residual Eigenvalues:
%   0.0898  0.0371  0.0041  0.0001  0.0000
% 
% Species-Environment Correlations (r):
%   0.9958  0.7278
% 
% Fraction of RESIDUAL variance explained (covariate removed):
% ------------------------------
% Canonical axes (total = 0.6340): 
%   0.4948  0.1392
% Cumulative: 
%   0.4948  0.6340
% 
% Residual axes  (total = 0.3660):
%   0.2508  0.1035  0.0114  0.0002  0.0000
% Cumulative: 
%   0.2508  0.3544  0.3658  0.3660  0.3660
% ------------------------------
% 
% (X & W variables standardized)
% ==================================================
% 
% 
% -> Since SUBSTRATE is correlated with DEPTH, we can see that SUBSTRATE
%    explains 55.5% of the variabiliity in fish abundance once the effect of
%    DEPTH has been removed.
% 
% -> 63% of the variablilily that remains after the covariate DEPTH has been
%   accounted for is explained by SUBSTRATE.
% 
% -> 0.6340 * (1 - 0.1253) = 0.5545



% Perform the same analysis using Distance-based RDA:
rdaDB = f_rdaDB(f_dis(y,'hel'),size(y,2),subRx,depth,1000,1);
% 
% Permuting residuals under a reduced model 999 times...
% 
% ==================================================
% Partial db-REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 5.1962    p    =  0.00300 (iter=1000) 
% R2  = 0.5545   R2adj =  0.42725 
% R2w = 0.1253 (covariate) 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   1.5948  0.4487
% Residual Eigenvalues:
%   0.8084  0.3338  0.0369  0.0007  0.0000  0.0000
% 
% Species-Environment Correlations (r):
%   0.9958  0.7278
% 
% Fraction of RESIDUAL variance explained (covariate removed):
% ------------------------------
% Canonical axes (total = 0.6340): 
%   0.4948  0.1392
% Cumulative: 
%   0.4948  0.6340
% 
% Residual axes  (total = 0.3660):
%   0.2508  0.1035  0.0114  0.0002  0.0000  0.0000
% Cumulative: 
%   0.2508  0.3544  0.3658  0.3660  0.3660  0.3660
% ------------------------------
% 
% (X & W variables standardized)
% ==================================================


% Show that the ordination scores from both methods offer the same depiction:
f_procrustes(rda.siteScores,rdaDB.siteScores,1,1000,1)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Forest Birds:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This example follows 'Case Study 1' in Leps & Smilauer (2003), but uses
% Hellinger dissimilarity instead of Euclidean distance.

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

% Hellinger transform the biotic data:
bio_2_H = f_hellinger(bio_2);

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

% Peform partial RDA:
rda = f_rda(bio_2_H,env(:,tar),env(:,1),1000,1);
% 
% Permuting residuals under a reduced model 999 times...
% 
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 1.6733    p    =  0.00100 (iter=1000) 
% R2  = 0.3176   R2adj =  0.10431 
% R2w = 0.0941 (covariate) 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.0385  0.0191  0.0142  0.0107  0.0081  0.0069  0.0055  0.0043  0.0035  0.0021
% Residual Eigenvalues:
%   0.0314  0.0220  0.0196  0.0171  0.0163  0.0135  0.0119  0.0102  0.0084  0.0079  0.0072  0.0057  0.0050  0.0045  0.0042  0.0041  0.0035  0.0031  0.0027  0.0022  0.0019  0.0016  0.0014  0.0012  0.0008  0.0007  0.0005  0.0004  0.0002  0.0001  0.0000
% 
% Species-Environment Correlations (r):
%   0.8669  0.7967  0.8082  0.7469  0.7324  0.6920  0.5892  0.5972  0.6568  0.4844
% 
% Fraction of RESIDUAL variance explained (covariate removed):
% ------------------------------
% Canonical axes (total = 0.3505): 
%   0.1196  0.0594  0.0440  0.0332  0.0252  0.0214  0.0171  0.0133  0.0109  0.0065
% Cumulative: 
%   0.1196  0.1790  0.2230  0.2562  0.2814  0.3028  0.3199  0.3332  0.3440  0.3505
% 
% Residual axes  (total = 0.6495):
%   0.0975  0.0682  0.0608  0.0531  0.0505  0.0419  0.0369  0.0315  0.0260  0.0244  0.0224  0.0178  0.0155  0.0139  0.0131  0.0127  0.0109  0.0096  0.0084  0.0070  0.0060  0.0049  0.0043  0.0037  0.0026  0.0021  0.0016  0.0012  0.0005  0.0004  0.0001
% Cumulative: 
%   0.0975  0.1657  0.2265  0.2797  0.3301  0.3721  0.4090  0.4405  0.4665  0.4910  0.5133  0.5312  0.5467  0.5606  0.5737  0.5864  0.5972  0.6068  0.6152  0.6222  0.6281  0.6330  0.6373  0.6410  0.6436  0.6457  0.6473  0.6485  0.6490  0.6493  0.6495
% ------------------------------
% 
% (X & W variables standardized)
% ==================================================
% 
% 
% -> these 10 environmental variables account for 32% of the total variability
% in the bird assemblage data; that is 32% in addition to the 9% already
% explained by the covariate ALTITUDE.
% 
% -> 35% of the variablilily that remains after the covariate ALTITUDE has been
% accounted for is explained by these 10 environmental variables.
% 
% -> i.e., 0.3505 * (1 - 0.0941) = 0.31752



% Perform the same analysis using Distance-based RDA:
rdaDB = f_rdaDB(f_dis(bio_2,'hel'),size(bio_2,2),env(:,tar),env(:,1),1000,1);
% 
% Permuting residuals under a reduced model 999 times...
% 
% ==================================================
% Partial db-REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 1.6733    p    =  0.00100 (iter=1000) 
% R2  = 0.3176   R2adj =  0.10431 
% R2w = 0.0941 (covariate) 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   1.6170  0.8037  0.5953  0.4486  0.3410  0.2891  0.2315  0.1793  0.1471  0.0880
% Residual Eigenvalues:
%   1.3182  0.9226  0.8227  0.7184  0.6828  0.5669  0.4990  0.4266  0.3520  0.3306  0.3024  0.2414  0.2103  0.1874  0.1770  0.1720  0.1468  0.1296  0.1131  0.0942  0.0806  0.0657  0.0584  0.0496  0.0351  0.0290  0.0212  0.0168  0.0064  0.0050  0.0014
% 
% Species-Environment Correlations (r):
%   0.8669  0.7967  0.8082  0.7469  0.7324  0.6920  0.5892  0.5972  0.6568  0.4844
% 
% Fraction of RESIDUAL variance explained (covariate removed):
% ------------------------------
% Canonical axes (total = 0.3505): 
%   0.1196  0.0594  0.0440  0.0332  0.0252  0.0214  0.0171  0.0133  0.0109  0.0065
% Cumulative: 
%   0.1196  0.1790  0.2230  0.2562  0.2814  0.3028  0.3199  0.3332  0.3440  0.3505
% 
% Residual axes  (total = 0.6495):
%   0.0975  0.0682  0.0608  0.0531  0.0505  0.0419  0.0369  0.0315  0.0260  0.0244  0.0224  0.0178  0.0155  0.0139  0.0131  0.0127  0.0109  0.0096  0.0084  0.0070  0.0060  0.0049  0.0043  0.0037  0.0026  0.0021  0.0016  0.0012  0.0005  0.0004  0.0001
% Cumulative: 
%   0.0975  0.1657  0.2265  0.2797  0.3301  0.3721  0.4090  0.4405  0.4665  0.4910  0.5133  0.5312  0.5467  0.5606  0.5737  0.5864  0.5972  0.6068  0.6152  0.6222  0.6281  0.6330  0.6373  0.6410  0.6436  0.6457  0.6473  0.6485  0.6490  0.6493  0.6495
% ------------------------------
% 
% (X & W variables standardized)
% ==================================================

% Show that the ordination scores from both methods offer the same depiction:
f_procrustes(rda.siteScores,rdaDB.siteScores,1,1000,1)
