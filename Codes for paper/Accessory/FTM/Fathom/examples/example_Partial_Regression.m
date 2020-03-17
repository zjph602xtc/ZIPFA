% Example of performing a partial multipe regression.
% 
% -----Author:-----
% by David L. Jones, Aug-2016
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PARTIAL REGRESSION:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This example follows Section 10.3.5 Partial Regressin in Legendre & Legendre
% (1998)

% Import data:
raw = f_importCSV('ecothau.csv');

% Parse data:
Y = raw.dat(:,1);
X = raw.dat(:,2:4);
W = [raw.dat(:,5:6) raw.dat(:,5).^2];
clear raw;

% Regress Y on W, collect residuals:
temp = f_mregress(W,Y,0,0,0);
Yr   = temp.resid;

% Regress each X on W, collect residuals:
[nr,nc] = size(X);    % get size of X
Xr      = NaN(nr,nc); % preallocate
for i=1:nc
   temp = f_mregress(W,X(:,i),0,0,0);
   Xr(:,i) = temp.resid;
end

% Perform partial regression:
model = f_mregress(Xr,Yr,1000,1,1);
% =====================================================================
%  Multiple Linear Regression via QR Factorization:
% ---------------------------------------------------------------------
% R2            R2adj            F-stat        para-p           perm-p 
% ---------------------------------------------------------------------
% 0.31969       0.19213       2.50625       0.08330       0.10900       
% ---------------------------------------------------------------------
% 
% ---------------------------------------------------------------------
% Variable      b             t-stat        parametric-p  permutation-p
% ---------------------------------------------------------------------
%     intercept 0.00000       0.00000       1.00000       0.80800       
%             1 -0.89572      -2.26817      0.03517       0.06400       
%             2 -1.33929      -1.83036      0.08293       0.09400       
%             3 0.53663       0.52267       0.60724       0.64800       
% 
% -> the R^2 here is the 'partial R^2'
% 
% 
% Show equation for partial multiple regression:
% Yr_fit = 0 - 0.90Xr1 - 1.34Xr2 + 0.54Xr3; (R^2 = 0.32);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 PARTIAL RDA:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rda = f_rda(Y,X,W,1000,1,0);
% Permuting residuals under a reduced model 999 times...
% 
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 3.7456    p    =  0.0300 (iter=1000) 
% R2  = 0.3777   R2adj =  0.2611 
% R2w = 0.1853 (covariate) 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------

% Regress Y on W, collect residuals:
temp = f_rda(Y,W,[],0,0,0);
Yr   = temp.res;

% Regress X on W, collect residuals:
temp = f_rda(X,W,[],0,0,0);
Xr   = temp.res;

% Regress Yr on Xr;
rda = f_rda(Yr,Xr,[],1000,1,0);

% Permuting residuals under a full model 999 times...
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 2.8964    p    =  0.0620 (iter=1000) 
% R2  = 0.3519   R2adj =  0.2304 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% --------------------------------------------------

