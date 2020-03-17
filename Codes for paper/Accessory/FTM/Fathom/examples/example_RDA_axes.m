% Example of testing the significance of RDA axes
% 
% by David L. Jones, May-2017
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Legendre et al. 2011 - Dune:                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example follows the example in 'Test runs.txt' provided in the supplement
% to Legendre et al. (2011).
% 
% Legendre, P., J. Oksanen, and C. J. F. ter Braak. Testing the significance of
%  canonical axes in redundancy analysis. Methods in Ecology & Evolution 2:
%  269-277. (See 'test.axes.cov.R')

% Clear workspace:
clz;

% Import data:
Y = f_importCSV('dune_Y.csv',2);
X = f_importCSV('dune_X.csv',2);


% Test RDA Axes:
result = f_rdaAxes(Y.dat,X.dat,0,1000,1,1,1);
% 
% ==================================================
%           Test Significance of RDA Axes:          
% --------------------------------------------------
% 
%     Axis    Eigenvalue       F         p  
%     ____    __________    _______    _____
%     1       19.403         6.2284    0.001
%     2       12.703         4.0777    0.006
%     3       3.3475         1.0746    0.723
%     4       1.9406        0.62293    0.858
% 
% # iterations = 1000
% --------------------------------------------------
% 
% -> only the first 2 canonical axes are significant
% 
% 
% Results from Legendre et al. (2011):
% > res.cov.2.2$test.axes
%        Eigenvalue F.marginal F.forward P-marginal P-forward
% Axis.1  19.403173  6.2283811 4.4969915      0.001     0.002
% Axis.2  12.703136  4.0776821 3.6631420      0.001     0.001
% Axis.3   3.347531  1.0745511 1.0317059      0.736     0.731
% Axis.4   1.940593  0.6229265 0.6229265      0.866     0.866




% Test RDA Axes with Covariables:
idxX = logical([1 0 1 1]);
idxW = logical([0 1 0 0]);
f_rdaAxes(Y.dat,X.dat(:,idxX),X.dat(:,idxW),1000,1,1);
%
% ==================================================
%           Test Significance of RDA Axes:          
% --------------------------------------------------
% 
%     Axis    Eigenvalue       F         p  
%     ____    __________    _______    _____
%     1       12.734         4.0877    0.006
%     2       3.4959         1.1222    0.734
%     3       2.1203        0.68062    0.843
% 
% # iterations = 1000
% --------------------------------------------------
% 
% -> only the first canonical axis is significant
% 
% 
% Results from Legendre et al. (2011):
% >res.cov.2.3$test.axes
% # Eigenvalue F.marginal F.forward P-marginal P-forward
% # Axis.1  12.734367  4.0877074 3.6491343      0.002     0.002
% # Axis.2   3.495853  1.1221622 1.0734545      0.666     0.664
% # Axis.3   2.120327  0.6806207 0.6806207      0.811     0.811













