% Example of assessing complete spatial randomness (CSR) using the Clark-Evans
% index
% 
% -----Author:-----
% by David L. Jones, Jan-2016
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Load California Redwood data:
Y = load('redwood.dat');
Y = -Y; % reverse the sign

% Plot the data:
figure;
plot(Y(:,1),Y(:,2),'bo');

% Test the null hypothesis that the redwood tree spatially distributed according
% to a random model.
% 
result = f_clarkEvans(Y,1000,1,1)
% 
% ==================================================
%       Clark-Evans - Complete Spatial Randomness:
% --------------------------------------------------
% CE index = 0.7839  p =  0.00100 
% Z stat   = -3.2548  p =  0.00057 
% No. of permutations = 1000 
% --------------------------------------------------
% 
% Based on the p-values for both the CE index and Z-stat, we reject the null
% hypothesis of complete spatial randomness for these data (at alpha = 0.05).
% Since CE < 1, we conclude these data exhibit a significant aggregated pattern.
