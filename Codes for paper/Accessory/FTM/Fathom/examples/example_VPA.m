% Example of performing Variation Partitioning Analysis (VPA)
% 
% -----Author:-----
% by David L. Jones, Sep-2016
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% This example follows Section 10.5 Variation Partitioning in Legendre &
% Legendre (2012: p.574)

% Clear workspace:
clz;

% Import data:
raw = f_importCSV('ecothau.csv');

% Parse data:
Y = raw.dat(:,1);
X = raw.dat(:,2:4);
W = [raw.dat(:,5:6) raw.dat(:,5).^2];
clear raw;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Traditional RDA:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variation Partitioning:
vpa = f_vpa(Y,X,W,[],1000,1,1);
% 
% ==================================================
%       Variation Partitioning Analysis (VPA):        
% --------------------------------------------------
% 
%     Fraction      R2adj        p  
%     ________    _________    _____
%     'A'           0.11827    0.142
%     'B'           0.26339      NaN
%     'C'         0.0096727    0.376
%     'D'           0.60867      NaN
% 
% --------------------------------------------------
% Fractions for each explanatory matrix:
%    X = 0.382  
%    W = 0.273  
% --------------------------------------------------

% Save PDF:
f_pdf('ecothau_VPA')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Distance-based RDA:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variation Partitioning:
vpaDB = f_vpaDB(f_dis(Y,'euc'),size(Y,2),X,W,[],1000,1,1);
% 
% ==================================================
%     Variation Partitioning Analysis (db-VPA):       
% --------------------------------------------------
% 
%     Fraction      R2adj        p  
%     ________    _________    _____
%     'A'           0.11827    0.158
%     'B'           0.26339      NaN
%     'C'         0.0096727    0.423
%     'D'           0.60867      NaN
% 
% --------------------------------------------------
% Fractions for each explanatory matrix:
%    X = 0.382  
%    W = 0.273  
% --------------------------------------------------
