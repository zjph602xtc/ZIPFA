% Example of Repeated Measures RDA
% 
% by David L. Jones, May-2017
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Repeated Measures Partial RDA - Ohraz:                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example follows Case Study 5 from Leps & Smilauer (2003) and demonstrates
% how to use partial RDA to examine a repeated measures BACI experiment.

% Clear workspace:
clz;

% Import and parse data:
raw   = f_importCSV('ohraz_Y.csv',1);
Y.dat = raw.dat(:,3:end);
Y.txt = raw.txt(3:end);
raw   = f_importCSV('ohraz_X.csv',1);
yr    = raw.dat(:,1);
plt   = raw.dat(:,2);
mow   = raw.dat(:,3);
fert  = raw.dat(:,4);
remov = raw.dat(:,5);
clear raw;

% Create dummy codes for plt:
plt_Rx = f_dummy(plt);

% Create interaction terms:
Yr_M = f_interaction(yr,mow,0);
Yr_F = f_interaction(yr,fert,0);
Yr_R = f_interaction(yr,remov,0);

% partial RDA with permutations restricted within groups:
verb = 1;
stnd = 1;
perm = 0;
simu = 0;
grp = plt;

C1 = f_rda(Y.dat,[yr Yr_M Yr_F Yr_R],plt_Rx,1000,verb,stnd,perm,simu,grp);
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 5.3777    p    =  0.0010 (iter=1000) 
% R2  = 0.1294   R2adj =  0.0912 
% R2w = 0.4614 (covariate) 
% --------------------------------------------------


C2 = f_rda(Y.dat,[Yr_M Yr_F Yr_R],[yr plt_Rx],1000,verb,stnd,perm,simu,grp);
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 2.7569    p    =  0.0010 (iter=1000) 
% R2  = 0.0498   R2adj =  0.0188 
% R2w = 0.5411 (covariate) 
% --------------------------------------------------


C3 = f_rda(Y.dat,Yr_F,[yr Yr_M Yr_R plt_Rx],1000,verb,stnd,perm,simu,grp);
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 4.3957    p    =  0.0020 (iter=1000) 
% R2  = 0.0264   R2adj =  0.0161 
% R2w = 0.5644 (covariate) 
% --------------------------------------------------


C4 = f_rda(Y.dat,Yr_M,[yr Yr_F Yr_R plt_Rx],1000,verb,stnd,perm,simu,grp);
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 2.5012    p    =  0.0200 (iter=1000) 
% R2  = 0.0150   R2adj =  0.0046 
% R2w = 0.5758 (covariate) 
% --------------------------------------------------


C5 = f_rda(Y.dat,Yr_R,[yr Yr_M Yr_F plt_Rx],1000,verb,stnd,perm,simu,grp);
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 1.3739    p    =  0.0910 (iter=1000) 
% R2  = 0.0083   R2adj =  -0.0023 
% R2w = 0.5826 (covariate) 
% --------------------------------------------------





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Repeated Measures Partial db-RDA - Ohraz:                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - repeat the analysis above, but use Distance-based instead of classical RDA

% partial db-RDA with permutations restricted within groups:
ncY  = size(Y.dat,2);
verb = 1;
stnd = 1;
skip = 0;
simu = 0;
grp = plt;

C1 = f_rdaDB(f_dis(Y.dat,'euc'),ncY,[yr Yr_M Yr_F Yr_R],plt_Rx,1000,verb,stnd,skip,simu,grp);
% 
% ==================================================
% Partial db-REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 5.3777    p    =  0.0010 (iter=1000) 
% R2  = 0.1294   R2adj =  0.0912 
% R2w = 0.4614 (covariate) 
% --------------------------------------------------





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Repeated Measures Partial RDA - Wine:                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This example follows Section 29.2 of Neter et al (1996) and shows how to
% perform a one-way repeated measures univariate ANOVA using partial RDA.

% Clear workspace:
clz;

% Import and parse data:
raw   = f_importCSV('wine.csv',1);
Y     = raw.dat(:,1);
judge = raw.dat(:,2);
wine  = raw.dat(:,3);
clear raw;

% Create dummy codes for plt:
J_Rx = f_dummy(judge);
W_Rx = f_dummy(wine);

% partial RDA with permutations restricted within groups:
verb = 1;
stnd = 1;
perm = 0;
simu = 0;
grp  = judge;

f_rda(Y,W_Rx,J_Rx,1000,verb,stnd,perm,simu,grp);
% 
% ==================================================
% Partial REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 57.5000    p    =  0.0010 (iter=1000) 
% R2  = 0.4416   R2adj =  0.3523 
% R2w = 0.5200 (covariate) 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% 
% (X & W variables standardized)
% 
% (Restricted permutations used)
% --------------------------------------------------
% 
% -> we reject the null hypothesis that the mean wine ratings are not
%    significantly different





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Repeated Measures Partial db-RDA - Wine:                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - repeat the analysis above, but use Distance-based instead of classical RDA

% partial db-RDA with permutations restricted within groups:
ncY  = size(Y,2);
verb = 1;
stnd = 1;
skip = 0;
simu = 0;
grp  = judge;

f_rdaDB(f_dis(Y,'euc'),ncY,W_Rx,J_Rx,1000,verb,stnd,skip,simu,grp);
% 
% ==================================================
% Partial db-REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 57.5000    p    =  0.0010 (iter=1000) 
% R2  = 0.4416   R2adj =  0.3523 
% R2w = 0.5200 (covariate) 
% --------------------------------------------------
% Response variable (Y) is univariate. 
% 
% (X & W variables standardized)
% 
% (Restricted permutations used)
% --------------------------------------------------

