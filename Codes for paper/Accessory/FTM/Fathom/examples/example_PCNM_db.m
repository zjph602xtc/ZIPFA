% This example illustrates the equivalence of using either traditional RDA or
% Distance-based RDA in PCNM-based analysis in spatial ecology.
% 
% by David L. Jones, Feb-2016
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Notes:-----
% This example follows the analysis of Oribatid mites presented in Borcard
% et al. (2004).

% -----Variables:-----
% bio = structure of abundances of 35 species of mites from 70 cores:
%  .num = # per 5 cm diameter, 7 cm deep core
%  .txt = species codes
% 
% env = structure of environmental variables with the following fields:
%  .x       = x coordinates (meters)
%  .y       = y coordinates (meters)
%  .density = bulk density  (grams/liter of dry matter)
%  .water   = water content (grams/liter of fresh matter)
%  .substra = substrate (1 = Sphagn_1, 2 = Sphagn_2, 3 = Sphagn_3, 4 = Sphagn_4,
%             5 = Lignlitt, 6 = Barepeat, 7 = Interface)   
%  .schrub  = schrub coverage (density classes ranging from 1-3)
%  .micro   = microtopography (1 = blanket, 2 = hummock)

% Load data:
load oribatid_mites.mat

% Dummy code qualitative variables:
env.Qsubstra = f_dummy(env.substra,1);
env.Qschrub  = f_dummy(env.schrub,1);
env.Qmicro   = f_dummy(env.micro,1);

% Create PCNM's:
pcnm     = f_pcnm([env.x env.y]);
pcnm.txt = cellstr(num2str([1:size(pcnm.evects,2)]')); % create text labels
% 
% Use PCNM's to model the spatial structure of the biotic data; since the
% preliminary RDA indicated a significant linear gradient in the data, use
% linearly de-trended data for subsequent analyses (i.e., the residuals from the
% preliminary analysis to perform a partial RDA).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Traditional RDA:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Hellinger transform the species data:
bio.H = f_hellinger(bio.num);

% Perform preliminary RDA on X,Y coordinates to identify linear trends:
rda_XY = f_rda(bio.H,[env.x env.y],0,1000,1);
% 
% Permuting residuals under a full model 999 times...
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 13.2798    p    =  0.0010 (iter=1000) 
% R2  = 0.2839   R2adj =  0.2625 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.1058  0.0062
% Residual Eigenvalues:
%   0.0936  0.0354  0.0299  0.0167  0.0149  0.0120  0.0102  0.0084  0.0077  0.0069  0.0056  0.0044  0.0042  0.0037  0.0034  0.0032  0.0032  0.0024  0.0023  0.0021  0.0018  0.0016  0.0013  0.0012  0.0012  0.0011  0.0008  0.0007  0.0006  0.0006  0.0005  0.0003  0.0002  0.0002  0.0001
% 
% Species-Environment Correlations (r):
%   0.8308  0.5157
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.2839): 
%   0.2683  0.0156
% Cumulative: 
%   0.2683  0.2839
% 
% Residual axes  (total = 0.7161):
%   0.2373  0.0898  0.0758  0.0423  0.0377  0.0305  0.0258  0.0214  0.0196  0.0174  0.0141  0.0112  0.0105  0.0095  0.0085  0.0082  0.0080  0.0061  0.0059  0.0053  0.0046  0.0039  0.0033  0.0031  0.0030  0.0028  0.0020  0.0019  0.0016  0.0014  0.0013  0.0008  0.0006  0.0005  0.0002
% Cumulative: 
%   0.2373  0.3271  0.4030  0.4452  0.4830  0.5135  0.5393  0.5606  0.5802  0.5976  0.6117  0.6229  0.6335  0.6430  0.6515  0.6597  0.6677  0.6738  0.6797  0.6850  0.6897  0.6936  0.6969  0.7000  0.7030  0.7058  0.7078  0.7097  0.7113  0.7128  0.7140  0.7149  0.7154  0.7159  0.7161
% ------------------------------
% 
% (X variables standardized)
% ==================================================
% 
% -> the X,Y coordinates explain a significant portion (28.39%) of the
%    variability in the biotic data.


% Residuals serve as the de-trended data:
bio.dt = rda_XY.res;

% RDA using the de-trended data and the PCNM's:
rda_PCNM = f_rda(bio.dt,pcnm.evects,0,1000,1);
% 
% Permuting residuals under a full model 999 times...
% 
% ==================================================
% REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 2.0202    p    =  0.0010 (iter=1000) 
% R2  = 0.4860   R2adj =  0.2454 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   0.0612  0.0256  0.0127  0.0069  0.0062  0.0047  0.0037  0.0028  0.0025  0.0020  0.0018  0.0015  0.0013  0.0013  0.0009  0.0007  0.0005  0.0004  0.0003  0.0002  0.0000  0.0000
% Residual Eigenvalues:
%   0.0382  0.0209  0.0150  0.0098  0.0081  0.0072  0.0067  0.0051  0.0047  0.0045  0.0032  0.0029  0.0025  0.0022  0.0020  0.0018  0.0015  0.0014  0.0012  0.0009  0.0009  0.0007  0.0007  0.0005  0.0005  0.0004  0.0004  0.0003  0.0002  0.0002  0.0002  0.0001  0.0001  0.0001  0.0000
% 
% Species-Environment Correlations (r):
%   0.8212  0.8635  0.7202  0.7290  0.6197  0.6819  0.6314  0.6760  0.6344  0.6410  0.6090  0.5259  0.5838  0.5314  0.4740  0.4555  0.4309  0.3562  0.2954  0.3624  0.0998  0.0958
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.4860): 
%   0.2166  0.0908  0.0449  0.0245  0.0221  0.0167  0.0130  0.0101  0.0089  0.0071  0.0065  0.0052  0.0046  0.0044  0.0032  0.0026  0.0017  0.0014  0.0009  0.0005  0.0001  0.0001
% Cumulative: 
%   0.2166  0.3074  0.3523  0.3768  0.3990  0.4157  0.4287  0.4387  0.4477  0.4547  0.4612  0.4664  0.4710  0.4754  0.4786  0.4812  0.4830  0.4843  0.4853  0.4858  0.4860  0.4860
% 
% Residual axes  (total = 0.5140):
%   0.1354  0.0740  0.0532  0.0345  0.0286  0.0255  0.0237  0.0180  0.0167  0.0158  0.0114  0.0103  0.0088  0.0078  0.0071  0.0063  0.0052  0.0050  0.0043  0.0033  0.0032  0.0026  0.0023  0.0019  0.0018  0.0014  0.0013  0.0011  0.0009  0.0007  0.0005  0.0004  0.0003  0.0002  0.0001
% Cumulative: 
%   0.1354  0.2094  0.2626  0.2971  0.3258  0.3513  0.3750  0.3930  0.4097  0.4255  0.4370  0.4473  0.4561  0.4639  0.4710  0.4774  0.4825  0.4875  0.4919  0.4951  0.4984  0.5010  0.5033  0.5053  0.5070  0.5084  0.5097  0.5108  0.5116  0.5123  0.5129  0.5133  0.5136  0.5139  0.5140
% ------------------------------
% 
% (X variables standardized)
% ==================================================
% 
% -> the PCNM's explain 48.6% of the variability in the detrended biotic data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Distance-based RDA:                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Perform preliminary RDA on X,Y coordinates to identify linear trends:
rdaDB_XY = f_rdaDB(f_dis(bio.num,'hel'),size(bio.num,2),[env.x env.y],0,1000,1);
% 
% Permuting residuals under a full model 999 times...
% 
% ==================================================
% db-REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 13.2798    p    =  0.0010 (iter=1000) 
% R2  = 0.2839   R2adj =  0.2625 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   7.2980  0.4249
% Residual Eigenvalues:
%   6.4562  2.4429  2.0632  1.1506  1.0260  0.8298  0.7027  0.5809  0.5329  0.4728  0.3842  0.3050  0.2866  0.2579  0.2319  0.2235  0.2181  0.1655  0.1605  0.1454  0.1264  0.1074  0.0900  0.0840  0.0819  0.0753  0.0551  0.0508  0.0443  0.0390  0.0343  0.0229  0.0157  0.0135  0.0050
% 
% Species-Environment Correlations (r):
%   0.8308  0.5157
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.2839): 
%   0.2683  0.0156
% Cumulative: 
%   0.2683  0.2839
% 
% Residual axes  (total = 0.7161):
%   0.2373  0.0898  0.0758  0.0423  0.0377  0.0305  0.0258  0.0214  0.0196  0.0174  0.0141  0.0112  0.0105  0.0095  0.0085  0.0082  0.0080  0.0061  0.0059  0.0053  0.0046  0.0039  0.0033  0.0031  0.0030  0.0028  0.0020  0.0019  0.0016  0.0014  0.0013  0.0008  0.0006  0.0005  0.0002
% Cumulative: 
%   0.2373  0.3271  0.4030  0.4452  0.4830  0.5135  0.5393  0.5606  0.5802  0.5976  0.6117  0.6229  0.6335  0.6430  0.6515  0.6597  0.6677  0.6738  0.6797  0.6850  0.6897  0.6936  0.6969  0.7000  0.7030  0.7058  0.7078  0.7097  0.7113  0.7128  0.7140  0.7149  0.7154  0.7159  0.7161
% ------------------------------
% 
% (X variables standardized)
% ==================================================
% 
% -> the X,Y coordinates significantly explained 28.4% (or 26.3% if you use
% R^2adj) of the variance in the species data.


% Residual matrix serves as the de-trended data:
bio.Gres = rdaDB_XY.Gres;

% db-RDA using the de-trended data and the PCNM's:
skip       = 2; % specify input is a symmetric residual matrix
rdaDB_PCNM = f_rdaDB(bio.Gres,size(bio.num,2),pcnm.evects,0,1000,1,1,skip);

% ==================================================
% db-REDUNDANCY ANALYSIS:
% --------------------------------------------------
%  F  = 2.0202    p    =  0.0010 (iter=1000) 
% R2  = 0.4860   R2adj =  0.2454 
% --------------------------------------------------
% 
% Canonical Eigenvalues:
%   4.2207  1.7685  0.8745  0.4775  0.4312  0.3263  0.2528  0.1960  0.1738  0.1380  0.1265  0.1008  0.0893  0.0866  0.0621  0.0507  0.0338  0.0271  0.0183  0.0104  0.0028  0.0012
% Residual Eigenvalues:
%   2.6375  1.4414  1.0363  0.6730  0.5581  0.4976  0.4620  0.3498  0.3261  0.3084  0.2224  0.2011  0.1723  0.1520  0.1386  0.1231  0.1009  0.0973  0.0843  0.0641  0.0633  0.0510  0.0452  0.0378  0.0347  0.0266  0.0248  0.0211  0.0171  0.0136  0.0107  0.0086  0.0061  0.0040  0.0023
% 
% Species-Environment Correlations (r):
%   0.8212  0.8635  0.7202  0.7290  0.6197  0.6819  0.6314  0.6760  0.6344  0.6410  0.6090  0.5259  0.5838  0.5314  0.4740  0.4555  0.4309  0.3562  0.2954  0.3624  0.0998  0.0958
% 
% Fraction of variance explained:
% ------------------------------
% Canonical axes (total = 0.4860): 
%   0.2166  0.0908  0.0449  0.0245  0.0221  0.0167  0.0130  0.0101  0.0089  0.0071  0.0065  0.0052  0.0046  0.0044  0.0032  0.0026  0.0017  0.0014  0.0009  0.0005  0.0001  0.0001
% Cumulative: 
%   0.2166  0.3074  0.3523  0.3768  0.3990  0.4157  0.4287  0.4387  0.4477  0.4547  0.4612  0.4664  0.4710  0.4754  0.4786  0.4812  0.4830  0.4843  0.4853  0.4858  0.4860  0.4860
% 
% Residual axes  (total = 0.5140):
%   0.1354  0.0740  0.0532  0.0345  0.0286  0.0255  0.0237  0.0180  0.0167  0.0158  0.0114  0.0103  0.0088  0.0078  0.0071  0.0063  0.0052  0.0050  0.0043  0.0033  0.0032  0.0026  0.0023  0.0019  0.0018  0.0014  0.0013  0.0011  0.0009  0.0007  0.0005  0.0004  0.0003  0.0002  0.0001
% Cumulative: 
%   0.1354  0.2094  0.2626  0.2971  0.3258  0.3513  0.3750  0.3930  0.4097  0.4255  0.4370  0.4473  0.4561  0.4639  0.4710  0.4774  0.4825  0.4875  0.4919  0.4951  0.4984  0.5010  0.5033  0.5053  0.5070  0.5084  0.5097  0.5108  0.5116  0.5123  0.5129  0.5133  0.5136  0.5139  0.5140
% ------------------------------
% 
% (X variables standardized)
% ==================================================
% 
% -> the PCNM's explain 48.6% of the variability in the detrended biotic data



rda_PCNM = f_rda(bio.dt,pcnm.evects(:,[1 3 4 6 7 10 11 20]),0,1000,1);