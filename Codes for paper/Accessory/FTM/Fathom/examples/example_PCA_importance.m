% Example of determining variable importance in PCA
% by David L. Jones, Mar-2016


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: '.../examples/legendre_gallagher.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This is an artificial data set comprising 19 sample sites and 9 species:
% Y       = abundance values of species for each site
% yLabels = cell array of species names
% sLabels = cell array of site names
% 
% (Source: Figure 3 in Legendre, P., and E. E. Gallagher. 2001.
% Ecologically meaningful transformations for ordination of species data.
% Oecologia 129: 271-280.)

% Load the data:
load legendre_gallagher.mat;

% PCA of raw data using covariance matrix:
pca = f_pca(Y,0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                PCA: Rank indicators via Scaled Eigenvectors:                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imp_01 = f_pca_importance(pca,1,1,yLabels);
% 
%     Frac_Explained    Descriptor
%     ______________    __________
%        0.29929        'sp2'     
%        0.29876        'sp4'     
%        0.26963        'sp3'     
%       0.060862        'sp5'     
%       0.060423        'sp1'     
%      0.0068585        'sp9'     
%      0.0024408        'sp8'     
%      0.0010316        'sp7'     
%     0.00070619        'sp6'     
% 
% (NOTE: Eigenvectors scaled by sqrt(eigenvalue))

% Show proportions account for 100%:
sum(imp_01.expl)
% ans = 1



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Alternative method:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imp_02 = f_pca_importance(pca,2,1,yLabels);
% 
%     Length    Descriptor
%     ______    __________
%     28.588    'sp2'     
%     28.563    'sp4'     
%     27.135    'sp3'     
%     12.892    'sp5'     
%     12.845    'sp1'     
%     4.3277    'sp9'     
%     2.5817    'sp8'     
%     1.6784    'sp7'     
%     1.3887    'sp6'     
% 
% (NOTE: Eigenvectors scaled by % variation explained)
