% Example of using the modified version of PERMANOVA described in Anderson et
% al. (2017). This method is used to test hypotheses concerning central tendency
% while accounting for unbalanced designs with heterogeneous multivariate
% dispersions.
% 
% -----Author:-----
% by David L. Jones, May-2017
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         NORWEGIAN CONTINENTAL SHELF:                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This example follows the ecological application presented in section 4 of
% Anderson et al. (2017).
% 
% Anderson, M. J., D. C. I. Walsh, K. R. Clarke, R. N. Gorley, and E.
% Guerra-Castro. 2017. Some solutions to the multivariate Behrens-Fisher problem
% for dissimilarity-based analyses. Aust. N. Z. J. Stat. 59(1): 57-79.

% Clear workspace:
clz;

% Import & parse data:
raw   = f_importCSV('Norway.csv',1);
Y     = raw.dat(:,3:end);
Y_txt = raw.txt(3:end);
grp   = raw.dat(:,1);
clear raw;

% Create Jaccard dissimilarity matrix:
yDis = f_dis(Y,'jac');

% Test the null hypothesis of no differences in centroids among groups:
f_permanova(yDis,grp,1000,1);
% 
% ==================================================
%                Modified PERMANOVA:                
% --------------------------------------------------
%  F  = 13.5099    p    =  0.0010 (iter=1000) 
% --------------------------------------------------
% 
% -> reject the null hypothesis at alpha=0.05


% Perform pair-wise tests to see which groups differed:
f_permanovaPW(yDis,grp,1000,1);
% ----------------------------------------------------------
% Results of pair-wise comparisons between each factor level:
% ==========================================================
%          t:     p:     p_bon: p_ds:  p_holm:
% 1 vs. 2: 3.6287 0.0010 0.0100 0.0100 0.0100 
% 1 vs. 3: 3.7141 0.0010 0.0100 0.0100 0.0100 
% 1 vs. 4: 5.0670 0.0010 0.0100 0.0100 0.0100 
% 1 vs. 5: 5.0187 0.0010 0.0100 0.0100 0.0100 
% 2 vs. 3: 2.6151 0.0010 0.0100 0.0100 0.0100 
% 2 vs. 4: 4.3762 0.0010 0.0100 0.0100 0.0100 
% 2 vs. 5: 4.2595 0.0010 0.0100 0.0100 0.0100 
% 3 vs. 4: 2.7733 0.0010 0.0100 0.0100 0.0100 
% 3 vs. 5: 2.9003 0.0010 0.0100 0.0100 0.0100 
% 4 vs. 5: 3.1949 0.0010 0.0100 0.0100 0.0100 
% ----------------------------------------------------------
% t      = t-statistic 
% p      = unadjusted p-value 
% p_bon  = Bonferroni adjusted p-value 
% p_ds   = Dunn-Sidak adjusted p-value 
% p_holm = Holms adjusted p-value 
% 
% -> all pairs of groups were significantly different in terms of the location
% of their centroids, regardless of whether their dispersions differed or not.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Tikus Coral Assemblages:                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The file 'tikus.mat' contains data on coral assemblages from Tikus
% Island, Indonesia taken along 10 replicate 30 m line transects during
% each of 6 years (Warwick et al., 1990). The 1982 El Niño triggered a
% disturbance event which resulted in widespread coral bleaching in this
% region. As such, the 1981, 1983, and 1985 samples were collected before,
% during, and after the El Niño disturbance event, repectively.
% 
% Warwick, R. M., K. R. Clarke, and Suharsono. 1990. A statistical analysis
% of coral community responses to the 1982-1983 El Niño in the Thousand
% Islands, Indonesia. Coral Reefs 8:171-179.
% 
% -----Variables:----
% coral     = percent cover of 75 species in each transect
% coral_txt = cell array of corresponding column labels
% yr        = year (1 = 1981, 3 = 1983, 4 = 1984, 5 = 1985, 7 = 1987, 8 = 1988)
% yr_txt    = cell array of corresponding row labels

% Clear workspace:
clz;

% Load data:
load tikus.mat

% Apply a mild data transformation to down-weight the influence of the most
% abundant species:
coral_2 = f_normal(coral,'2');

% Create a Bray-Curtis dissimilarity matrix from the transformed data:
disBC = f_dis(coral_2,'bc'); % Bray-Curtis

% 1981 = before 
% 1983 = during
% 1985 = after disturbance event
% Create index to just these years:
idx = (yr==1 | yr==3 | yr==5);

% Test the null hypothesis that there is no difference in location (i.e., in the
% composition and abundance of species) among the 3 regions, irrespective of
% differences in multivariate dispersion:

% Quantify multivariate dispersion for each of these years:
f_permanova(disBC(idx,idx),yr(idx),1000,1);
% 
% ==================================================
%                Modified PERMANOVA:                
% --------------------------------------------------
%  F  = 5.3514    p    =  0.0010 (iter=1000) 
% --------------------------------------------------
% 
% -> reject the null hypothesis at alpha=0.05

% Perform pair-wise tests to see which groups differed:
f_permanovaPW(disBC(idx,idx),yr(idx),1000,1);
% 
% ----------------------------------------------------------
% Results of pair-wise comparisons between each factor level:
% ==========================================================
%          t:     p:     p_bon: p_ds:  p_holm:
% 1 vs. 2: 2.1980 0.0010 0.0030 0.0030 0.0030 
% 1 vs. 3: 2.7544 0.0010 0.0030 0.0030 0.0030 
% 2 vs. 3: 2.1232 0.0010 0.0030 0.0030 0.0030 
% ----------------------------------------------------------
% t      = t-statistic 
% p      = unadjusted p-value 
% p_bon  = Bonferroni adjusted p-value 
% p_ds   = Dunn-Sidak adjusted p-value 
% p_holm = Holms adjusted p-value 
% 
% 
% -> all pairs of groups were significantly different in terms of the location
% of their centroids, regardless of whether their dispersions differed or not.
