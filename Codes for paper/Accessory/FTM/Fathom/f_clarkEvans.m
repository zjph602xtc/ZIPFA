function result = f_clarkEvans(X,iter,verb,plt)
% - Clark-Evans index of spatial randomness
%
% USAGE: result = f_clarkEvans(X,iter,verb,plt)
%
% X    = two-column matrix of spatial coordinates
% iter = # iterations for permutation test                         (default = 0)
% verb = send results to display                                   (default = 1)
% plt  = optionally create a plot                                  (default = 0)
%
% result = structure of result with the following fields:
%  .CE     = Clark-Evans index
%  .p_perm = permutation-based p-value for CE
%  .Z      = Z-statistic
%  .p_para = parametric-based p-value for Z

% -----Notes:-----
% CE = 1: complete spatial randomness (CSR)
% CE > 1: uniform   (NN's are LARGER than expected)
% CE < 1: clustered (NN's are SMALLER than expected)
%
% NN's = nearest neighbor distances
%
% The Clark-Evans index uses the nearest neighbor distance (NN) among sites
% to evaluate the null hypothesis of complete spatial randomness (CSR). The
% index is based on the ratio of the observed average NN distance to that
% expected for a Poisson distribution (= CSR).
%
% The Z-stat is a test statistic based on (R_obs - R_exp)/SE, which provides a
% parametric method of evaluating the null hypothesis of CSR.

% -----References:-----
% Clarke, P. J. and F. C. Evans. 1954. Distance  to nearest neighbor as a
%   measure of spatial relationships in populations. Ecology 35: 445-453.
% Corral-Rivas, J. J., C. Wehenkel, H. A. Castellanos-Bocaz, B. Vargas-Larreta,
%   and U. Dieguez-Aranda. 2010. A permutation test of spatial randomness:
%   application to nearest neighbor indices in forest stands.

% -----Author:-----
% by David L. Jones, Nov-2014
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.
%
% Jan-2016: updated; permutation test now based on CE instead of Z following
% Corral-Rivas et al. (2010), which provides a more meaningful plot

% -----Set defaults & check input:-----
if (nargin <  2), iter = 0; end % default no permutation test
if (nargin <  3), verb = 1; end % send output to display by default
if (nargin <  4), plt  = 0; end % no plot by default
% -------------------------------------

% Get sample size:
n = size(X,1);

% Euclidean distance matrix:
dis = f_dis(X,'euc');

% Replace diagonal with NaN's:
dis(logical(eye(size(dis)))) = NaN;

% Get distance to nearest neighbor (NN):
nn = nanmin(dis);

% Observed average NN distance:
R_ave = mean(nn);

% Get index to convex hull:
[~,area] = convhull(X(:,1),X(:,2));

% Get expected density:
den = n/area;

% Expected average NN distance:
R_exp = 1/(2 * sqrt(den));

% Ratio of observed to expected:
CE = R_ave/R_exp;

% Z-stat:
SE = sqrt((4-pi)/(n*4*pi*den));
Z  = (R_ave - R_exp)/SE;

if CE>=1
   p_para = 1 - normcdf(Z);
else
   p_para = normcdf(Z);
end

% Test statistic:
% c = (R_ave - R_exp) / ( 0.2613 / (sqrt(n*(n/area))) );

% Optional permutation test:
if (iter>0)
   CE_perm = [CE;zeros(iter-1,1)]; % initialize
   for i=2:iter
      result_perm = f_clarkEvans(f_shuffle(X,3),0,0,0);
      CE_perm(i)  = result_perm.CE;                   % collected permuted CE's
   end
   
   % Determine which tail of the distribution to assess:
   if CE >= 1;
      p_perm = sum(CE_perm >= CE)/iter;
   else
      p_perm = sum(CE_perm <= CE)/iter;
   end
   
   % Plot histogram of permuted R-ratios:
   if (plt>0)
      figure;
      [bin,xBin] = hist(CE_perm,100); % get 100 bins
      bin  = bin/iter;               % convert absolute to relative
      h(1) = bar(xBin,bin);
      title('Frequency Histogram of Permuted \it{CE}''s')
      xTxt = sprintf('Values (n = %d)',iter);
      xlabel(xTxt)
      ylabel('Relative Frequency')
      box on;
      grid on;
      hold on;
      h(2) = plot(CE_perm(1),0,'ro','MarkerFaceColor','r'); % plot observed R
      txt = {'permuted \it{CE}' 'observed \it{CE}'};
      legend(h,txt);
   end
else
   p_perm = NaN;
end

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('      Clark-Evans - Complete Spatial Randomness:\n');
   fprintf('--------------------------------------------------\n');
   fprintf('CE index = %-3.4f  p =  %3.5f \n',CE,p_perm);
   fprintf('Z stat   = %-3.4f  p =  %3.5f \n',Z,p_para);
   fprintf('No. of permutations = %d \n',iter);
   fprintf('--------------------------------------------------\n\n');
end

% Wrap results up into a structure:
result.CE     = CE;     % Clark-Evans index
result.p_perm = p_perm; % permutation-based p-value for CE
result.Z      = Z;      % Z-statistic
result.p_para = p_para; % parametric-based p-value for Z
