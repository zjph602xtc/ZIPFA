function result = f_permanova(yDis,grp,iter,verb)
% - one-way (modified) PERMANOVA
%
% USAGE: result = f_permanova(yDis,grp,iter,verb)
%
% yDis  = square symmetric dissimilarity matrix derived from response variables
% grp   = matrix of integers specifying group membership for objects in yDis
%         (column-wise)
% iter  = # iterations for permutation test  (default = 0)
% verb  = 1: send results to display; 2: also create plot          (default = 1)
%
% result.F  = F-statistic
% result.p  = permutation-based p-value
%
% See also: f_permanovaPW, f_npManova

% -----Notes:-----
% This function performs a one-way PERMANOVA using the modified pseudo-F
% statistic (F2) described in Anderson et al. (2017) to account for unbalanced
% designs with heterogeneous multivariate dispersions.

% -----References:-----
% Anderson, M. J. 2000. NPMANOVA: a FORTRAN computer program for non-parametric
%  multivariate analysis of variance (for any two-factor ANOVA design) using
%  permutation tests. Dept. of Statistics, University of Auckland.
%  (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
% Anderson, M. J. 2001. A new method for non-parametric multivariate
%   analysis of variance. Austral Ecology 26: 32-46.
% Anderson, M. J. 2002. DISTML v.2: a FORTRAN computer program to calculate a
%   distance-based multivariate analysis for a linear model. Dept. of Statistics
%   University of Auckland. (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
% Anderson, M. J., D. C. I. Walsh, K. R. Clarke, R. N. Gorley, and E.
%   Guerra-Castro. 2017. Some solutions to the multivariate Behrens-Fisher
%   problem for dissimilarity-based analyses. Aust. N. Z. J. Stat. 59(1): 57-79.
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.

% -----Author:-----
% by David L. Jones, May-2017
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), iter = 0; end % default iterations for permutation test
if (nargin < 4), verb = 1; end % send output to display by default

% Check size of input:
n = size(yDis,1);
if n ~= size(grp,1), error('yDis & GRP need same # of rows'); end
if ~(size(grp,2)==1), error('GRP must be a column vector!');  end
% -------------------------------------

% Gower's centered matrix:
A   = -0.5*(yDis.^2);
I   = eye(n,n);
uno = ones(n,1);
G   = (I-(1/n)*(uno*uno'))*A*(I-(1/n)*(uno*uno'));

% Get logical indices for each group:
uGrp   = unique(grp,'stable'); % get unsorted list of unique groups
nGrp   = numel(uGrp);          % # of groups
idxGrp = false(n,nGrp);        % preallocate
sGrp   = NaN(nGrp,1);
for i=1:nGrp
   idxGrp(:,i) = (grp == uGrp(i)); % get logical index to this group
   sGrp(i)     = sum(idxGrp(:,i)); % get size of this group
end

% Create Hat matrix:
[Q1,~] = qr([uno idxGrp(:,1:end-1)],0); H = Q1*Q1'; % compute Hat-matrix via QR

% Get within-group dispersions:
R_dia = diag((I-H)*G*(I-H)); % extract diagonals of residual G
V     = NaN(nGrp,1);         % preallocate
for i=1:nGrp
   V(i) = sum(R_dia(idxGrp(:,i))) / (sGrp(i)-1); % eq. 4 (Anderson et al. 2017)
end

% Modified F-statistic:
F = trace(H*G*H) / sum( (1-(sGrp/n)) .* V );     % eq. 2 (Anderson et al. 2017)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PERMUTATION TEST:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (iter>0)
   if (verb>0)
      fprintf('\nPermuting raw data under a full model %d times...\n',iter-1);
   end
   Fperm = [F;zeros(iter-1,1)]; % initialize
   for j = 2:iter
      idxS = randperm(n)'; % get shuffled index to rows/cols
      
      % Get within-group dispersions:
      R_dia = diag((I-H)*G(idxS,idxS)*(I-H)); % permuted residual diagonals
      V     = NaN(nGrp,1);                    % preallocate
      for i=1:nGrp
         V(i) = sum(R_dia(idxGrp(:,i))) / (sGrp(i)-1);
      end
            
      % Get permuted F-ratio:
      Fperm(j) = trace(H*G(idxS,idxS)*H) / sum( (1-(sGrp/n)) .* V );
   end
   p = sum(Fperm>=F)/iter; % convert counts to a p-value
   
   % Optionally plot histogram of permuted F-ratios:
   if (verb>1)
      figure;
      [bin,xBin] = hist(Fperm,100);  % get 100 bins
      bin  = bin/iter;               % convert absolute to relative
      h(1) = bar(xBin,bin);
      title('Frequency Histogram of Permuted F-ratios')
      xTxt = sprintf('Values (n = %d)',iter);
      xlabel(xTxt)
      ylabel('Relative Frequency')
      box on;
      grid on;
      hold on;
      h(2) = plot(F,0,'ro','MarkerFaceColor','r'); % plot observed delta value
      txt = {'permuted \it{F}' 'observed \it{F}'};
      legend(h,txt);
   end
   
else
   p = NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Send output to display:
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('               Modified PERMANOVA:                \n');
   fprintf('--------------------------------------------------\n');
   fprintf(' F  = %-3.4f    p    =  %-3.4f (iter=%d) \n',F,p,iter);
   fprintf('--------------------------------------------------\n');
end

% Wrap results up into a structure:
result.F = F;
result.p = p;
