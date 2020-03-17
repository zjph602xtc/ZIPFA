function result = f_mantel(yDis,xDis,rank,iter,plt,simu)
% - standardized Mantel statistic for 2 symmetric distance matrices
%
% Usage: result = f_mantel(yDis,xDis,rank,iter,plt);
%
% yDis = symmetric distance matrix (permuted in randomization test)
% xDis = symmetric distance matrix (MODEL MATRIX in hypothesis tests)
% rank = 1: rank distances in yDis, 2: rank both                     (default = 0)
% iter = number of iterations to use for 1-tailed randomization test (default = 0)
% plt  = plot histogram of permuted r-values                         (default = 0)
%
% result = structure of results with the following fields:
%  .r = standardized Mantel statistic
%  .p = permutation-based p-value
%
% -----Notes:-----
% A) xDis can be a model matrix for hypothesis testing.
% B) The permutation test shuffles the objects making up yDis.
% C) yDis & xDis must be derived INDEPENDENTLY.
% D) Don't rank xDis if it's a model matrix.
%
% SEE ALSO: f_modelMatrix, f_anosim, f_bioenv, f_procrustes, f_correlogram

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam. xv + 853 pp. (page 552)

% -----Author:-----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 18-Mar-02: make proper call to f_shuffle for distance matrix & corrected
%            calculation of p-value
% 27-Mar-02: made x the permuted matrix
% May-2008: changed | to ||
% Mar-2011: f_transform replaced with f_stnd
% Feb-2015: now optionally plots permutation distribution
% Apr-2016: edited code for readability; no longer calls f_shuffle; added
%           support for simu; results now returned as structure; rank option
%           more flexible

% -----Set defaults & check input:-----
if (nargin < 3), rank = 0; end % don't rank by default
if (nargin < 4), iter = 0; end % don't perform randomization test
if (nargin < 5), plt  = 0; end % default no plot
if (nargin < 6), simu = 0; end % internal flag: not a simulation run

% Get # rows:
nr = size(yDis,1);

% Skip during a SIMULATION RUN:
if (simu==0)
   % Check for symmetric distance matrices:
   if (f_issymdis(yDis) == 0) || (f_issymdis(xDis) == 0)
      error('yDis & xDis must be square symmetric distance matrices!');
   end
   
   % Check input:
   if ~isequal(nr,size(xDis,1))
      error('yDis & xDis must be equal in size!')
   end
   
else
   % Check plot options:
   if (plt>0)
      error('When SIMU=1, set PLT to 0!')
   end
end
% -------------------------------------

% Optionally rank distances:
switch rank
   case 0
      % don't rank
   case 1
      yDis = f_ranks(yDis);
   case 2
      yDis = f_ranks(yDis);
      xDis = f_ranks(xDis);
   otherwise
      error('RANK must be 0, 1, or 2!')
end

% Unwrap lower tridiagonal, standardize:
yVec = f_stnd(f_unwrap(yDis,0));
xVec = f_stnd(f_unwrap(xDis,0));
n    = numel(yVec); % get number of elements

% Take sum of cross-products, divide by n-1:
r = (sum(yVec.*xVec))/(n-1);

% -----Permutation Test:-----
if (iter>0)
   perm_r = [r; NaN(iter-1,1)]; % preallocate
   
   for i = 2:iter % observed value is considered a permutation
      idxS      = randperm(nr)'; % get index to shuffled rows/cols
      perm_r(i) = (sum(f_stnd(f_unwrap(yDis(idxS,idxS),0)).*xVec))/(n-1); % collect randomized stat
   end
   
   if (r>=0)
      p = sum(perm_r >= r)/iter; % right-tail test
   else
      p = sum(perm_r <= r)/iter; % left-tail test
   end
else
   p = NaN;
end
% ---------------------------

% Optionally plot histogram of permuted stats:
if (plt>0)
   figure;
   [bin,xBin] = hist(perm_r,100); % get 100 bins
   bin  = bin/iter;               % convert absolute to relative
   h(1) = bar(xBin,bin);
   title('Frequency Histogram of Permuted \it{r}')
   xTxt = sprintf('Values (n = %d)',iter);
   xlabel(xTxt)
   ylabel('Relative Frequency')
   box on;
   grid on;
   hold on;
   h(2) = plot(r,0,'ro','MarkerFaceColor','r'); % plot observed delta value
   txt = {'permuted \it{r}' 'observed \it{r}'};
   legend(h,txt);
end

% if (nargout==0)
%    fprintf('\nr = %3.4f \np = %3.4f (%3.0f iterations) \n',r,p,iter);
% end

% Wrap results up into a structure:
result.r = r;
result.p = p;
