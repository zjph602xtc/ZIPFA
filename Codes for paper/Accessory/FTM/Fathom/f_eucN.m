function D = f_eucN(X,d,sym)
% - pair-wise Euclidean distances normalized to 0-1 scale
%
% X   = input matrix (row = observations, cols = variables)
% d   = maximum possible discrepancy             (if empty, calculate from data)
% sym = return distances as a square symmetric distance matrix     (default = 1)
%
% D   = normalized Euclidean distances
%
% SEE ALSO: f_dis

% -----Notes:-----
% If the maximum possible discrepancy (d) is the same for all columns of X,
% simply provide it as a scalar. If it varies by column, provide it as a row
% vector. If the value of 'd' is unknown, set d = [] and it will be calculated
% from the data.

% -----References:-----
% Barrett, P. 2005. Euclidean distance: raw, normalized, and double-scaled
%  coefficients. Advanced Projects R&D, The Technical Whitepaper Series 6; 26 pp.
%  http://www.pbarrett.net/techpapers/euclid.pdf [Accessed 12-Aug-2015]

% -----Author:-----
% by David L. Jones, Aug-2015
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set default and check input:-----
if (nargin <  2), d     = []; end % default set d to empty
if (nargin <  3), sym   = 1;  end % default symmetric distance matrix

% Get size of input:
[n,p] = size(X);

% Calculate maximum discrepancy for each column:
if isempty(d)
   d = max(X(:)) - min(X(:));
end

% Expand d to a row vector:
if isscalar(d)
   d = repmat(d,1,p);
end

% Check compatibility of d & X:
if size(d,2) ~= size(X,2)
   error('d  & X should have the SAME # of columns!')
end
if size(d,1)~=1
   error('d should be a scalar or row vector!')
end
% --------------------------------------

D   = NaN(n-1,1); % preallocate
idx = [0 0];      % initialize

% Scale data to 0-1 range:
X = X ./ repmat(d,n,1);

% Euclidean Distance (eq. 7.34/7.35 in Legendre & Legendre, 1998):
for i = 1:n-1 % repeat for each row except the last
   blk    = numel(i+1:n);              % get size of this block of distances
   r      = repmat(X(i,:),blk,1);      % extract this row, replicate
   R      = X(i+1:n,:);                % extract remaining rows
   idx    = idx(end)+1:idx(end)+blk;   % get index to this block of distances
   D(idx) = (sqrt( sum( (r-R).^2 ,2) ))./sqrt(p); % sum by row produces a column
end

% Optionally wrap distances up into a square symmetric matrix:
if (sym>0)
   D = f_rewrap(D);
end
