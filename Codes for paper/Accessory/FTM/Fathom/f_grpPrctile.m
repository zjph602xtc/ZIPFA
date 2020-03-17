function result = f_grpPrctile(X,grp,p)
% - returns the percentiles of X (column-wise) separately for groups
%
% USAGE: result = f_grpPrctile(X,grp,p);
%
% X   = column vector of input data
% grp = column vector of integers specifying group memberhip
% p   = percentiles
%
% 
% SEE ALSO: prctile, f_grpMean

% -----Author:-----
% by David L. Jones, June-2015
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

grp = grp(:); % col vector
p   = p(:)';  % row vector

if (size(X,2) > 1)
   error('X must be a column vector!')
end

if (size(X,1) ~= size(grp,1))
   error('# of rows in X and GRPS must be  equal !');
end

uGrp = f_unique(grp); % unique groups, unsorted
nGrp = numel(uGrp);   % # of groups
n    = numel(p);      % # of percentiles
prc  = nan(nGrp,n);   % preallocate


for i=1:nGrp
   idx      = (grp == uGrp(i));   % get indices for members group i
   prc(i,:) = prctile(X(idx,:),p);
end

% Rename for output:
result = prc;