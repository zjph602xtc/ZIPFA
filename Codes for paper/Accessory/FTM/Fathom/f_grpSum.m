function gSum = f_grpSum(x,grp)
% - returns the sum of X (column-wise) separately for groups
%
% USAGE: gSum = f_grpSum(x,grp);
%
% x   = column vector or matrix of input data
% grp = column vector of integers specifying group memberhip
%
% gSum  = sum of each group (row-wise) for each variable in X (col-wise)
%
% SEE ALSO: f_grpMean, f_grpRel, f_grpOutlier, f_grpSize

% -----Notes:-----
% This function is used to calculate the sum of each variable in matrix X
% (column-wise) separately for each group specified by GRP. The number of
% rows in gSum will equal the number of unique groups specified in GRP.

% -----Author:-----
% by David L. Jones, March-2015
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

grp     = grp(:); % col vector
[nr,nc] = size(x);

if (nr ~= size(grp,1))
   error('# of rows in X and GRP must be  equal !');
end

uGrp = unique(grp,'stable'); % unique groups, unsorted
nGrp = numel(uGrp);          % # of groups

gSum = zeros(nGrp,nc); % preallocate
for i=1:nGrp
   idx = grp == uGrp(i);         % get indices for members group i
   gSum(i,:) = nansum(x(idx,:)); % SUM for members of this group
end
