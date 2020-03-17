function [bDis,bIdx] = f_bootDis(dis,grp)
% - bootstrap resampling with replacement of a symmetric distance matrix
%
% USAGE: [bDis,bIdx] = f_bootDis(dis,grp);
%
% dis = square symmetric distance matrix
% grp = vector of integers specifying group membership       (default = all 1's)
%
% bDis = bootstrapped distance matrix
% bIdx = bootstrapped indices
%
% SEE ALSO: f_boot, f_shuffle

% -----Author:-----
% by David L. Jones, Jul-2015
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults and check input:-----
if (nargin < 2), grp = []; end % default grp is empty

n = size(dis,1); % # of rows

if ~isequal(n,size(dis,2))
   error('DIS must be a square matrix!')
end

% Check grouping vector:
if isempty(grp)
   grp = ones(n,1);
else
   grp = grp(:); % force as col vector
end
if (n ~= numel(grp))
   error('# of rows in DIS and GRP must be equal !');
end
% -----------------------------------------

uGrp = unique(grp,'stable'); % get list of unique groups
nGrp = numel(uGrp);          % get # unique groups
bIdx = [];                   % initialize
% Bootstrap each group separately:
for i=1:nGrp
   idx        = find(grp == uGrp(i));                  % get index to rows of this group
   idxB       = sort(idx(unidrnd(numel(idx),numel(idx),1))); % bootstrap the indices
   dis(idx,:) = dis(idxB,:);                           % bootstrap rows
   dis(:,idx) = dis(:,idxB);                           % bootstrap cols
   bIdx       = [bIdx; idxB];                          % collect bootstrapped indices
end

% Rename for output:
bDis = dis;
