function y = f_shuffle(x,method,grp)
% - randomly sorts vector, matrix, or symmetric distance matrix
%
% Usage: y = f_shuffle(x,method,grp)
%
% -----Input/Output:-----
% x      = vector, matrix, or symmetric distance matrix
%
% method = type of permutation to perform
%          (default = 1 for a regular matrix or vector)
%          (default = 2 for a symmetric distance matrix)
%          1: unrestricted permutation
%          2: unrestricted, rows & cols are permuted the same way
%          3: permutation restricted to within columns of matrix
%          4: permute order of rows only (works across the matrix)
%          5: permutation restricted to within groups defined by grp
%          6: permutation within groups for symmetric distance matrix
%
% grp    = optional vector of integers specifying group membership for
%          restricted permutation
%
% y      = random permutation of x
%
% SEE ALSO: f_boot, f_bootDis, f_grpBoot, f_randRange

% -----Notes:-----
% When permuting a symmetric distance matrix, care must be taken to shuffle the
% objects making up the matrix and not the tridiagonal (see references below)
%
% To initialize RAND to a different state, call the following at the beginning
% of your routine (but not within a loop):
% >> rand('twister',sum(100*clock));


% -----References (for permutation of distance matrix):-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%  Elsevier Science BV, Amsterdam. xv + 853 pp. [page 552]
% Sokal, R. R. and F. J. Rohlf. 1995. Biometry - The principles and
%  practice of statistics in bioligical research. 3rd ed. W. H.
%  Freeman, New York. xix + 887 pp. [page 817]

% ----- Author: -----
% by David L. Jones, Mar-2002
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% 31-Mar-2002: added restricted permutation via grouping vector
% 18-Apr-2002: added switch-case handling of method options
%              and column-restricted permutation
% Dec-2007:    edited documentation
% Jan-2008:    changed & to &&; initialize RAND each time, changed 'find' to
%              logical indexing in 'case 5', preallocation of Y.
% Apr-2008:    Don't initialize RAND each time.
% Oct-2010:    Method 3 is now done internally
% Jun-2015:    Added method 6, re-wrote method 5

% -----Set Defaults:-----
if (nargin<2) && (f_issymdis(x)==0), method = 1; end; % default for non-symmetric matrices
if (nargin<2) && (f_issymdis(x)==1), method = 2; end; % default for symmetric matrices
if (nargin<3),                       grp = [];   end; % default GRP is empty

if (method==5) && (isempty(grp))
   error('Method 5 requires a grouping vector!')
end

if (method==6) && ( (f_issymdis(x)==0) || (isempty(grp)) )
   error('Method 6 requires a symmetric distance matrix & grouping vector!')
end

% Check GRP:
if ~isempty(grp)
   grp = grp(:); % force as column vector
   if ~isequal(size(x,1),size(grp,1));
      error('X & GRP are not of compatible sizes!');
   end
end
% -----------------------

[nr,nc] = size(x);
y = zeros(nr,nc); % preallocate

switch method
   case 1 % Permutation of a regular vector or matrix:
      y = x(randperm(length(x(:))));
      y = reshape(y,nr,nc);
      
   case 2 % Permutation of rows then colums,in the same way
      if (nr~=nc)
         error('Method 2 requires a square matrix')
      end
      i = randperm(nr); % get permuted indices
      y = x(i,:);       % permute rows
      y = y(:,i);       % permute cols
      
   case 3 % Permutation restricted to columns
      %    for i = 1:nc
      %       y(:,i) = f_shuffle(x(:,i),1);
      %    end
      for i = 1:nc
         y(:,i) = x(randperm(nr),i);
      end
      
   case 4 % Permute order of rows only (works across the matrix)
      i = randperm(nr); % get permuted indices
      y = x(i,:);       % permute rows
      
   case 5 % Permutation restricted to groups:      
      uGrp = unique(grp,'stable'); % unique groups
      nGrp = numel(uGrp);           % # of groups
      
      for i = 1:nGrp
         idx      = find(grp == uGrp(i));      % get index to rows of this group
         idxS     = idx(randperm(numel(idx))); % shuffle the indices
         y(idx,:) = x(idxS,:);
      end
                  
   case 6 % Permutation of symmetric distance matrix, restricted to groups:
      uGrp = unique(grp,'stable'); % get list of unique groups
      nGrp = numel(uGrp);          % get # unique groups
      
      % Process each group separately:
      for i=1:nGrp
         idx      = find(grp == uGrp(i));      % get index to rows of this group
         idxS     = idx(randperm(numel(idx))); % shuffle the indices
         x(idx,:) = x(idxS,:);                 % permute rows
         x(:,idx) = x(:,idxS);                 % permute cols
      end
      y = x; % rename for output:
   
   otherwise
      error('Unknown permutation method!');
end
