function I = f_interaction(A,B,trim)
% - create an interaction term between factor A and factor B
% 
% USAGE: I = f_interaction(A,B,trim)
% 
% A    = column vector of factor A
% B    = column vector of factor B
% trim = trim last column to avoid a singular matrix (default = 1)
% 
% I = A*B interaction term
% 
% SEE ALSO: f_designMatrix

% -----Author:-----
% by David L. Jones, May-2017
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.


% -----Set defaults & check input:-----
if (nargin < 3), trim = 1; end % default to trimmed codes

% Force input as column vectors:
A = A(:);
B = B(:);

% Check size:
[nrA,ncA] = size(A);
[nrB,ncB] = size(B);

if ~isequal(nrA,nrB)
   error('A and B must have same #rows!')
end
% -------------------------------------

% Interaction Terms: (AxB):
I = []; % initialize
for i = 1:ncA
   I = [I repmat(A(:,i),1,ncB) .* B];
end
 
if trim>0 % trim last column to avoid a singular matrix
   I(:,end) = [];
end
