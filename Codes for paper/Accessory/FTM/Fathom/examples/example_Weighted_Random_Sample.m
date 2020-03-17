% Example of creating a weighted random sample
% 
% -----Author:-----
% by David L. Jones, Aug-2015
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.


X   = [3 6 9];              % list of all possible values
w   = [0.3 0.1 0.2];        % corresponding weighting factors
idx = f_wtRndSamp(10000,w); % index to random sample
R   = X(idx);               % apply index to create weighted random sample

% Tabulate results:
tabulate(R)
% 
%   Value    Count   Percent
%       1        0      0.00%
%       2        0      0.00%
%       3     5032     50.32%
%       4        0      0.00%
%       5        0      0.00%
%       6     1672     16.72%
%       7        0      0.00%
%       8        0      0.00%
%       9     3296     32.96%


% Note the weights normalized to sum = 1 are the same as the porportions of each
% value in weighted random sample:
(w./sum(w))'
% 
% ans =
% 
%           0.5
%       0.16667
%       0.33333
