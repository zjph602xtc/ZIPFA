function Y = f_filterSinclair(X,rep,win)
% - filter/smooth columns of a matrix via 11-point Moving Median + Average
%
% USAGE: Y = f_filterSinclair(X,rep,win);
%
% X = input matrix (rows = observations, cols = variables)
% rep = # of times to filter data                                  (default = 1)
% win = size of filter window                                     (default = 11)
% 
% Y = filtered/smoothed data
%
% SEE ALSO: f_filterMA

% -----Notes:-----
% This function uses the method of Sinclair et al. (1998) to filter/smooth
% time series data (or remove outliers) by first applying an 11-point
% moving median, then an 11-point moving average. This method is also
% recommended by Sanborn & Telmer (2003). Each column of data is smoothed
% separately by the same type of filter.

% -----References:-----
% Sanborn, M. and K. Telmer. 2003. The spatial resolution of LA-ICP-MS line
%  scans across heterogeneous materials such as fish otoliths and zoned
%  minerals. J. Anal. At. Spectrom. 18: 1231-1237.
% Sinclair, D. J., , L. P. J. Kinsley, & M. T. McCulloch. 1998. High
%  resolution analysis of trace elements in corals by laser ablation ICP-MS.
%  Geochim. Cosmochim. Acta 62:1889-1901.

% -----Author(s):-----
% by David L. Jones, Jun-2011
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin < 2), rep = 1;  end % don't repeat by default
if (nargin < 3), win = 11; end % default 11 point window size

% Set # points on either side of current point:
span = ceil((win-1)/2);
% -------------------------------------

Y = X;

for i = 1:rep
   
   % Filter (moving median):
   Y = f_filterMA(Y,span,1,0);
   
   % Smooth (moving average):
   Y = f_filterMA(Y,span,1,1);
end
