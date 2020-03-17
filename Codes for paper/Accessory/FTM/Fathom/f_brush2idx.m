function [idx,XY] = f_brush2idx(ax)
% - create index to rows of data currently selected by 'brushing'
% 
% USAGE: idx = f_brush2idx(ax)
% 
% ax  = whole number indicating which axis to examine            (default = 1)
% idx = index to rows of brushed data
% dat = X,Y coordinates of brushed data
% 
% SEE ALSO: brush

% -----Author:-----
% by David L. Jones, Jul-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% May-2015: added an input parameter to specify which axis to examine, in case
%           the data in the figure have been plotted using separate calls to
%           'plot' (i.e., the figure is made up of multiple 'child' axes)

% -----Set defaults & check input:-----
if (nargin <  1), ax = []; end % default ax is empty
% 
% Determine # of axes:
n = numel(get(gca,'Children')); % % get # of axes

if (n>1) && isempty(ax)
   error('Please specify AX, since there are more than one!')
end

if (ax > n)
   error('AX specifies too many axes!')
end

if (ax < 1)
   error('AX specifies an invalid axis!')
end
% -------------------------------------

% Get index to brushed data for this axis:
try
   h   = findall(gca,'tag','Brushing');
   XY  = get(h, {'Xdata','Ydata'});
   idx = find(~isnan(XY{ax,1}))'; % Puts a '1' at indices that aren't NAN
catch
   error('Please BRUSH some data points first!')
end
