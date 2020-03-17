function result = f_isWhole(X)
% - determine whether X is a whole number
% 
% USAGE: result = f_isWhole(X)
% 
% X      = input
% result = 1: all values are whole numbers
%          0: one or more values are not whole numbers

% -----Author:-----
% by David L. Jones, Jun-2015
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

result = (mod(X,floor(X))==0);