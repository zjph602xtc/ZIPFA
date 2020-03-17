function idx = f_wtRndSamp(n,w)
% - get index to create a weighted random sample
% 
% USAGE: idx = f_wtRndSamp(n,w);
%  
% n = size of random sample
% w = row or col vector of weights to apply during the sampling procedure
% 
% idx = index to apply to the original data to create a weighted random sample
% 
% SEE ALSO: randi, rand

% -----Notes:-----
% If X represents the original data you wish to sample, the # of elements in 'w'
% should equal that of X.

% -----References:-----
% This function uses the Matlab 'HISTC' function to tally random numbers that
% occur within intervals of the uniform probability distribution that serve as
% weighting factors, see:
% http://stackoverflow.com/questions/2977497/weighted-random-numbers-in-matlab 

% -----Author:-----
% by David L. Jones, Aug-2015
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

w = w(:)';     % force as row vector
w = w./sum(w); % normalize weights to sum = 1:

% Get indices to weighted random sample of original data:
[~, idx] = histc(rand(n,1),cumsum([0 w]));
