function result = f_mvabund(Y,grp,verb)
% - GLM-based analysis of multivariate abundance data
%
% USAGE: result = f_mvabund(Y,grp,verb);
%
% Y    = matrix of response variables (rows = obs, cols = species)
% grp  = column vector of integers specifying group membership for rows in Y
% verb = send output to display                                    (default = 0)
% 
% result = structure of results with the following fields:
%  .LR = sum of likelihood ratios
%  .p  = p-value

% -----References:-----
% Wang, Y., U. Naumann, S. T. Wright, and D. I. Warton. 2012. mvabund - an R
%  package for model-based analysis of multivariate abundance data. Methods in
%  Ecology & Evolution 3: 471-474.

% -----Author:-----
% by David L. Jones, May-2016
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.


% -----Set defaults & check input:-----
if (nargin < 3), verb = 0;    end % default don't send output to display

[nr,nc] = size(Y); % # obs, # variables

grp = grp(:); % force as column vector
if size(Y,1) ~= nr, error('Y & GRP need same # of rows'); end
% -------------------------------------

% Autocreate variable labels:
Y_txt = f_num2cell(1:nc)';

% Define working directory:
DIR = '/Users/djones/Dropbox/djones/work/mWork/Fathom/mvabund/';

% Export data for R:
row_txt = f_num2cell(1:size(Y,1))';                   % autocreate row labels
f_exportR([grp Y],row_txt,['grp' Y_txt],[DIR 'mvabund_IN.csv'])

% Build command script:
cmd = ['R CMD BATCH ' DIR 'run_mvabund.R ' DIR 'R_log.txt'];

% Run command script:
system(cmd);

% Load results:
temp = load([DIR 'mvabund_OUT.dat']);
LR   = temp(1);
p    = temp(2);

% -----Display table of results:-----
if (verb>0)
   fprintf('\n==================================================\n');
   fprintf('          mvabund - GLM-based GEE Test:           \n');
   fprintf('--------------------------------------------------\n');
   fprintf(' LR = %-3.4f  p = %-3.4f (iter=%d) \n',LR,p,1000);
   fprintf('--------------------------------------------------\n');
   fprintf('\n')
end

% Wrap results up into a structure:
result.LR = LR;
result.p  = p;
