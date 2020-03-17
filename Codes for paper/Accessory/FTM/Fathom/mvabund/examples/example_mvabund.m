% Examples of running the 'mvabund for R' package from MATLAB
% 
% -----Author:-----
% by David L. Jones, May-2016
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Gee et al. (1985)                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Follow Warton's (2011) analysis of nematode data from Gee et al. (1985)
% 
% Data from Gee et al. (1985) provided by the mvabund package for R as 'solberg'
% Y     = abundances of 53 taxa of nematodes from 12 experimental ponds
% Y_txt = corresponding column labels
% grp   = nutrient levels (1 = control, 2 = low, 3 = high)

% Clear workspace:
clz

% Import data:
raw = f_importCSV('solberg.csv',1);

% Parse data:
Y     = raw.dat(:,2:end);
Y_txt = raw.txt(2:end);
grp   = raw.dat(:,1);
clear raw;
fname = 'solberg';
% saver

% Get index to taxa that occur in at least 4 ponds
idx = (sum(Y>0))>=4;

% Show # taxa that match criterion:
sum(idx)
% ans = 26

% Remove taxa that don't match criterion
Y(:,~idx)   = [];
Y_txt(~idx) = [];

% Show data
f_table([grp Y],'%d','t')
% 
% 1	0	0	1	38	0	0	0	4	0	3	0	11	0	0	2	9	2	0	1	4	1	2	2	11	0	2
% 1	1	0	1	8	0	0	3	7	1	0	1	9	2	0	17	4	4	0	1	4	1	0	0	27	1	3
% 1	2	0	0	7	0	1	5	14	0	1	1	12	1	1	7	9	11	1	1	2	0	0	0	21	0	1
% 1	1	0	1	6	0	1	4	23	1	1	2	8	1	3	4	2	8	1	2	2	1	0	4	16	0	2
% 2	2	1	0	10	9	0	3	6	1	2	1	12	2	0	5	8	2	0	0	8	1	0	0	13	2	4
% 2	0	2	0	21	4	1	0	11	1	0	0	2	3	0	10	7	18	0	1	6	0	0	3	4	0	3
% 2	1	1	1	16	1	6	9	5	1	0	3	8	2	1	2	2	7	1	2	11	0	1	0	9	1	4
% 2	1	1	0	26	0	2	9	6	0	0	0	6	0	0	18	6	7	0	0	3	0	0	0	8	0	6
% 3	1	0	0	4	0	0	0	7	1	2	0	26	0	0	7	7	4	0	0	0	0	0	0	29	3	3
% 3	0	1	1	17	0	2	0	17	3	1	1	14	4	0	6	7	8	0	0	0	0	0	1	14	0	2
% 3	2	0	4	4	2	0	0	4	0	7	0	29	2	2	2	17	4	2	0	0	0	2	0	9	0	4
% 3	0	0	3	32	0	0	0	1	0	1	4	0	0	0	16	1	11	3	4	0	0	1	0	7	3	8
%  
% -> this matches Table 1 in Warton (2011)


% Use 'mvabund for R' to perform a GEE Test:
result = f_mvabund(Y,grp,1);

% ==================================================
%           mvabund - GLM-based GEE Test:           
% --------------------------------------------------
%  LR = 92.5898  p = 0.0460 (iter=1000) 
% --------------------------------------------------
