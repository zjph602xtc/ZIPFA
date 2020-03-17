% Example of creating custom box plots

% -----Author:-----
% by David L. Jones, Feb-2015
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Import data concerning benthic fauna from the Anclote River:
raw = f_importCSV('anclote.csv',1);

% Parse data:
grp = raw.dat(:,1);
tax = raw.dat(:,2);
den = raw.dat(:,3);
div = raw.dat(:,4);
clear raw;

% Create GRP labels:
txt = {'AM3' 'FF3' 'NF3' 'AM4' 'FF4' 'NF4' 'AM1' 'FF1' 'NF1' 'AM2' 'FF2' 'NF2'};

% Number of Taxa:
f_boxPlot(tax,grp,txt);

% Densities of individuals:
f_boxPlot(den,grp,txt);
% Force formatting of Y-axis labels:
set(gca,'Ytick',get(gca,'Ytick'),'YTickLabel',f_num2cell(get(gca,'Ytick')));

% Shannon Diversity:
f_boxPlot(div,grp,txt)
