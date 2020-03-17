function f_boxPlot(X,grp,txt)
% - create customized box plot
%
% USAGE: f_boxPlot(X,grp,txt)
%
% X   = column vector of data
% grp = grouping variable
% txt = cell array of unique group labels (default = autocreate)
% 
% SEE ALSO: f_grpPrctile

% -----NOTES:-----
% This function is used to construct a custom formatted box and whisker plot.
% The upper and lower edges of the red  box represents the 25th and 75th
% percentiles (= first/lower and third/upper quartiles). The vertical height of
% the box represents the 'interquartile range' (= difference between 1st and 3rd
% quartile). The horizontal black line within the red box represents the median.
% The whiskers represent the highest and lowest data points not considered
% outliers or extreme values. Black circles represent outliers, which are values
% +/-  1.5 to 3 times the interquartile range from the upper and lower
% quartiles. Red triangles represent data points that are more than 3 times the
% interquartile range from the upper and lower quartiles.

% -----Author:-----
% by David L. Jones, Jun-2015
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Check input & set defaults:-----
if (nargin <  3), txt = []; end % autogenerate labels

% Force column vectors:
X   = X(:);
grp = grp(:);

% Check size of input:
if numel(X) ~= numel(grp)
   error('X & GRP must have same # rows!')
end
% -------------------------------------

% Process groups:
uGrp = f_unique(grp); % get unique groups
nGrp = numel(uGrp);   % # of groups

% Autocreate group labels:
if isempty(txt)
   txt = f_num2cell((1:nGrp)');
end

% Create box plot:
H = boxplot(X,grp,'plotstyle','traditional','boxstyle','outline','colors','k',...
   'symbol','r.');

% Get percentiles:
P  = f_grpPrctile(X,grp,[25 75]);
q1 = P(:,1); % 1st (= lower) quartile
q3 = P(:,2); % 3rd (= upper) quartile

% Plot outliers & extremes:
hold on;

for i=1:nGrp
   
   % Extract data for this group:
   dat = X(grp==i);
   
   % Get index to outliers (after help in BOXPLOT.M):
   out = ( (dat>=(q3(i) + 1.5*(q3(i) - q1(i)))) & (dat<=(q3(i) + 3*(q3(i) - q1(i)))) ) |...
      ( (dat<=(q1(i) - 1.5*(q3(i) - q1(i)))) & (dat>=(q1(i) - 3*(q3(i) - q1(i)))) ) ;
   
   % Get index to extremes:
   ext = ( (dat>(q3(i) + 3*(q3(i) - q1(i)))) ) |...
      ( (dat<(q1(i) - 3*(q3(i) - q1(i)))) ) ;
   
   % Convert logical to row indices:
   out = find(out);
   ext = find(ext);
   
   % Plot outliers
   if ~isempty(out)
      plot(repmat(i,numel(out),1),dat(out),'ko','MarkerFaceColor','k')
   end
   
   % Plot extremes:
   if ~isempty(ext)
      plot(repmat(i,numel(ext),1),dat(ext),'rv','MarkerFaceColor','r')
   end
end

% -----Customize plot:-----
%
% Set group labels:
set(gca,'Xtick',1:nGrp,'XTickLabel',txt)

% Change line style of whiskers:
set(H(1:2,:),'LineStyle','-')

% Add filled box:
for i=1:nGrp
   xBox = get(H(5,i),'Xdata'); % get x coordinates of box
   yBox = get(H(5,i),'Ydata'); % get y coordinates of box
   fill(xBox,yBox,'r');        % plot box filled with red
   
   xMed = get(H(6,i),'Xdata'); % get x coordinates of median
   yMed = get(H(6,i),'Ydata'); % get y coordinates of median
   plot(xMed,yMed,'k-');       % replot median to show over box   
end
