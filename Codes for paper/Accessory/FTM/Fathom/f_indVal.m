function result = f_indVal(Y,grp,iter,verb,txt,alpha)
% - species indicator values analysis (IndVal)
%
% USAGE: result =  f_indVal(Y,grp,iter,verb,txt,alpha);
%
% Y    = matrix of response variables (rows = obs, cols = species)
% grp  = column vector of integers specifying group membership for rows in Y
% iter = # iterations for permutation test                          (default = 0)
% verb = display all (=1), significant (=2), or sorted (=3) results (default = 0)
% txt  = cell array of Y labels; if empty, autocreate
%        e.g., Ytxt = {'sp1' 'sp2' 'sp3'};
% alpha = significance level                                    ( default = 0.05)
%
% result = structure of results with the following fields:
%  .indVal = species indicator power values (ranges from 0-100%)
%  .p      = corresponding permutation-based p-value
%  .IV     = all indicator values
%  .G      = for each species (col), the group associated with each indVal
%  .txt    = Y labels
%  .col    = index to columns of Y for significant indicator values

% -----Notes:-----
% A:      specificity is maximum when a species occurs in only 1 group
% B:      fidelity is maximum when a species occurs in all sites defining a group
% indVal: the indicator value is maximum (= 1) when a species occurs in
%         only one group and is present in all sites comprising that group.
%
% An indicator species is a taxon that is characteristic of a group of
% sites because: (1) it mostly occurs only in one group; and (2) occurs in
% most of the sites comprising that group.
%
% indVal = this metric is a measure of the magniture of a species ability to
% indicate (characterize) groups of sites, while the corresponding p-value
% is a measure of the statistical significance of that metric.
%
% IV = this matrix simply returns all the indicator values used to
% determine maximum value used to determine indVal.
%
% G = this matrix provides the actual group designation associated with the
% maximum indicator value returned by indVal
%
% Note: Legendre & Legendre (2012) suggest interpreting p-values with
% caution if the same taxa are used to both: (1) define the groups and (2)
% caluclate species indicator values.
%
% De Caceres & Legendre (2009) suggest using sqrt(IndVal)

% -----References:-----
% De Caceres, M. and P. Legendre. 2009. Associations between species and groups
%   of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
% Dufrene, M. and P. Legendre. 1997. Species assemblages and indicator
%   species: the need for a flexible asymmetrical approach. Ecol. Monogr.
%   67(3): 345-366.
% Legendre, P., and L. Legendre. 2012. Numerical ecology, 3rd English
%   edition. Developments in Environmental Modelling, Vol. 24. Elsevier
%   Science BV, Amsterdam. xiv + 990 pp.

% -----Author:-----
% by David L. Jones, Apr-2013
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2016: added support for verb and txt; removed call to f_shuffle;
% Apr-2016: p-values no longer returned for IV (just indVal) following De Caceres
%           et al. (2009); S is now calculated by max, but only G is returned;
% May-2016: checks for negative presence or abundance values
% Mar-2017: now optimally sorts output (descending) by indVal when verb=3

% -----Set defaults & check input:-----
if (nargin < 3), iter  = 0;    end % no permutation test by default
if (nargin < 4), verb  = 0;    end % default don't send output to display
if (nargin < 5), txt   = [];   end % default Y labels
if (nargin < 6), alpha = 0.05; end % default significance level

[nr,nc] = size(Y); % # obs, # variables

grp = grp(:); % force as column vector
if nr ~= size(grp,1), error('Y & GRP need same # of rows'); end

% Check for negative presence or abundance values:
if any(Y<0)
   error('Y cannot contain negative presence or abundance values!')
end

% Autocreate variable labels:
if isempty(txt)
   txt = f_num2cell(1:nc)';
end

% Check labels:
txt = txt(:)'; % force as row vector
if size(txt,2) ~= nc
   error('# cols of TXT must match Y!')
end

% Check ITER and VERB:
if (iter==0) && (verb>0)
   error('Specify a value for ITER for VERB>0!)')
end
% -------------------------------------

% Process groups:
uGrp = unique(grp,'stable'); % unique groups, unsorted
nGrp = length(uGrp);         % # unique groups
if (iter>0 && nGrp==1)
   error('There''s only 1 group, re-run with iter=0!')
end

% Calculate Specificity & Fidelity:
Yb = double(Y>0);  % convert to binary (presence/absence)
A  = NaN(nGrp,nc); % preallocate
B  = A;
for i = 1:nGrp % repeat for each group
   idx    = (grp==uGrp(i));                     % get index to rows of this group
   n      = sum(idx);                           % # sites in this group
   A(i,:) = mean(Y(idx,:),1);                   % mean abundance of this group
   B(i,:) = sum(Yb(idx,:),1) ./ repmat(n,1,nc); % relative frequency within group
end

% Convert to relative abundance across groups:
A = A ./ repmat(sum(A,1),nGrp,1);

% Calculate indicator values:
IV = (A .* B)*100;

% Replace NaN's with 0:
IV(isnan(IV)) = 0;

% Keep maximum values:
[indVal,S] = max(IV,[],1);

% Identify group corresponding to indVal:
G = uGrp(S)';

%-----Permutation test:-----
if iter>0
   fprintf('\nPermuting the data %d times...\n',iter-1);
   
   perm_indVal = [indVal; NaN(iter-1,nc)];         % preallocate with observed statistic
   for i = 2:iter
      idxS             = randperm(nr)';             % get shuffled index to rows
      temp             = f_indVal(Y(idxS,:),grp,0); % permute obs (rows)
      perm_indVal(i,:) = temp.indVal;               % permuted indicator values
   end
   
   % Get permuted stats >= to observed stat
   p = sum((perm_indVal - repmat(indVal,iter,1)) >= 0)/iter; % convert counts to probability
else
   p = NaN;
end
%-----------------------------

% Index columns of significant indicators:
col = find(p<=alpha);

% Sort species (desending) by indicator power values:
[~,key] = sort(indVal,2,'descend'); % get sort key

% -----Display table of results:-----
if (verb>0)
   if (verb==3)
      idx = key;
   elseif (verb>1) && ~isempty(col)
      idx = col; % index to indicator species
   else
      idx = 1:nc; % index to all species
   end
   
   % Create a table:
   varName = {'indVal' 'G' 'p' 'txt' 'col'}; % define variable names
   T = table(indVal(idx)', G(idx)', p(idx)', txt(idx)',idx', 'VariableNames', varName);
   
   % Display Table:
   % -----Send output to display:-----
   fprintf('\n==================================================\n');
   if (verb>1) && ~isempty(col)
      fprintf('Significant Indicator Species Values (IndVal):\n');
   else
      fprintf('Indicator Species Values (IndVal):\n');
   end
   fprintf('--------------------------------------------------\n');
   fprintf('indVal = indicator power values (ranges from 0-100%%) \n');
   fprintf('G      = corresponding group \n')
   fprintf('p      = p-value \n');
   fprintf('txt    = species labels\n');
   fprintf('col    = index to columns of Y\n')
   fprintf('--------------------------------------------------\n');
   fprintf('\n')
   disp(T) % display table
end
% -----------------------------------

% Wrap results up into a structure:
result.indVal = indVal;
result.p      = p;
result.IV     = IV;
result.G      = G;
result.txt    = txt;
if (verb==3)
   result.col = key;
else
   result.col    = col; % index to columns of significant indicators
end
