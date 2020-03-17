function result = f_simper(X,grp,dis,verb,txt,lim)
% - similarity percentages (SIMPER) & variable importance
%
% USAGE: result = f_simper(X,grp,'dis',verb,txt,lim);
%
% X    = input matrix (rows = observations, cols = variables)
% grp  = column vector of whole numbers specifying group memberhip
%         e.g., grp = [1 1 2 2 1 1 2]';
% dis = 'bc' (Bray-Curtis dissimilarity) or 'euc' (Euclidean distance)
% verb = optionally send results to display  (default = 1)
% txt  = cell array of column labels; if empty, autocreate
%         e.g., txt = {'tem' 'sal' 'dep'};
% lim  = limit display to most important variables (default = empty)
%
% result = structure of results with the following fields:
%  .avDisTot = total average among-group distance
%  .avDis    = average among-group distance by variable
%  .avDis_SD = discrimination power of each variable
%  .expl     = percent contribution to avDisTot
%
% SEE ALSO: f_simper

% -----Notes:-----
% This function has been tested against Primer 6 for Windows and provides
% similar results.
%
% If LIM = 90, for example, only the most important variables will be displayed
% in the command window, i.e., only those variables that cumulatively explain up
% to 90% of the observed distances among groups.

% -----References:-----
% Clarke, K. R. and R. M. Warwick. 1994. Change in marine communities: an
%  approach to statistical analysis and interpretation. Natural Environment
%  Research Council, UK, 144 pp. [See page 7-3, eq. 7.1]

% -----Author:-----
% by David L. Jones, Apr-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Feb-2016: updated to calculate results for all possible pairs of groups
%           defined in GRP; now supports the Euclidean distance.

% -----Check input & set defaults:-----
if (nargin < 4), verb =  1; end % default send output to display
if (nargin < 5), txt  = f_num2cell((1:size(X,2))'); end % default labels
if (nargin < 6), lim  = []; end % default set display limit to 90%

dis = lower(dis); % force lowercase
if ~isequal(dis,'bc') && ~isequal(dis,'euc')
   error('DIS must specify ''bc'' or ''euc''!')
end

grp   = grp(:)'; % force as a row vector (so indexing below works)
[n,p] = size(X); % get # rows, # variables
if size(grp',1) ~= n, error('X and GRP must have same # of rows!'); end

% If labels are not cell arrays, try forcing them:
if iscell(txt)<1, txt = num2cell(txt); end

% Make sure labels are of compatible size:
txt = txt(:); % force cell array into a column
if size(txt,1) ~= p
   error('Size of TXT doesn''t match # of X variables!')
end

% Check limit:
if ~isempty(lim)
   if (lim <= 0 || lim >100)
      error('LIM must range from 0-100%!')
   end
end
% -------------------------------------

uGrp = unique(grp); % unique groups
nGrp = numel(uGrp); % # of unique groups

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 AMONG GROUPS:                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% List all possible pairs:
idxPW   = combnk(1:nGrp,2);   % get index to all pairwise combinations
PW      = uGrp(idxPW);        % get list of all pairwise combinations
[~,key] = sortrows(PW,[1 2]); % get sort key
PW      = PW(key,:);          % apply sort key to pairwise list
nr      = size(PW,1);         % get # rows of pairwise list

% Process each pair separately:
for j = 1:nr
   idx1 = find(grp==PW(j,1)); % get index to rows of group 1
   idx2 = find(grp==PW(j,2)); % get index to rows of group 2
   nr_1 = numel(idx1);        % get # rows in group 1
   blk  = numel(idx2);        % get # rows in group 2
   
   D    = NaN(nr_1*blk,1);    % preallocate
   d    = NaN(nr_1*blk,p);
   idx  = [0 0];              % initialize
   
   % Modified pair-wise dissimilaity, only inter- vs. intra-group site comparisons:
   for i = 1:nr_1 % repeat for each row in group 1
      r   = repmat(X(idx1(i),:),blk,1); % extract this row of group 1, replicate
      R   = X(idx2,:);                  % extract all rows of group 2
      idx = idx(end)+1:idx(end)+blk;    % get index to this block of distances
      switch dis
         case 'bc'  % Bray-Curtis Dissimilarity:
            D(idx)   = sum(abs(r-R),2)./ sum(r+R,2);      % dissimilarity among groups (all variables)
            d(idx,:) = abs(r-R)./ repmat(sum(r+R,2),1,p); % dissimilarity separate by variables
         case 'euc' % Euclidean Distance:
            D(idx)   = sum( (r-R).^2 ,2);    % distance among groups (all variables)
            d(idx,:) = (r-R).^2;             % distance separate by variable
      end
   end
   
   % Calculate statistics:
   avDisTot   = mean(D);                    % total average among-group dissimilarity
   avDis      = mean(d)';                   % average among-group dissimilarity by variables
   SD         = std(d)';                    % standard deviation by variable
   avDis_SD   = avDis./SD;                  % discrimination power of each variable
   expl       = NaN(p,1);                   % preallocate
   idxN       = ~isnan(SD);                 % prevent divide by NaN
   expl(idxN) = (avDis(idxN)/avDisTot)*100; % percent contribution to avDisTot
   
   % Sort by 'percent contribution:
   [~, keyS] = sort(expl,1,'descend'); % get sort key
   
   % Collect results in a structure:
   result.avDisTot{j} = avDisTot;
   result.avDis{j}    = avDis(keyS);
   result.avDis_SD{j} = avDis_SD(keyS);
   result.expl{j}     = expl(keyS);
   result.cum{j}      = cumsum(expl(keyS));
   result.var{j}      = txt(keyS);
end

% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   switch dis
      case 'bc'
         fprintf('Bray-Curtis Dissimilarity based  SIMPER:\n');
      case 'euc'
         fprintf('Euclidean Distance based SIMPER:\n');
   end
   fprintf('--------------------------------------------------\n');
   fprintf('avDis    = average among-group distance \n');
   fprintf('avDis_SD = discrimination power \n');
   fprintf('expl     = percent contribution (expl)\n');
   fprintf('cum      = cumulative percent contribution \n');
   fprintf('var      = variable\n');
   % fprintf('--------------------------------------------------\n');
   
   % Tabulate each pair separately:
   for j = 1:nr
      
      fprintf('\n\n==================================================\n');
      fprintf('  %d vs. %d: \n', PW(j,1),PW(j,2));
      switch dis
         case 'bc'
            fprintf('  (Total Average Dissimilarity = %1.4f)\n',result.avDisTot{j});
         case 'euc'
            fprintf('  (Total Average Distance = %1.4f)\n',result.avDisTot{j});
      end
      fprintf('==================================================\n');
      
      % Create Table:
      if isempty(lim)
         idxF = numel(result.expl{j});
      else
         idxF = find(result.cum{j}>=lim); % get index to last important variable
         idxF = idxF(1);                % select only the first
      end
      varName = {'avDis' 'avDis_SD' 'expl' 'cum' 'var'}; % define variable names
      T = table(result.avDis{j}(1:idxF),result.avDis_SD{j}(1:idxF),...
         result.expl{j}(1:idxF),result.cum{j}(1:idxF),result.var{j}(1:idxF),...
         'VariableNames',varName);
      fprintf('\n')
      disp(T) % display table
      fprintf('--------------------------------------------------\n');
   end
end
% ---------------------------------
