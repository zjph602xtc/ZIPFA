function result = f_npDisp(yDis,grp,ctr,iter,verb,plt,simu)
% - homogeneity of multivariate dispersion (= NP-DISP, PERMDISP)
%
% USAGE: result = f_npDisp(yDis,grp,ctr,iter,verb,plt)
%
% yDis = square symmetric dissimilarity matrix derived from response variables
% grp  = column vector of integers specifying group membership for rows in yDis
% ctr  = use spatial median ('S') or centroid ('C')              (default = 'S')
% verb = optionally send results to display                        (default = 0)
% iter = # iterations for permutation test                         (default = 0)
% plt  = plot histogram of permuted statistics                     (default = 0)
%
% result = structure of results with the following fields:
%  .z       = residuals (= distance b/n each obs & its group's central tendency)
%  .grp     = grouping vector
%  .gLabels = group labels
%  .scores  = scaled PCoA eigenvectors
%  .evals   = eigenvalues
%  .C       = coordinates of each group's (row-wise) spatial median or centroid
%  .gMean   = mean of the residuals within each group
%  .F       = F-statistic
%  .p       = corresponding permutation-based p-value
%  .type    = 'spatial median' or 'centroid'
%
% SEE ALSO: f_npDispPlot

% -----Notes:-----
% Currently this code only supports 1 grouping variable
%
% SIMULATION RUN: (simu=1) returns only F statistic & p-value

% -----References:-----
% Anderson, M. J. 2006. Distance-based tests for homogeneity of multivariate
%   dispersions. Biometrics 62: 245-253
% Anderson, M. J., K. E. Ellingsen, and B. H. McArdle. 2006. Multivariate
%   dispersion as a measure of beta diversity. Ecology Letters 9(6): 683-693
% Anderson, M. J. and D. C. Walsh. 2013. PERMANOVA, ANOSIM, and the Mantel test
%   in the face of heterogeneous dispersions: What null hypothesis are you
%   testing? Ecological Monographs 83(4): 557-574. [see Sim_Functions.R]

% -----Author:-----
% by David L. Jones, Feb-2010
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Oct-2010: unique replaced with f_unique
% Feb-2011: uGrp, nGrp now defined earlier
% Dec-2011: removed 'Homogeneity' from documentation and table
% Oct-2012: moved last fprintf statement to proper position
% Feb-2015: replaced repmat with nan, added support for Cochran's test
% Jun-2016: overhaul: improve efficiency; perform significance test by permuting
%           the residuals

% -----Check input & set defaults:-----
if (nargin < 3), ctr  = 'S'; end % default use spatial median
if (nargin < 4), iter = 0;   end % default no permutation test
if (nargin < 5), verb = 0;   end % do not send output to display by default
if (nargin < 6), plt  = 0;   end % default no plot
if (nargin < 7), simu = 0;   end % internal flag: not a simulation run

% Check options during SIMULATION RUN:
if (simu>0)
   if (verb>0) || (plt>0)
      error('When SIMU=1, set VERB & PLT to 0!')
   end
end

% Check ctr:
ctr = lower(ctr); % force lower case
if ~isequal(ctr,'s') && ~isequal(ctr,'c')
   error('CTR must be ''S'' or ''C''')
end

n = size(yDis,1); % # observations

% Skip during a SIMULATION RUN:
if (simu<1)
   if n ~= size(grp,1), error('yDis & GRP need same # of rows'); end;
   
   if (f_issymdis(yDis) == 0)
      error('Input yDIS must be a square symmetric distance matrix');
   end
   
   if (~size(grp,2)==1)
      error('Only 1 grouping vector is currently supported!')
   end
end
% -------------------------------------

% Set up groups:
uGrp = unique(grp,'stable'); % unique groups, unsorted
nGrp = numel(uGrp);          % number of unique groups

% Get Principal Coordinates (scaled), keep negative eigenvalues:
temp  = f_pcoa(yDis,0,1,1);
U     = temp.scores;
evals = temp.evals;
clear temp;

% Get column indices & preallocate:
idxN  = (evals<0); % negative eigenvalues
idxP  = ~idxN;     % positive eigenvalues
neg   = any(idxN); % negative eigenvalues present (=1) or absent (=0)

% Preallocate:
z = NaN(n,1);
C = NaN(nGrp,size(U,2));

% -----Get residuals for each group:-----
for i = 1:nGrp
   idxG = (grp==uGrp(i)); % get index of rows of this group
   
   % Get spatial median or centroid of each group:
   if isequal(ctr,'s')
      C(i,:) = median(U(idxG,:));
   else
      C(i,:) = mean(U(idxG,:));
   end
   
   % Get Euclidean distance of each obs to its group center (= residuals):
   r_pos = sqrt(sum((U(idxG,idxP) - repmat(C(i,idxP),sum(idxG),1)).^2,2));
   
   % Handle negative eigenvalues:
   if (neg>0)
      % Get Euclidean distance of each obs to its group center (= residuals):
      r_neg = sqrt(sum((U(idxG,idxN) - repmat(C(i,idxN),sum(idxG),1)).^2,2));
      
      % Subtract IMAGINARY distance from REAL distance:
      z(idxG) = sqrt(abs(r_pos.^2 - r_neg.^2)); % ABS() avoids imaginary distances when r_neg > r
   else
      z(idxG) = r_pos;
   end
end
% ---------------------------------------

% Create ANOVA design matrix via dummy coding:
Rx = zeros(n,nGrp-1);
for i=1:(nGrp-1)
   Rx(grp==uGrp(i),i) = 1;
end

% Get observed F-statistic:
F = sub_anova(z,Rx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Permutation Test:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (iter>0)
   if (simu<1) % Show messages if not a simulation run
      fprintf('\nPermuting residuals %d times...\n',iter-1);
   end
   
   % Center scores witin each group to remove location effects (A&W, 2013):
   Uc = NaN(size(U)); % preallocate
   for i = 1:nGrp
      idxG = (grp==uGrp(i)); % get index of rows of this group
      Uc(idxG,:) = U(idxG,:) - repmat(C(i,:),sum(idxG),1);
   end
   
   F_perm = [F;zeros(iter-1,1)]; % initialize
   for j=2:iter
      idxS   = randperm(n)'; % get shuffled index to rows
      U_perm = Uc(idxS,:);   % permute centered scores
      
      % -----Get residuals for each group:-----
      z_perm = NaN(n,1); % preallocate
      C_perm = NaN(nGrp,size(U_perm,2));
      
      for i = 1:nGrp
         idxG = (grp==uGrp(i)); % get index of rows of this group
         
         % Get spatial median or centroid of each group:
         if isequal(ctr,'s')
            C_perm(i,:) = median(U_perm(idxG,:));
         else
            C_perm(i,:) = mean(U_perm(idxG,:));
         end
         
         % Get Euclidean distance of each obs to its group center (= residuals):
         r_pos = sqrt(sum((U_perm(idxG,idxP) - repmat(C_perm(i,idxP),sum(idxG),1)).^2,2));
         
         % Handle negative eigenvalues:
         if (neg>0)
            % Get Euclidean distance of each obs to its group center (= residuals):
            r_neg = sqrt(sum((U_perm(idxG,idxN) - repmat(C_perm(i,idxN),sum(idxG),1)).^2,2));
            
            % Subtract IMAGINARY distance from REAL distance:
            z_perm(idxG) = sqrt(abs(r_pos.^2 - r_neg.^2)); % ABS() avoids imaginary distances when r_neg > r
         else
            z_perm(idxG) = r_pos;
         end
      end
      % ---------------------------------------
      
      F_perm(j) = sub_anova(z_perm,Rx); % collected permuted F-ratios
   end
   p = sum(F_perm>=F)/iter;             % convert counts to a p-value
   
   % Optionally plot histogram of permuted F-ratios:
   if (plt>0)
      figure;
      [bin,xBin] = hist(F_perm,100); % get 100 bins
      bin = bin/iter;                   % convert absolute to relative
      h(1) = bar(xBin,bin);
      title('Frequency Histogram of Permuted F Statistics')
      xTxt = sprintf('Values (n = %d)',iter);
      xlabel(xTxt)
      ylabel('Relative Frequency')
      box on;
      grid on;
      hold on;
      h(2) = plot(F,0,'ro','MarkerFaceColor','r'); % plot observed delta value
      txt = {'permuted \it{F}' 'observed \it{F}'};
      legend(h,txt);
   end
else
   p = NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stop here for SIMULATION RUNS:
if (simu>0)
   result.F = F;
   result.p = p;
   return;
end

% Calculate group means:
gMean = NaN(nGrp,1); % preallocate
for i=1:nGrp
   idxG     = (grp==uGrp(i)); % get index to rows of this group
   gMean(i) = mean(z(idxG));
end

% Wrap results up into a structure:
result.z       = z;
result.grp     = grp;
result.gLabels = cellstr(num2str(uGrp)); % group labels, unsorted
result.scores  = U;
result.evals   = evals;
result.C       = C;
result.gMean   = gMean;
result.F       = F;
result.p       = p;
if isequal(ctr,'s')
   result.type = 'spatial median';
else
   result.type = 'centroid';
end

% -----Send output to display:-----
if (verb>0)
   fprintf('\n====================================================\n');
   fprintf(' NP-DISP: Homogeneity of Multivariate Dispersion  ');
   fprintf('\n----------------------------------------------------\n\n');
   fprintf('F = %3.2f  p = %-3.3f (iter=%d) \n\n',F,p,iter);
   fprintf('# Pos Eigenvalues = %d\n', sum(idxP));
   fprintf('# Neg Eigenvalues = %d\n\n', sum(idxN));
   txt = ['Average distance to ' result.type ':\n'];
   fprintf(txt);
   for i = 1:nGrp
      fprintf('Group %d = %2.4f\n',uGrp(i),gMean(i));
   end
   fprintf('-----------------------------------------------------\n');
end
% ---------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 SUBFUNCTION:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = sub_anova(z,Rx)
% - ANOVA F-statistic
% z  =  response variable
% Rx = ANOVA design matrix

Y  = z - mean(z); % center response variables
n   = size(Y,1);  % # observations
uno = ones(n,1);  % intercept term

% Calculate fitted values and residuals:
[Q,~] = qr([uno Rx],0); H = Q*Q'; % Hat-matrix for Rx
Yfit  = H*Y;                      % fitted values
Yres  = Y - Yfit;                 % residuals

% McArdle & Anderson (2001):
SSr = trace(Yfit'*Yfit); % SS_regression
SSe = trace(Yres'*Yres); % SS_error

% L&L,1998 (eq. 11.19); L,2007 (eq. 3):
dfR = size(Rx,2);                  % df_regression
F   = (SSr/dfR) / (SSe/(n-dfR-1)); % F-stat
