function result = f_rda(Y,X,W,iter,verb,stnd,perm,simu,grp)
% - Redundancy Analysis (RDA)
%
% USAGE: result = f_rda(Y,X,W,iter,verb,stnd)
%
% Y    = matrix of response variables (rows = sites, cols = species)
%
% X    = (1) matrix of explanatory variables, or
%        (2) ANOVA design matrix specified by dummy coding
%
% W    = matrix of covariables                                        (0 = none)
% iter = # iterations for permutation test                         (default = 0)
% verb = 1: send results to display; 2: also create plot           (default = 1)
% stnd = optionally standardize X & W                              (default = 1)
%
% result = structure returning the following values:
%  .F          = F statistic
%  .p          = p-value for F stat
%  .R2         = Coefficient of determination (R^2)
%  .R2adj      = adjusted R^2
%  .R2w        = R^2 of covariate
%  .dfR        = df regression
%  .dfE        = df error
%  .dfT        = df total
%  .SSr        = SS regression
%  .SSr_W      = SS regression of covariable
%  .SSe        = SS error
%  .SSt        = SS total
%  .MSr        = MS regression
%  .MSe        = MS error
%  .X          = standardized explanatory variables
%  .Y          = centered response variables
%  .fit        = fitted values of Y
%  .res        = residual values of Y
%  .siteScores = site scores        (= WA scores)
%  .fitScores  = fitted site scores (= LC scores)
%  .resScores  = residual scores
%  .U          = canonical eigenvectors
%  .Ures       = residual eigenvectors
%  .evals      = canonical eigenvalues
%  .evalsRes   = residual eigenvalues
%  .canVar     = fraction of variance explained by canonical axes (= R^2)
%  .canCorr    = Species-Environ correlations (r) of each canonical axis
%  .B          = Matrix of regression coefficients produced by Y\X
% SEE ALSO: f_rdaPlot, f_rdaDB, f_rdaAnova, f_rdaAIC, f_rdaStepwise, f_npManova

% -----Notes:-----
% If X or W are a design matrix of dummy codes set STND=1 and do not include an
% intercept term, otherwise the F-ratio will be incorrect.
%
% To permute residuals: specify + value for ITER
% To permute raw data;  specify - value for ITER
%
% If Y is a column vector (= single response variable), RDA is simply multiple
% linear regression analysis and the PCA steps are skipped.
%
% When a matrix of covariates is provided (W), the R2 and R2adj provide the
% 'semipartial' versions of these statistics.

% -----Scaling:-----
% In this implementation of RDA the eigenvectors are scaled to lengths of 1,
% creating DISTANCE biplots which preserve the distances among sites. These are
% interpreted as follows:
%
%     A) Distances among sites approximate their Euclidean distance.
%     B) Projecting sites onto a Y arrow approximates the value of that
%        site along the variable.
%     C) Angles among Y variables are meaningless.
%     D) Angles b/n Y and X arrows reflect their correlations.
%     E) Distances among centroids and between centroids and sites
%        approximate Euclidean distance.

% -----RUN MODES:-----
% PERMUTATION RUN: (perm=1) returns only basic statistics so the calling
%                  function can calculate permutation based p-values
%
% SIMULATION RUN:  (simu=1) returns only basic statistics & p-value

% -----RESTRICTED PERMUTATION:-----
% grp: specifies the groups that permutations should be restricted within

% -----References:-----
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% Legendre, P. 2007. Studying beta diversity: ecological variation partitioning
%   by multiple regression and canonical analysis. Journal of Plant Ecology,
%   doi:10.1093/jpe/rtm001
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.
% Peres-Neto, P. R., P. Legendre, S. Dray, & D. Borcard. 2006. Variation
%   partitioning of species data matrices: estimation and comparison of
%   fractions. Ecology 87(10): 2614-2625.
% Schielzeth, H. 2010. Simple means to improve the interpretability of
%   regression coefficients. Methods in Ecology & Evolution 1(2): 103-113.

% -----Author:-----
% by David L. Jones, Mar-2003
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Oct-2003: updated documentation, fixed 'Check scaling', added verbPerm to
%           reduce clutter when run by f_rdaStepwise
% Jan-2008: updated documentation
% Feb-2008: updated documentation; changed | to ||, changed & to &&;
%           preallocated corrX_fit, corrY_fit,canonCorr; num2Str canged to
%           num2str; moved .corrX, .corrY, and .biplotX and all plotting
%           routines to 'f_rdaPlot' (as a result the 5th command-line option has
%           been removed so adjust legacy code calling this)
% Mar-2008: overhaul: F-stat relies on SS, not permuted eigenvectors;
%           minimized calculations during permutation runs; trim residual axes
%           with 0 variance; added R2adj; added ANOVA stats; PCA is skipped for
%           univarite response data; make sure W ~= X
% May-2008: added SSe, SSt to output when perm>0
% Mar-2013: modified to be compatible with changes made to f_pca
% May-2014: now outputs scores and canVar for univariate data for plotting
% Feb-2016: permute residuals or raw data; optional plot of permuted F stats;
%           now returns R2w; updated permutation test; no longer calls
%           f_shuffle; updated the way residual axes with variance = 0 are
%           trimmed; also trims residual scores
% Mar-2016: 'stnd' option now allows you to NOT standardize X & W
% Apr-2016: verpPerm replaced with SIMU; added 'simulation mode';  no covariate
%           now specified by W=[] instead of W=0
% May-2017: now supports restricted permutations

% -----Changes made by Joshua P. Kilborn-----
% Feb-2018: Added the B matrix to the output results structure

% -----Set defaults & check input:-----
if (nargin < 3), W    = []; end % no covariables by default
if (nargin < 4), iter = 0;  end % no permutation test by default
if (nargin < 5), verb = 1;  end % send output to display by default
if (nargin < 6), stnd = 1;  end % default standardize X & W
if (nargin < 7), perm = 0;  end % internal flag: not a permutation run
if (nargin < 8), simu = 0;  end % internal flag: not a simulation run
if (nargin < 9), grp  = []; end % permutation group is empty by default

[n,ncY] = size(Y);   % # sites, # species
ncX     = size(X,2); % # predictors
if n   ~= size(X,1), error('Y & X need same # of rows'); end;

% Set W = [] for no covariables:
if isequal(0,W), W = []; end

% Check for covariables:
if isempty(W)
   partial = 0; % no partial analysis
   ncW     = 0; % # covariables
else
   if n ~= size(W,1), error('Y & W need same # of rows!');    end;
   if isequal(X,W),   error('X & W cannot be the same matrix!'); end
   
   partial = 1;         % do partial analysis
   ncW     = size(W,2); % # covariables
end

% Check permutation method:
if iter<0
   iter    = abs(iter);
   permRes = 0; % permute raw data if # iterations is negative
else
   permRes = 1; % permute residuals
end

% Check options during PERMUTATION RUN:
if (perm>0)
   if (verb>0) || (simu>0)
      error('When PERM=1, VERB & SIMU should be 0!')
   end
end

% Check options during SIMULATION RUN:
if (simu>0)
   if (verb>0) || (perm>0)
      error('When SIMU=1, set VERB & PERM to 0!')
   end
end
% -------------------------------------

if (perm<1) % skip during a permutation run: [***]
   Y = f_center(Y);     % center response variables           (L&L,1998: p.580)
   if (stnd>0)
      X = f_stnd(X);    % standardize explanatory variables
      if ~isempty(W)
         W = f_stnd(W); % standardize covariables
      end
   end
end

% -----Setup groups for restricted permutations:-----
if ~isempty(grp)
   uGrp = unique(grp,'stable'); % unique groups
   nGrp = numel(uGrp);          % # of groups
   
   % Create logical index to each group:
   idxGrp = false(n,nGrp);      % preallocate
   for i = 1:nGrp
      idxGrp(:,i) = (grp == uGrp(i)); % get logical index to this group
   end
   % Get size of each group:
   sGrp = sum(idxGrp);
end
% ---------------------------------------------------

% =========================================================================
if (partial<1) % No covariables:
   B            = X\Y;      % regression coefficients via QR
   Yfit         = X*B;      % fitted values of Y
   Yres         = Y - Yfit; % residuals
   
   result.B     = B;        % Output to results
   
else
   % Variables + covariables:
   XW           = [X W];
   B_XW         = XW\Y;
   Yfit_XW      = XW*B_XW;
   Yres         = Y - Yfit_XW;
   
   % Covariables:
   B_W          = W\Y;
   Yfit_W       = W*B_W;
   Yres_W       = Y - Yfit_W;       % needed for permutation test
   Yfit         = Yfit_XW - Yfit_W; % partial out covariables
   
   result.B_XW  = B_XW;     % For interaction term of X and W
   result.B_W   = B_W;      % For just W
end
% =========================================================================

% McArdle & Anderson (2001):
SSr = trace(Yfit'*Yfit); % Sum-of-Squares of regression model
SSe = trace(Yres'*Yres); % Sum-of-Squares error

if (partial>0)
   SSr_W = trace(Yfit_W'*Yfit_W); % Sum-of-Squares of covariable
end

% L&L,1998 (eq. 11.19); L,2007 (eq. 3):
dfR = ncX;         % Degrees-of-Freedom regression model
dfE = n-ncX-ncW-1; % Degrees-of-Freedom error
MSr = SSr/dfR;     % Mean Squares regression model
MSe = SSe/dfE;     % Mean Squares error
F   = MSr/MSe;     % Observed F-stat

% Coefficient of determination:
SSt   = trace(Y'*Y);                    % Sum-of-Squares total
R2    = SSr/SSt;                        % Legendre,2007 (eq. 5)
R2adj = 1 - ((1-R2)*((n-1)/(n-ncX-1))); % Legendre,2007 (eq. 6)

% Variability explained by covariate:
if (partial>0)
   R2w = SSr_W/SSt;
end

% Stop here for PERMUTATION RUNS:
if (perm>0)
   result.F     = F;
   result.p     = NaN; % <- no p-value
   result.R2    = R2;
   result.R2adj = R2adj;
   result.SSe   = SSe;
   result.SSt   = SSt;
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PERMUTATION TEST:                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (iter>0)
   Fperm = [F;zeros(iter-1,1)]; % initialize
   
   if (partial<1) % NO COVARIABLES:
      % Permutation of residuals under a full model (L&L,1998 p.608):
      if (simu<1) % Show messages if not a simulation run
         if (permRes>0)
            fprintf('\nPermuting residuals under a full model %d times...\n',iter-1);
         else
            fprintf('\nPermuting raw data under a full model %d times...\n',iter-1);
         end
      end
      
      for i = 2:iter % observed value is considered a permutation
         idxS = randperm(n)'; % get shuffled index to rows
         if (permRes>0)
            temp = f_rda(Yres(idxS,:),X,[],0,0,0,1,0); % permute rows of residuals
         else
            temp = f_rda(Y(idxS,:)   ,X,[],0,0,0,1,0); % permute rows of raw data
         end
         Fperm(i) = temp.F;
      end
      
   else % WITH COVARIABLES:
      % Permutation of residuals under a reduced model (L&L,1998 p.609-610):
      
      if (simu<1) % Show messages if not a simulation run
         fprintf('\nPermuting residuals under a reduced model %d times...\n',iter-1);
      end
      
      for i = 2:iter % observed value is considered a permutation
         if isempty(grp)
            idxS = randperm(n)'; % get shuffled index to rows
         else
            % Permutations restricted within groups:
            idxS = NaN(n,1); % preallocate
            for j = 1:nGrp
               idxS(idxGrp(:,j)) = randperm(sGrp(j))'; % permute this group
            end
         end
         % Add permuted residuals to form new Y:
         temp     = f_rda(Yfit_W + Yres_W(idxS,:),X,W,0,0,0,1,0);
         Fperm(i) = temp.F;
      end
   end
   p = sum(Fperm>=F)/iter; % convert counts to a p-value
   
   % Optionally plot histogram of permuted F-ratios:
   if (verb>1)
      figure;
      [bin,xBin] = hist(Fperm,100);  % get 100 bins
      bin  = bin/iter;               % convert absolute to relative
      h(1) = bar(xBin,bin);
      title('Frequency Histogram of Permuted F-ratios')
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
   result.F     = F;
   result.p     = p; % <- return p-value
   result.R2    = R2;
   result.R2adj = R2adj;
   result.SSe   = SSe;
   result.SSt   = SSt;
   return
end

if ncY==1 % Univariate response data:
   siteScores = Y;
   fitScores  = Yfit;
   resScores  = Yres;
   U          = NaN;
   Ures       = NaN;
   evals      = NaN;
   evalsRes   = NaN;
   canVar     = 1;
   canonCorr  = NaN;
   
else % Multivariate response data:
   % PCA:
   pca      = f_pca(Yfit); % fitted Y
   pcaRes   = f_pca(Yres); % residuals
   U        = pca.evects;
   evals    = pca.evals(:)' *(n-1);
   Ures     = pcaRes.evects;
   evalsRes = pcaRes.evals(:)' * (n-1);
   clear pca pcaRes;
   
   % Project in canonical space:
   siteScores = Y*U;       % site scores        (in space of Y)         (eq.11.12)
   fitScores  = Yfit*U;    % fitted site scores (in space of X)         (eq.11.13)
   resScores  = Yres*Ures; % residual scores    (in space of residuals) (fig.11.2)
   
   % Trim unnecessary axes:
   s          = min([ncX ncY (n-1)]); % # non-zero canonical eigenvalues
   sRes       = min([ncY (n-1)]);     % # non-zero residual eigenvalues
   siteScores = siteScores(:,1:s);
   fitScores  = fitScores(:,1:s);
   resScores  = resScores(:,1:sRes);
   U          = U(:,1:s);
   Ures       = Ures(:,1:sRes);
   evals      = evals(1:s);
   evalsRes   = evalsRes(1:sRes);
   
   % Only keep residual axes with variances larger than 0 (L&L,1998:p.590)
   idx              = abs(evalsRes) < 0.000001;
   evalsRes(idx)    = [];
   Ures(:,idx)      = [];
   resScores(:,idx) = [];
   
   % Species-Environment correlation (r) of each exis:
   canonCorr = zeros(1,s); % preallocate
   for i = 1:s
      canonCorr(i) = f_corr(siteScores(:,i),fitScores(:,i));
   end
   
   % Proportion of variance:
   totVar = sum([evals evalsRes]); % total variance in Y (= total inertia)
   canVar = evals./totVar;         % fraction of variance explained by canonical axes
   resVar = evalsRes./totVar;      % fraction of variance explained by residual axes
end


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   if (partial<1)
      fprintf('REDUNDANCY ANALYSIS:\n');
   else
      fprintf('Partial REDUNDANCY ANALYSIS:\n');
   end
   fprintf('--------------------------------------------------\n');
   fprintf(' F  = %-3.4f    p    =  %-3.4f (iter=%d) \n',F,p,iter);
   fprintf('R2  = %-3.4f   R2adj =  %-3.4f \n',R2,R2adj);
   
   if (partial>0)
      fprintf('R2w = %-3.4f (covariate) \n',R2w);
   end
   
   fprintf('--------------------------------------------------\n');
   
   if (ncY==1) % Univariate response variable:
      fprintf('Response variable (Y) is univariate. \n');
      if (stnd>0)
         if isempty(W)
            fprintf('\n(X variables standardized)\n');
         else
            fprintf('\n(X & W variables standardized)\n');
         end
         if (~isempty(grp) && (iter>0))
            fprintf('\n(Restricted permutations used)\n');
         end
      end
      fprintf('--------------------------------------------------\n');
   else        % Multivariate response variable:
      
      fprintf('\nCanonical Eigenvalues:\n');
      fprintf('  %-3.4f',evals);
      fprintf('\nResidual Eigenvalues:\n');
      fprintf('  %-3.4f',evalsRes);
      
      fprintf('\n\nSpecies-Environment Correlations (r):\n');
      fprintf('  %-3.4f',canonCorr);
      
      if (partial<1)
         fprintf('\n\nFraction of variance explained:\n');
      else
         fprintf('\n\nFraction of RESIDUAL variance explained (covariate removed):\n');
      end
      fprintf('------------------------------\n');
      fprintf('Canonical axes (total = %-3.4f): \n',sum(canVar));
      fprintf('  %-3.4f',canVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %-3.4f',cumsum(canVar));
      fprintf('\n\nResidual axes  (total = %-3.4f):\n',sum(resVar));
      fprintf('  %-3.4f',resVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %-3.4f',cumsum(resVar));
      fprintf('\n------------------------------\n');
      if (stnd>0)
         if isempty(W)
            fprintf('\n(X variables standardized)\n');
         else
            fprintf('\n(X & W variables standardized)\n');
         end
      end
      if (~isempty(grp) && (iter>0))
         fprintf('\n(Restricted permutations used)\n');
      end
      fprintf('==================================================\n');
   end
end
% ---------------------------------


% -----Wrap results up into a structure:-----
result.F          = F;
result.p          = p;
result.R2         = R2;
result.R2adj      = R2adj;
if (partial>0)
   result.R2w     = R2w;
end

result.dfR        = dfR;
result.dfE        = dfE;
result.dfT        = n-1;

result.SSr        = SSr;
if (partial>0)
   result.SSr_W   = SSr_W;
end
result.SSe        = SSe;
result.SSt        = SSt;

result.MSr        = MSr;
result.MSe        = MSe;

result.X          = X;
result.Y          = Y;
result.fit        = Yfit;
result.res        = Yres;

result.siteScores = siteScores;
result.fitScores  = fitScores;
result.resScores  = resScores;

result.U          = U;
result.Ures       = Ures;
result.evals      = evals;
result.evalsRes   = evalsRes;

result.canVar     = canVar;
result.canCorr    = canonCorr;
% ------------------------------------------
