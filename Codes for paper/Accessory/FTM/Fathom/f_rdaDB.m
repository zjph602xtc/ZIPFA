function result = f_rdaDB(yDis,ncY,X,W,iter,verb,stnd,skip,simu,grp)
% - distance-based Redundancy Analysis (db-RDA)
%
% USAGE: result = f_rdaDB(yDis,ncY,X,W,iter,verb,stnd,skip);
%
% yDis  = square symmetric distance matrix derived from response variables
% ncY   = # columns of original transformed data used to derive yDis
%
% X     = (1) matrix of explanatory variables, or
%         (2) ANOVA design matrix specified by dummy coding
%
% W     = matrix of covariables                                   (default = [])
% iter  = # iterations for permutation test                        (default = 0)
%         + value: permute residual
%         - value: permute raw data
% verb  = 1: send results to display; 2: also create plot          (default = 1)
% stnd  = optionally standardize X & W                             (default = 1)
% skip  = 0: proceed normally                                      (default = 0)
%         1: permutation run (internal flag)
%         2: input is a symmetric residual matrix
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
%  .Y          = NaN
%  .G          = Gower's centered version of yDis (= variance-covariance matrix)
%  .Gfit       = fitted values of G
%  .Gres       = residual values of G
%  .siteScores = site scores        (= WA scores)
%  .fitScores  = fitted site scores (= LC scores)
%  .resScores  = residual scores
%  .U          = canonical eigenvectors
%  .Ures       = residual eigenvectors
%  .evals      = canonical eigenvalues
%  .evalsRes   = residual eigenvalues
%  .canVar     = fraction of variance explained by canonical axes (= R^2)
%  .canCorr    = Species-Environ correlations (r) of each canonical axis
%
% SEE ALSO: f_rdaDB_AIC,  f_rdaDB_Stepwise, f_rda, f_rdaPlot, f_mregress

% -----Notes:-----
% Distance-based Redundancy Analysis is a constrained form of PCoA. Therefore,
% it may be instructive to compare an ordination plot derived from PCoA, which
% depicts a maximum of the total variation in the response variable with a
% minimum number of axes, vs. db-RDA, which depicts a maximum of the total
% variation in the repsonse variable explained by the predictor variable(s) with
% a minimum number of axes.
%
% ITER = specify a POS value to permute residuals
% ITER = specify a NEG value to permute raw data
%
% STND = 0: If X or W are a design matrix of dummy codes set STND=1 and do not
%           include an intercept term, otherwise the F-ratio will be incorrect.
%
% SKIP = 0 : proceed normally (default)
% SKIP = 1 : permutation run (used as an internal flag for permutation tests)
%             (a) skip checking input is a symmetric distance matrix
%             (b) skip centering of input by Gower's method
%             (c) skip standardizing X & W
%             (d) only return F, R2, R2adj, SSe, SSt
% SKIP = 2 : use when input matrix is a square symmetric residual matrix (= Gres)
%             (a) skip checking input is a symmetric distance matrix
%             (b) skip centering of input by Gower's method
%
% --PERMUTATION METHODS--:
% 1) Permutation of raw data under a full model     (X, no covariates)
% 2) Permutation of residual under a full model     (X, no covariates)
% 3) Permutation of residuals under a reduced model (W, only covariates)
% 4) Permutation of residuals under a full model    (X + W)


% -----RUN MODES:-----
% PERMUTATION RUN: (skip=2) returns only basic statistics so the calling
%                  function can calculate permutation based p-values
%
% SIMULATION RUN:  (simu=1) returns only basic statistics & p-values

% -----RESTRICTED PERMUTATION:-----
% grp: specifies the groups that permutations should be restricted within

% -----References:-----
% Anderson, M. J. 2001. A new method for non-parametric multivariate
%   analysis of variance. Austral Ecology 26: 32-46.
% Anderson, M. J. 2002. DISTML v.2: a FORTRAN computer program to calculate a
%   distance-based multivariate analysis for a linear model. Dept. of Statistics
%   University of Auckland. (http://www.stat.auckland.ac.nz/PEOPLE/marti/)
% Legendre, P. & L. Legendre. 1998. Numerical ecology. 2nd English ed.
%   Elsevier Science BV, Amsterdam.
% Legendre, P. 2007. Studying beta diversity: ecological variation partitioning
%   by multiple regression and canonical analysis. Journal of Plant Ecology,
%   doi:10.1093/jpe/rtm001
% McArdle, B. H. and M. J. Anderson. 2001. Fitting multivariate models to
%   community data: a comment on distance-based redundancy analysis. Ecology
%   290-297.
% Schielzeth, H. 2010. Simple means to improve the interpretability of
%   regression coefficients. Methods in Ecology & Evolution 1(2): 103-113.
%
%  Most of the equations related to db-RDA are based on McArdle & Anderson (2001).

% -----Author:-----
% by David L. Jones, Mar-2008
%
% This file is part of the FATHOM Toolbox for Matlab and
% is released under the GNU General Public License, version 2.

% Mar-2010: sub_decompose no longer trims its results, eigenvectors are now
%           scaled
% Dec-2012: 'I' and 'uno' defined for permutation runs, too.
% Feb-2016: added support for permuting residuals or raw data; optional plot of
%           permuted F stats; internal flag 'perm' renamed 'skip'; users can now
%           provide a square symmetric residual (Gres) matrix as input; now
%           returns G, Gfit, and Gres; alternative formulas show equivalent
%           methods of calculating fitted values and residuals; now returns R2w;
%           updated permutation test; no longer calls f_shuffle; no long calls
%           f_eig from a subfunction; sub_decompose renamed f_eig; updated the
%           way residual axes with variance = 0 are trimmed; also trims residual
%           scores.
% Mar-2016: 'stnd' option now allows you to NOT standardize X & W
% Apr-2016: verpPerm replaced with SIMU; added 'simulation mode'; no covariate
%           now specified by W=[] instead of W=0
% May-2016: changed (ncY==-1) to (ncY==1)
% May-2017: now supports restricted permutations

% -----Set defaults & check input:-----
if (nargin < 4),  W    = []; end % no covariables by default
if (nargin < 5),  iter = 0;  end % no permutation test by default
if (nargin < 6),  verb = 1;  end % send output to display by default
if (nargin < 7),  stnd = 1;  end % default standardize X & W
if (nargin < 8),  skip = 0;  end % default proceed normally
if (nargin < 9),  simu = 0;  end % internal flag: default not a simulation run
if (nargin < 10), grp  = []; end % permutation group is empty by default

n    = size(yDis,1); % # sites
ncX  = size(X,2);    % # predictors
if n ~= size(X,1), error('yDis & X need same # of rows'); end

% Check input if this isn't a permutation run:
if (skip==0)
   if (f_issymdis(yDis) == 0)
      fprintf('\nUse SKIP = 2 if yDis is a residual matrix, otherwise...\n');
      error('Input yDIS must be a square symmetric distance matrix!');
   end
end

% Set W = [] for no covariables:
if isequal(0,W), W = []; end

% Check for covariables:
if isempty(W)
   partial = 0; % no partial analysis
   ncW     = 0; % # covariables
else
   if n ~= size(W,1), error('yDis & W need same # of rows!');    end;
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
if (skip==1)
   if (verb>0) || (simu>0)
      error('When SKIP=1, set VERB & SIMU to 0!')
   end
end

% Check options during SIMULATION RUN:
if (simu>0)
   if (verb>0) || (skip==1)
      error('When SIMU=1, set VERB==0 & SKIP=0 or 2!')
   end
end
% -------------------------------------

% Skip optional standardization if permutation run:
if (stnd>0) && (~skip==1)
   X = f_stnd(X); % standardize explanatory variables
   if ~isempty(W)
      W = f_stnd(W); % standardize covariables
   end
end

% Centering and Intercept terms:
I   = eye(n,n);
uno = ones(n,1);

if (skip==0)
   % Response variable:
   A = -0.5*(yDis.^2);
   G = (I-(1/n)*(uno*uno'))*A*(I-(1/n)*(uno*uno')); % Gower's centered matrix
else
   % yDis is already 'G' during a permutation run:
   G = yDis;
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
   [Q1,~]  = qr([uno X],0); H = Q1*Q1';      % Hat-matrix for X
   Gfit    = (H*G*H);                        % fitted values for X
   Gres    = (I-H)*G*(I-H);                  % residuals for X
   SSr     = trace(Gfit);                    % SS regression
   SSe     = trace(Gres);                    % SS error
   
else
   % Variables + covariables:
   [Q1,~]  = qr([uno X W],0); H_XW = Q1*Q1'; % Hat-matrix for XW
   Gfit_XW = (H_XW*G*H_XW);                  % fitted values for XW
   Gres    = (I-H_XW)*G*(I-H_XW);            % residuals for XW
   SSr_XW  = trace(Gfit_XW);                 % SS regression
   SSe     = trace(Gres);                    % SS error
   
   % Covariables:
   [Q1,~]  = qr([uno W],0); H_W = Q1*Q1';    % Hat-matrix for W
   Gfit_W  = (H_W*G*H_W);                    % fitted values for W
   Gres_W  = (I-H_W)*G*(I-H_W);              % residuals for W
   SSr_W   = trace(Gfit_W);                  % SS regression for W
   SSr     = SSr_XW - SSr_W;                 % partial out covariables
   Gfit    = (H_XW - H_W)*G*(H_XW - H_W);    % fitted value for X
end

% Alternative formulas produce different versions of
% the 'G' matrices, but the trace values are equivalent:
%
% --No covariables:--
% Gres   = G - Gfit;
% Gres   = (I-H)*G*(I-H);
%
% --Covariables:--
% Gres   = G - Gfit_XW
% Gres   = (I-H_XW)*G*(I-H_XW);
% Gres_W = G - Gfit_W;
% Gres_W = (I-H_W)*G*(I-H_W);
% Gfit   = Gfit_XW - Gfit_W;
% Gfit   = (H_XW - H_W)*G*(H_XW - H_W);
% =========================================================================

% L&L,1998 (eq. 11.19); L,2007 (eq. 3)
dfR = ncX;         % Degrees-of-Freedom regression model
dfE = n-ncX-ncW-1; % Degrees-of-Freedom error
MSr = SSr/dfR;     % Mean Squares regression model
MSe = SSe/dfE;     % Mean Squares error
F   = MSr/MSe;     % Observed F-stat; same as F = (SSr/ncX)/(SSe/(n-ncX-ncW-1))

% Coefficient of determination:
SSt   = trace(G);                       % SS total
R2    = SSr/SSt;                        % Legendre,2007 (eq. 5)
R2adj = 1 - ((1-R2)*((n-1)/(n-ncX-1))); % Legendre,2007 (eq. 6)

% Variability explained by covariate:
if (partial>0)
   R2w = SSr_W/SSt;
end

% Stop here for PERMUTATION RUNS:
if (skip==1)
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
         idxS = randperm(n)'; % get shuffled index to rows/cols
         if (permRes>0)
            temp = f_rdaDB(Gres(idxS,idxS),ncY,X,[],0,0,0,1,0); % permute residuals
         else
            temp = f_rdaDB(G(idxS,idxS),   ncY,X,[],0,0,0,1,0); % permute raw data
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
            idxS = randperm(n)'; % get shuffled index to rows/cols
         else
            % Permutations restricted within groups:
            idxS = NaN(n,1); % preallocate
            for j = 1:nGrp
               idxS(idxGrp(:,j)) = randperm(sGrp(j))'; % permute this group
            end
         end
         % Add permuted residuals to form new G:
         temp     = f_rdaDB(Gfit_W + Gres_W(idxS,idxS),ncY,X,W,0,0,0,1,0);
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

% Eigenvalue decomposition:
[U,evals]       = sub_eig(Gfit); % fitted Y
[Ures,evalsRes] = sub_eig(Gres); % residuals

% % ***CONVERT EIGENVALUES [DLJ]:*** <- don't do this
% evals    = evals/(n-1);
% evalsRes = evalsRes/(n-1);

% Project PCoA's in canonical space:
siteScores = G*U;       % site scores        (in space of Y) (eq.11.12)
fitScores  = Gfit*U;    % fitted site scores (in space of X) (eq.11.13)
resScores  = Gres*Ures; % residual scores    (in space of residuals) (fig.11.2)

% Scale axes to sqrt of their eigenvalue:
siteScores = siteScores .* repmat(abs((evals.^0.5)),size(siteScores,1),1);
fitScores  = fitScores  .* repmat(abs((evals.^0.5)),size(fitScores,1),1);
resScores  = resScores  .* repmat(abs((evalsRes.^0.5)),size(resScores,1),1);

% Trim unnecessary axes:
s          = min([ncX ncY (n-1)]); % # non-zero canonical eigenvalues
% sRes       = min([ncY (n-1)]);     % # non-zero residual eigenvalues
siteScores = siteScores(:,1:s);
fitScores  = fitScores(:,1:s);
% resScores  = resScores(:,1:sRes);
U          = U(:,1:s);
% Ures       = Ures(:,1:sRes);
evals      = evals(1:s);
% evalsRes   = evalsRes(1:sRes);

% Only keep residual axes with variances larger than 0 (L&L,1998:p.590)
% idx              = abs(evalsRes) < 0.000001;
idx              = abs(evalsRes) < sqrt(eps);
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


% -----Send output to display:-----
if (verb>0)
   fprintf('\n==================================================\n');
   if (partial<1)
      fprintf('db-REDUNDANCY ANALYSIS:\n');
   else
      fprintf('Partial db-REDUNDANCY ANALYSIS:\n');
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
         if (isempty(W))
            fprintf('\n(X variables standardized)\n');
         else
            fprintf('\n(X & W variables standardized)\n');
         end
      end
      if (~isempty(grp) && (iter>0))
         fprintf('\n(Restricted permutations used)\n');
      end
      fprintf('--------------------------------------------------\n');
   else        % Multivariate response variable:
      
      fprintf('\nCanonical Eigenvalues:\n');
      fprintf('  %3.4f',evals);
      fprintf('\nResidual Eigenvalues:\n');
      fprintf('  %3.4f',evalsRes);
      
      fprintf('\n\nSpecies-Environment Correlations (r):\n');
      fprintf('  %3.4f',canonCorr);
      
      if (partial<1)
         fprintf('\n\nFraction of variance explained:\n');
      else
         fprintf('\n\nFraction of RESIDUAL variance explained (covariate removed):\n');
      end
      fprintf('------------------------------\n');
      fprintf('Canonical axes (total = %-3.4f): \n',sum(canVar));
      fprintf('  %3.4f',canVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %3.4f',cumsum(canVar));
      fprintf('\n\nResidual axes  (total = %-3.4f):\n',sum(resVar));
      fprintf('  %3.4f',resVar);
      fprintf('\n');
      fprintf('Cumulative: \n');
      fprintf('  %3.4f',cumsum(resVar));
      fprintf('\n------------------------------\n');
      if (stnd>0)
         if (isempty(W))
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
%
result.SSr        = SSr;
if (partial>0)
   result.SSr_W   = SSr_W;
end
result.SSe        = SSe;
result.SSt        = SSt;

result.MSr        = MSr;
result.MSe        = MSe;

result.X          = X;
result.Y          = NaN;
result.G          = G;
result.Gfit       = Gfit;
result.Gres       = Gres;

result.siteScores = siteScores;
result.fitScores  = fitScores;
result.resScores  = resScores;

result.U          = U;
result.Ures       = Ures;
result.evals      = evals;
result.evalsRes   = evalsRes;

result.canVar     = canVar;
result.canCorr    = canonCorr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                SUBFUNCTION:                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [evects,evals] = sub_eig(M)
% - eigenvalue decomposition (after f_eig)

% Eigenanalysis:
[U,D] = eig(M);
evals = diag(D);

% Discard imaginary components:
U     = real(U);
evals = real(evals);

% Sort by eigenvalues, descending:
sorted = flipud(sortrows([U' evals],size(U',1)+1));
evects = sorted(:,(1:end-1))';
evals  = sorted(:,end)'; % row vector of eigenvalues (1 for each axis)
