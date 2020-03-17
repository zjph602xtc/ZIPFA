function result = f_rdaAxes(Y,X,W,iter,verb,stnd,res)
% - test significance of canonical axes for RDA
%
% USAGE: result = f_rdaAxes(Y,X,W,iter,verb,stnd,res)
%
% Y    = matrix of response variables (rows = sites, cols = species)
%
% X    = (1) matrix of explanatory variables, or
%        (2) ANOVA design matrix specified by dummy coding
%
% W    = matrix of covariables                                        (0 = none)
% iter = # iterations for permutation test                         (default = 0)
% verb = send results to display                                   (default = 1)
% stnd = optionally standardize X & W                              (default = 1)
% res  = permute residuals (= 1) or create new residuals (= 0)     (default = 1)
%
% result = structure of results with the following fields:
%  .F     = F-statistics
%  .p     = permutation-based p-values
%  .evals = eigenvalues
%  .txt   = 'RDA + permuted residuals' or 'RDA + new residuals'
%
% SEE ALSO: f_rdaDB_Axes, f_rda

% -----Notes:-----
% This function follows Legendre et al. (2011) and uses the marginal test to
% determine the significance of each canonical axes from an RDA. This function
% has been tested against their 'test.axes.cov.R' and provides similar results.

% -----References:-----
% Legendre, P., J. Oksanen, and C. J. F. ter Braak. Testing the significance of
%  canonical axes in redundancy analysis. Methods in Ecology & Evolution 2:
%  269-277. (See 'test.axes.cov.R')

% -----Author:-----
% by David L. Jones, May-2017
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 3), W    = []; end % no covariables by default
if (nargin < 4), iter = 0;  end % no permutation test by default
if (nargin < 5), verb = 1;  end % send output to display by default
if (nargin < 6), stnd = 1;  end % default standardize X & W
if (nargin < 7), res  = 1;  end % default permute residuals

% Set W = [] for no covariables:
if isequal(0,W), W = []; end

% Get size of input:
n   = size(Y,1);

% Get rank of X and W (after 'test.axes.cov.R')
tol         = sqrt(eps);
[~,evalsX]  = f_svd(cov(X));
mX          = sum(evalsX>tol);
[~,evalsXW] = f_svd(cov([X W]));
mXW         = sum(evalsXW>tol);
clear tol evalsX evalsXW;

% Check for covariables:
if isempty(W)
   covar = 0; % no covariables
   a     = 2;
else
   covar = 1; % covariables
   a     = 1;
end
% -------------------------------------

% Turn warning off:
warning('off','MATLAB:rankDeficientMatrix');

% Center/standardize data:
Y = f_center(Y);     % center response variables           (L&L,1998: p.580)
if (stnd>0)
   X = f_stnd(X);    % standardize explanatory variables
   if ~isempty(W)
      W = f_stnd(W); % standardize covariables
   end
end

% Residualize X on W when covariables present:
if (covar>0)
   [~,X] = sub_rda(X,W,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 TEST AXIS I:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,SSe,evals,LC] = sub_rda(Y,X,0);          % RDA
nEvals             = numel(evals);            % get # eigenvalues
F                  = NaN(numel(evals),1);     % preallocate
p                  = NaN(numel(evals),1);
F(1)               = (evals(1)*(n-1-mX))/SSe; % eq 13 (Legendre et al. 2011)

% -----Permutation Test - Shuffle Raw Data:-----
if (covar==0) % Permute raw data when no covariables:
   if (iter>0)
      Fperm = [F(1);zeros(iter-1,1)]; % initialize
      for i = 2:iter
         idxS = randperm(n)'; % get shuffled index to rows
         [~,~,SSEperm,evalsPerm] = sub_rda(Y(idxS,:),X,0);
         Fperm(i)  = (evalsPerm(1)*(n-1-mX))/SSEperm;
      end
      p(1) = sum(Fperm>=F(1))/iter; % convert counts to a p-value
      
   else
      p(1) = NaN;
   end
end
% ----------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             TEST REMAINING AXES:                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get denominator of F-statistic:
[~,~,SSe_XW] = sub_rda(Y,[X W],1);

for i=a:nEvals
   % Get F-statistic for each eigenvalue:
   F(i) = (evals(i)*(n-1-mXW))/SSe_XW;
   
   % Residualize X on W + Previous Axes:
   if (i>1)
      [~,Xres] = sub_rda(X,[W LC(:,1:(i-1))],1);
   else
      Xres = X;
   end
   
   % -----Permutation Test:-----
   if (iter>0)
      Fperm = [F(i);zeros(iter-1,1)]; % initialize
      if (res>0) % --Permute Residuals:--
         [Yfit,Yres] = sub_rda(Y,[W LC(:,1:(i-1))],1);
         for j = 2:iter
            idxS              = randperm(n)';              % get shuffled index to rows
            Yperm             = Yfit + Yres(idxS,:);       % permute residuals
            SSTperm           = trace(Yperm'*Yperm);       % permuted Sum-of-Squares Total
            [~,~,~,evalsPerm] = sub_rda(Yperm,Xres,0);     % reduced model (might be 'rank deficient')
            YfitPerm          = sub_rda(Yperm,[X W],1);    % full model
            SSRperm           = trace(YfitPerm'*YfitPerm); % permuted Sum-of-Squares Regression
            Fperm(j)          = (evalsPerm(1)*(n-1-mX))/(SSTperm - SSRperm);
         end
      else % --Create New Residuals from Permuted Y:--
         for j = 2:iter
            idxS              = randperm(n)';                % get shuffled index to rows
            [~,Yres]          = sub_rda(Y(idxS,:),[W LC(:,1:(i-1))],1); % create new residuals
            SSTperm           = trace(Yres'*Yres);           % permuted Sum-of-Squares Total
            [~,~,~,evalsPerm] = sub_rda(Yres,Xres,0);        % reduced model (might be 'rank deficient')
            YfitPerm          = sub_rda(Yres,[X W],1);       % full model
            SSRperm           = trace(YfitPerm'*YfitPerm);   % permuted Sum-of-Squares Regression
            Fperm(j)          = (evalsPerm(1)*(n-1-mX))/(SSTperm - SSRperm);
         end
      end
      p(i) = sum(Fperm>=F(i))/iter; % convert counts to a p-value
   else
      p(i) = NaN;
   end
   % ---------------------------
end


% -----Display table of results:-----
if (verb>0)
   
   % Create a table:
   varName = {'Axis' 'Eigenvalue' 'F' 'p' }; % define variable names
   T = table((1:nEvals)', evals/(n-1), F, p, 'VariableNames', varName);
   
   % Display Table:
   % -----Send output to display:-----
   fprintf('\n==================================================\n');
   fprintf('          Test Significance of RDA Axes:          \n');
   fprintf('--------------------------------------------------\n');
   fprintf('\n')
   disp(T) % display table
   if (iter>0)
      fprintf('\n# iterations = %d\n',iter);
   end
   fprintf('--------------------------------------------------\n');
end
% -----------------------------------

% Wrap results up into a structure:
result.F     = F;
result.p     = p;
result.evals = evals/(n-1);   % convert eigenvalues for SSCP to COV matrix
if (res>0)
   result.txt = 'RDA + permuted residuals';
else
   result.txt = 'RDA + new residuals';
end

% Turn warning back on:
warning('on','MATLAB:rankDeficientMatrix');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 SUBFUNCTION:                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Yfit,Yres,SSe,evals,LC] = sub_rda(Y,X,skip)
% - mini-version of RDA after 'f_rda'

% Get size of input:
[n,ncY] = size(Y);
ncX     = size(X);

% Regression:
B    = X\Y;                 % regression coefficients via QR
Yfit = X*B;                 % fitted values of Y
Yres   = Y - Yfit;          % residuals
SSe    = trace(Yres'*Yres); % Sum-of-Squares error

% Optionally skip PCA:
if (skip>0)
   evals = [];
   LC    = [];
   return
else
   pca   = f_pca(Yfit);   % fitted Y
   U     = pca.evects;    % eigenvectors
   evals = pca.evals(:)'; % eigenvalues
   evals = evals'*(n-1);  % convert eigenvalues for COV matrix to SSCP
   LC    = Yfit*U;        % axes that are linear combinations of X (= LC)

   % Trim unnecessary axes:
   s     = min([ncX ncY (n-1)]); % # non-zero canonical eigenvalues
   LC    = LC(:,1:s);
   evals = evals(1:s);
end












