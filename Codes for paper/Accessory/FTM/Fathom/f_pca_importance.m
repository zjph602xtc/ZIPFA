function result = f_pca_importance(pca,scale,verb,txt)
% - descriptor variable importance measures for a PCA
%
% USAGE: result = f_pca_importance(pca,scale,verb,txt)
%
% pca   = structure of results obtained from f_pca
% scale = method used to scale eigenvectors                        (default = 1)
%         1 : scaled to length equal to sqrt(eigenvalue) = 'PCA scaling 2'
%         2 : scaled by percent variation explained
% verb   = send results to display                                 (default = 1)
% txt    = cell array of descriptor labels; if empty, autocreate
%           e.g., txt = {'temp' 'sal' 'depth'};
%
% result = structure of results with the following fields
%  .Usc  = scaled eigenvectors
%  .expl = fraction of total variation occupied by each descriptor for SCALE=1
%  .len  = length of scaled descriptor eigenvector for SCALE=2
%  .txt  = descriptor labels
%
% SEE ALSO: f_pca

% -----NOTES:-----
% Usually the eigenvetors produced by a PCA are scaled to lengths = 1. However,
% if we scale them by the eigenvalues or the percent variation explained of each
% eigenvalue, the lengths of the scaled eigenvectors can be used to determine
% the relative importance of the original descriptor variables.

% -----References:-----
% Legendre, P. and E. D. Gallagher. 2001. Ecologically meaningful
%   transformations for ordination of species data. Oecologia 129: 271-280.
% Legendre, P., and L. Legendre. 2012. Numerical ecology, 3rd English edition.
%  Developments in Environmental Modelling, Vol. 24. Elsevier Science BV,
%  Amsterdam. xiv + 990 pp.
% 
% Leggallfigs.m by Eugene.Gallagher@umb.edu, 16-Mar-2002

% -----Author:-----
% by David L. Jones, Mar-2016
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% -----Set defaults & check input:-----
if (nargin < 2), scale = 1;  end % default scale eigenvector by sqrt(eigenvalue)
if (nargin < 3), verb  = 1;  end % send output to display by default
if (nargin < 4), txt   = []; end % default descriptor labels

% Autocreate txt labels:
if isempty(txt)
   txt = cellstr(num2str((1:size(pca.evects,2))'))';
end

% If labels are not cell arrays, try forcing them:
if iscell(txt)<1, txt = num2cell(txt); end

% Make sure labels are of compatible size:
txt = txt(:)';
if size(txt,2) ~= size(pca.evects,2)
   error('Size mismatch b/n TXT and PCA.EVECTS!')
end
% -------------------------------------

switch scale
   case 1
      % Extract variables:
      U = pca.evects;      % eigenvectors
      V = diag(pca.evals); % eigenvalues as a diagonal matrix
      
      % Scale eigenvectors to length equal to sqrt of eigenvalue, (i.e., their
      % standard deviation) (L&L:2012 eq.9.8):
      Usc = U*(V.^0.5);
      
      % Get fraction of variation explained for each descriptor:
      expl = Usc.^2/sum(Usc.^2); % Table 2 in L&G:2001 and Leggallfigs.m
      
      % Sort descriptors according to % variation explained:
      [~,key] = sortrows(expl,-1);
      
      % Display table of results:
      if (verb>0)
         varName = {'Frac_Explained' 'Descriptor' }; % define variable names
         T = table(expl(key),txt(key)','VariableNames',varName);
         fprintf('\n')
         disp(T) % display table
         fprintf('\n(NOTE: Eigenvectors scaled by sqrt(eigenvalue))\n\n')
      end
      
      % Wrap results up into a structure
      result.Usc  = Usc;
      result.expl = expl;
      result.txt  = txt;
      
   case 2
      % Extract variables:
      U       = pca.evects;    % eigenvectors
      PC_expl = pca.expl(:,1); % percent variation explained by each PC axis
      
      % Scale eigenvectors by % var explained:
      Usc = U .* repmat(PC_expl',size(U,1),1);
      
      % Get length of each scaled vector:
      len = sqrt(sum(Usc.^2,2)); % sqrt of sum of squared differences from 0
      
      % Sort descriptors according to length of scaled vectors:
      [~,key] = sortrows(len,-1);
            
      % Display table of results:
      if (verb>0)
         varName = {'Length' 'Descriptor' }; % define variable names
         T = table(len(key),txt(key)','VariableNames',varName);
         fprintf('\n')
         disp(T) % display table
         fprintf('\n(NOTE: Eigenvectors scaled by %% variation explained)\n\n')
      end
      
      % Wrap results up into a structure
      result.Usc = Usc;
      result.len = len;
      result.txt = txt;
      
   otherwise
      error('SCALE must be 1 or 2!')
end
