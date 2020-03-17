function result = f_vpa(Y,X,W,Z,iter,verb,plt)
% - variation partitioning analysis
%
% USAGE: result = f_vpa(Y,X,W,Z,iter,verb,plt)
% Y    = matrix of response variables (rows = sites, cols = species)
% X    = 1st matrix of explanatory variables
% W    = 2nd matrix of explanatory variables
% Z    = 3rd matrix of explanatory variables                      (default = [])
% iter = # iterations for permutation test                         (default = 0)
% verb = optionally send results to display                        (default = 1)
% plt  = plot Venn diagram                                         (default = 0)
%
% result = structure of results with the following fields:
%  .A   = A;
%  .B   = B;
%  .C   = C;
%  .ABC = ABC;
%  .D   = D;
%
% SEE ALSO: f_vpaDB, f_rda

% -----References:-----
% Legendre, P., and L. Legendre. 2012. Numerical ecology, 3rd English edition.
%   Developments in Environmental Modelling, Vol. 24. Elsevier Science BV,
%   Amsterdam. xiv + 990 pp. [See page 573]
%
% 'varpart3.R' from the 'Vegan package for R' was consulted for help with
% partitioning three explanatory matrices.

% -----Author:-----
% by David L. Jones, Sep-2016
%
% This file is part of the 'FATHOM TOOLBOX FOR MATLAB' and
% is released under the GNU General Public License, version 2.

% Apr-2017: now allows up to 3 explanatory matrices

% -----Set defaults & check input:-----
if (nargin < 4), Z    = []; end % default only 2 explanatory matrices
if (nargin < 5), iter = 0;  end % no permutation test by default
if (nargin < 6), verb = []; end % send output to display by default
if (nargin < 7), plt  = 0;  end % default don't create Venn diagram

% Set Z=0 to Z=[]:
if isequal(0,Z), Z = []; end

% Define # explanatory matrices:
if isempty(Z)
   nEXP = 2;
else
   nEXP = 3;
end

% Check size of input:
n = size(Y,1); % get # rows
if (n ~= size(X,1)) || (n ~= size(W,1))
   error('Y, X, and W must have same # rows!')
end
if (nEXP==3) && (n ~= size(Z,1))
   error('Y and Z must have same # rows!')
end

% Set options for RDA:
nIter = 0; % RDA does not perform permutation test
vb    = 0; % RDA does not display output
stnd  = 0; % RDA does not standardize X or W
perm  = 1; % RDA returns only basic statistics
simu  = 0; % RDA returns only basic statistics & p-value
% -------------------------------------

% Center response variables (L&L,1998: p.580)
Y = f_center(Y);

switch nEXP
   case 2
      % Standardize explanatory variables:
      X = f_stnd(X);
      W = f_stnd(W);
      
      % Collect adjusted R-squared values:
      temp = f_rda(Y,[X W],0,nIter,vb,stnd,perm,simu);
      ABC  = temp.R2adj;
      
      temp = f_rda(Y,X,0,nIter,vb,stnd,perm,simu);
      AB   = temp.R2adj;
      
      temp = f_rda(Y,W,0,nIter,vb,stnd,perm,simu);
      BC   = temp.R2adj;
      
      % Partition the variation:
      D = 1   - ABC;
      A = ABC - BC;
      C = ABC - AB;
      B = AB  - A;
      
      % Perform Partial RDA to calculate p-values:
      temp   = f_rda(Y,X,W,iter,vb,stnd,0,1);
      p.A    = temp.p;
      
      temp   = f_rda(Y,W,X,iter,vb,stnd,0,1);
      p.C    = temp.p;
            
      % -----Display table of results:-----
      if (verb>0)
         % Create a table:
         varName      = {'Fraction','R2adj' 'p'}; % define variable names
         rowName      = {'A' 'B' 'C' 'D'}';       % define row names
         tbl_R2adj    = [A B C D]';               % tabulate R2adj values
         tbl_p        = NaN(numel(tbl_R2adj),1);  % tabulate p-values
         tbl_p([1 3]) = [p.A p.C]';
         T            = table(rowName,tbl_R2adj,tbl_p,'VariableNames', varName);
         fprintf('\n==================================================\n');
         fprintf('      Variation Partitioning Analysis (VPA):        \n');
         fprintf('--------------------------------------------------\n');
         fprintf('\n')
         disp(T) % display table
         fprintf('\n--------------------------------------------------\n');
         fprintf('Fractions for each explanatory matrix:\n')
         fprintf('   X = %-3.3f  \n',AB);
         fprintf('   W = %-3.3f  \n',BC);
         fprintf('--------------------------------------------------\n');
      end
      % -----------------------------------
      
      % Create Venn Diagram:
      if (plt>0)
         figure;
         f_venn([AB BC],B,'FaceColor',{[1 1 1]*0.5,[1 1 1]*0.5},...
            'FaceAlpha',{0.6,0.6},'EdgeColor','black', 'LineWidth',{1.5,1.5},...
            'LineStyle',{'-','--'});
         box on;
         % set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
         axis equal;
         axis(axis*1.1)
         
         % Add square:
         axisVar = axis;
         cen(1) = axisVar(1) + (axisVar(2) - axisVar(1))/2;
         cen(2) = axisVar(3) + (axisVar(4) - axisVar(3))/2;
         rectangle('Position', [cen - 0.5 1 1])
         axis tight;
         axis(axis*1.5)
      end
      
      % Wrap results up into a structure:
      result.A   = A;
      result.B   = B;
      result.C   = C;
      result.ABC = ABC;
      result.D   = D;
      
   case 3
      % Standardize explanatory variables:
      X = f_stnd(X);
      W = f_stnd(W);
      Z = f_stnd(Z);
      
      % Collect adjusted R-squared values:
      temp   = f_rda(Y,X,0,nIter,vb,stnd,perm,simu);
      ADFG   = temp.R2adj;
      
      temp   = f_rda(Y,W,0,nIter,vb,stnd,perm,simu);
      BDEG   = temp.R2adj;
      
      temp   = f_rda(Y,Z,0,nIter,vb,stnd,perm,simu);
      CEFG   = temp.R2adj;
      
      temp   = f_rda(Y,[X W],0,nIter,vb,stnd,perm,simu);
      ABDEFG = temp.R2adj;
      
      temp   = f_rda(Y,[X Z],0,nIter,vb,stnd,perm,simu);
      ACDEFG = temp.R2adj;
      
      temp   = f_rda(Y,[W Z],0,nIter,vb,stnd,perm,simu);
      BCDEFG = temp.R2adj;
      
      temp   = f_rda(Y,[X W Z],0,nIter,vb,stnd,perm,simu);
      ABCDEFG = temp.R2adj;
      
      % Partition the variation:
      H = 1       - ABCDEFG;
      A = ABCDEFG - BCDEFG;
      B = ABCDEFG - ACDEFG;
      C = ABCDEFG - ABDEFG;
      
      D = ACDEFG  - CEFG   - A;
      E = ABDEFG  - ADFG   - B;
      F = BCDEFG  - BDEG   - C;
      G = ADFG    - A      - D  - F;
      
      % Perform Partial RDA to calculate p-values:
      temp   = f_rda(Y,X,[W Z],iter,vb,stnd,0,1);
      p.A    = temp.p;
      
      temp   = f_rda(Y,W,[X Z],iter,vb,stnd,0,1);
      p.B    = temp.p;
      
      temp   = f_rda(Y,Z,[X W],iter,vb,stnd,0,1);
      p.C    = temp.p;
            
      % -----Display table of results:-----
      if (verb>0)
         % Create a table:
         varName    = {'Fraction','R2adj' 'p'};           % define variable names
         rowName    = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'}'; % define row names
         tbl_R2adj  = [A B C D E F G H]';                 % tabulate R2adj values
         tbl_p      = NaN(numel(tbl_R2adj),1);            % tabulate p-values
         tbl_p(1:3) = [p.A p.B p.C]';
         T          = table(rowName,tbl_R2adj,tbl_p,'VariableNames', varName);
         fprintf('\n==================================================\n');
         fprintf('      Variation Partitioning Analysis (VPA):        \n');
         fprintf('--------------------------------------------------\n');
         fprintf('\n')
         disp(T) % display table
         fprintf('\n--------------------------------------------------\n');
         fprintf('Fractions for each explanatory matrix:\n')
         fprintf('   X = %-3.3f  \n',ADFG);
         fprintf('   W = %-3.3f  \n',BDEG);
         fprintf('   Z = %-3.3f  \n',CEFG);
         fprintf('--------------------------------------------------\n');
      end
      % -----------------------------------
      
      % Create Venn Diagram:
      if (plt>0)
         figure;
         f_venn([ADFG BDEG CEFG],[D+G F+G E+G G],'FaceColor',{[1 1 1]*0.5,[1 1 1]*0.5,[1 1 1]*0.5},...
            'FaceAlpha',{0.6,0.6,0.6},'EdgeColor','black', 'LineWidth',{1.5,1.5,1.5},...
            'LineStyle',{'-','--',':'},'ErrMinMode','TotalError');
         box on;
         % set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',[],'YTick',[])
         axis equal;
         axis(axis*1.1)
         
         % Add square:
         axisVar = axis;
         cen(1) = axisVar(1) + (axisVar(2) - axisVar(1))/2;
         cen(2) = axisVar(3) + (axisVar(4) - axisVar(3))/2;
         rectangle('Position', [cen - 0.5 1 1])
         axis tight;
         axis(axis*1.5)
      end
      
      % Wrap results up into a structure:
      result.A = A;
      result.B = B;
      result.C = C;
      result.D = D;
      result.E = E;
      result.F = F;
      result.G = G;
      result.H = H;
end
