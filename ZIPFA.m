function [Ufit,Vfit,itr,finaltau,Likelihood,CVLikelihood] = ZIPFA(X, k, varargin)
% ZIPFA   To conduct the zero inflated zero factor analysis.
%	[Ufit,Vfit] = ZIPFA(X, k)
%
%	[Ufit,Vfit] = ZIPFA(X, k, tau, 'cut', 1, 'display', false, ...)
%
%	[Ufit,Vfit,itr,finaltau,Likelihood] = ZIPFA(X, k)
%
%	[Ufit,Vfit,itr,finaltau,Likelihood,CVLikelihood] = ZIPFA(X, k, 'missing', missingmat)
%
% -X: The matrix to be decomposed.
% -k: The number of factors.
% -tau (0.1): The initial guess for tau.
% -'cut' (0.8): Whether to delete columns that has more than 100('cut')% zeros. 'Cut' = 1, if no filtering.
% -'tolLnlikelihood' (5e-4): The max percentage of log likelihood differences in two iterations.
% -'iter' (20): Max iterations.
% -'tol' (1e-4): Percentage of l2 norm change of [tau beta] in ZIP regression.
% -'maxiter' (100): Max iterations in ZIP regression.
% -'initialtau' (iteration'): Choose the initial value of tau at the beginning of EM iteration in ZIP regression.
%       'stable': estimate tau from fitted beta in last round;
%       'initial': always use the initially assigned tau in 'tau' or 'initial';
%           Use the default tau = 0.1 if 'initial' is empty.
%       'iteration': use fitted tau in last round.
% -'Madj' (true): Whether adjust for relative library size M.
% -'display' (true): Display the fitting procedure.
% -'missing' ([]): T/F matrix. If 'missing' is not empty, then CVLikelihood is likelihood of X with missing = T.
% -'rept' ([]): Which fold is in cross validation. If rept is empty, then do not save result in this function.
%
%	See also cv_ZIPFA, EMzeropoisson_mat

% Tianchen Xu

[m,n0] = size(X);
p = inputParser;
addRequired(p,'k',@(x)validateattributes(x,{'numeric'},{'scalar','>=',1}));
addOptional(p,'tau',0.1,@isscalar);
addParameter(p,'cut',0.8,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}))
addParameter(p,'tolLnlikelihood',5e-4,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'iter',20,@(x)validateattributes(x,{'numeric'},{'scalar','>=',2}));
addParameter(p,'tol',1e-04,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'maxiter',100,@(x)validateattributes(x,{'numeric'},{'scalar','>=',2}));
addParameter(p,'initialtau','iteration',@(x)any(validatestring(x,{'stable','initial','iteration'})));
addParameter(p,'Madj',true,@islogical);
addParameter(p,'display',true,@islogical);
addParameter(p,'missing',[],@(x)validateattributes(x,{'logical'},{'size',[m,n0]}));
addParameter(p,'rept',[],@(x)isscalar(x)||isempty(x));

if length(dbstack)>1
    p.KeepUnmatched = true;
end
parse(p,k,varargin{:});

tau = p.Results.tau;
display = p.Results.display;
Madj = p.Results.Madj;
missing = p.Results.missing;
rept = p.Results.rept;

if (~isempty(missing))
    missing = missing(:,sum(X==0)/m<p.Results.cut);
end
X=X(:,sum(X==0)/m<p.Results.cut);
[m,n] = size(X);
if n~=n0
    warning('%d columns have been dropped due to the ''cut'' option.\n',n0-n)
end


Xt = X';
if any(isnan(X))
    if (~isempty(missing))
        error('Cannot assign ''missing'' when X already has missing values')
    end
    indexram = ~isnan(X(:));
    indexram_row = ~isnan(Xt(:));
else
    if (~isempty(missing))
        indexram = missing(:);
        indexram_row = missing';
        indexram_row = indexram_row(:);
    else
        indexram = 1:(m*n);
        indexram_row = indexram;
    end
end

XNA = X;
XNA(XNA==0) = nan;
lambdahat = mean(log(XNA),1,'omitnan');
for i=1:n
    XNA(isnan(XNA(:,i)),i) = round(exp(lambdahat(i)),0);
end

% initial value
[U,S,V] = svd(log(XNA));
Uold = U(:,1:k) * S(1:k,1:k);
Vold = V(:,1:k);

% iteration begins
Ufit = cell(p.Results.iter,1);
Vfit = cell(p.Results.iter,1);
Likelihood=[];
CVLikelihood=[];

if Madj
    ratio = repmat(sum(X,2,'omitnan')/median(sum(X,2,'omitnan')),n,1);
    ratio_row = reshape(ratio,m,[])';
    ratio_row = ratio_row(:);
end

for itr=1:p.Results.iter
    if display
        fprintf('\n ****************** \n Round %.0f \n ******************\n', itr)
    end
    %%% update U
    if Madj
        mi = ratio_row(indexram_row);
    else
        mi = 1;
    end
    Voldc = mat2cell(sparse(Vold),n,k);
    Voldc = repmat(Voldc, m, 1);
    dat = [Xt(:) blkdiag(Voldc{:})];
    dat = dat(indexram_row,:);
    if itr==1
        uuuu = EMzeropoisson_mat(dat,0, p.Results,'initial',[tau reshape(Uold',1,[])],'intercept',false,'m',mi);
    else
        uuuu = EMzeropoisson_mat(dat,0, p.Results,'initial',[uuuu(end,1) reshape(Ufit{itr-1}',1,[])],'intercept',false,'m',mi);
    end
    
    %  lnL_mat([Xt(:) blkdiag(Voldc{:})], uuuu, indexram_row, ratio_row, 'off')
    Unew = reshape(uuuu(end,2:end), k,[])';
    
    %%% update V
    if Madj
        mi = ratio(indexram);
    else
        mi = 1;
    end
    Uoldc = mat2cell(sparse(Uold),m,k);
    Uoldc = repmat(Uoldc, n, 1);
    dat = [X(:) blkdiag(Uoldc{:})];
    dat = dat(indexram,:);
    if itr==1
        uuuu = EMzeropoisson_mat(dat,0, p.Results,'initial',[tau reshape(Vold',1,[])],'intercept',false,'m',mi);
    else
        uuuu = EMzeropoisson_mat(dat,0, p.Results,'initial',[uuuu(end,1) reshape(Vfit{itr-1}',1,[])],'intercept',false,'m',mi);
    end
    Vnew = reshape(uuuu(end,2:end), k,[])';
    
    if Madj
        Likelihood=[Likelihood lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, indexram, ratio);];
        if ~isempty(missing) && nargout>5
            CVLikelihood=[CVLikelihood  lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, ~indexram, ratio);];
        end
    else
        Likelihood=[Likelihood lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, indexram, []);];
        if ~isempty(missing) && nargout>5
            CVLikelihood=[CVLikelihood  lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, ~indexram, []);];
        end
    end
    
    if display
        disp(Likelihood);
        disp(CVLikelihood);
    end
    
    %%% next step
    [U,S,V] = svd(Unew * Vnew');
    Uold = U(:,1:k) * S(1:k,1:k);
    Vold = V(:,1:k);
    
    %%% save answer
    Ufit{itr} = Uold;
    Vfit{itr} = Vold;
    
    if (itr~=1)
        if display
            fprintf('Max Ln likelihood diff = %.4g %% \n',full(100*(Likelihood(itr)-Likelihood(itr-1))/abs(Likelihood(itr-1))))
        end
        if (full((Likelihood(itr)-Likelihood(itr-1))/abs(Likelihood(itr-1))))<p.Results.tolLnlikelihood
            break
        end
    end
end

finaltau=uuuu(end,1);
Ufit = Ufit(~cellfun('isempty',Ufit));
Vfit = Vfit(~cellfun('isempty',Vfit));

if ~isempty(rept)
    save(['factor' num2str(k) 'fold' num2str(rept)])
end
end

function [likelihood] = lnL_mat(data, coefficient, subset, mi)
% note: data and mi are in full data size.
coefficient = coefficient(end,:);
tau = coefficient(1);
beta = coefficient(2:end)';
y = data(:,1);
n = length(y);

Dx = data(:,2:end);
lnlam = Dx * beta;
p = 1./ (1+exp(tau.*lnlam));
if isempty(mi)
    mi = ones(n,1);
end

like = y;
indexz = y==0;
indexnotz = y>=20;
indexnotzl = (y~=0) & (y<20);

like(indexz) = log(p(indexz)+(1-p(indexz)).*(poisspdf(0,mi(indexz).*exp(lnlam(indexz)))));
like(indexnotz) = log(1-p(indexnotz)) - mi(indexnotz).*exp(lnlam(indexnotz))...
    +y(indexnotz).*(log(mi(indexnotz))+lnlam(indexnotz))-(0.5.*log(2*pi*y(indexnotz))...
    +y(indexnotz).*(log(y(indexnotz))-1));
like(indexnotzl) = log(1-p(indexnotzl)) - mi(indexnotzl).*exp(lnlam(indexnotzl))...
    +y(indexnotzl).*(log(mi(indexnotzl))+lnlam(indexnotzl))-log(factorial(y(indexnotzl)));

if ~isempty(subset)
    like = like(subset);
end
likelihood = sum(like);
end

