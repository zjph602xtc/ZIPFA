function [cvsample,Allres] = cv_ZIPFA (X, k, varargin)
% cv_ZIPFA    To conduct cross validation on ZIPFA model.
%   [cvsample,Allres] = cv_ZIPFA (X, k)
%
%   [cvsample,Allres] = cv_ZIPFA (X, k, fold, tau)
%
%   [cvsample,Allres] = cv_ZIPFA (X, k, 'cut', 1)
%
%   [cvsample,Allres] = cv_ZIPFA (X, k, 'display', false, 'savemat', true, ...)
%
% -X: The matrix to be decomposed.
% -k: The number of factors. It can be a vector.
% -fold (10): The number of folds used in cross validation.
% -tau (0.1): The initial guess for tau.
% -'cut' (0.8): Whether to delete columns that has more than 100('cut')% zeros. 'Cut' = 1, if no filtering.
% -'tolLnlikelihood' (5e-4): The max percentage of log likelihood differences in two iterations.
% -'iter' (20): Max iterations.
% -'tol' (1e-4): Percentage of l2 norm change of [tau beta] in ZIP regression.
% -'maxiter' (100): Max iterations in ZIP regression.
% -'initialtau' ('iteration'): Choose the initial value of tau at the beginning of EM iteration in ZIP regression.
%       'stable': estimate tau from fitted beta in last round;
%       'initial': always use the initially assigned tau in 'tau' or 'initial';
%           Use the default tau = 0.1 if 'initial' is empty.
%       'iteration': use fitted tau in last round.
% -'Madj' (true): Whether adjust for relative library size M.
% -'display' (true): Display the fitting procedure. Info in ZIPFA will not be shown in 'Parallel' mode even 'Display' is true.
% -'savemat' (false): Whether to save ZIPFA results in all factor numbers and each fold.
% -'parallel' (true): Use parallel toolbox to accelerate.
%
%	See also ZIPFA, EMzeropoisson_mat

% Tianchen Xu

p = inputParser;
addOptional(p,'fold',10,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
addOptional(p,'tau',0.1,@isscalar);
addParameter(p,'cut',0.8,@(x)validateattributes(x,{'numeric'},{'scalar','>=',0,'<=',1}))
addParameter(p,'tolLnlikelihood',5e-4,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'iter',20,@(x)validateattributes(x,{'numeric'},{'scalar','>=',2}));
addParameter(p,'tol',1e-04,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'maxiter',100,@(x)validateattributes(x,{'numeric'},{'scalar','>=',2}));
addParameter(p,'initialtau','iteration',@(x)any(validatestring(x,{'stable','initial','iteration'})));
addParameter(p,'Madj',true,@islogical);
addParameter(p,'display',true,@islogical);
addParameter(p,'parallel',true,@islogical);
addParameter(p,'savemat',false,@islogical);

if any(isnan(X))
    error('There is Na''s in X !')
end
validateattributes(k,{'numeric'},{'vector','>=',1});
parse(p,varargin{:});

[m, n0] =size(X);
X=X(:,sum(X==0)/m<p.Results.cut);
[m, n] =size(X);
if n~=n0
   warning('%d columns have been dropped due to the ''cut'' option.\n',n0-n)
end
fold = p.Results.fold;


% generate assignment
cvsample=[];
for i=1:fold-1
    cvsample=[cvsample i*ones(1,floor(m*n/fold))];
end
cvsample=[cvsample fold*ones(1,m*n-(fold-1)*floor(m*n/fold))];
cvsample=datasample(cvsample,length(cvsample),'Replace',false);
cvsample=reshape(cvsample,m,n);
% check existing assignment
try
    ext_cvsample=evalin('base','cvsample');
    validateattributes(ext_cvsample, {'numeric'},{'size',size(X)});
    disp('Use existing cvsample!')
    cvsample = ext_cvsample;
end
%
Allres = [];
for nf = k
    Mlike=zeros(fold,1);
    if p.Results.parallel
        parfor rept=1:fold
            missing=cvsample;
            missing(missing==rept)=0;
            missing(missing~=0)=1;
            if p.Results.savemat
                save = rept;
            else
                save = [];
            end
            [~,~,~,esttau,Likelihood,MLikelihood]=ZIPFA(X, nf, 0, p.Results,'missing',logical(missing),'rept',save,'display',false);
            if p.Results.display
                fprintf('\n factor = %d; fold = %d; esttau = %.3g; Likelihood = %.5g; MLikelihood = %.4g \n',nf,rept,esttau,full(Likelihood(end)),full(MLikelihood(end)))
            end
            MLikelihood = MLikelihood(~isinf(MLikelihood));
            Mlike(rept)=MLikelihood(end);
        end
    else
        for rept=1:fold
            missing=cvsample;
            missing(missing==rept)=0;
            missing(missing~=0)=1;
            if p.Results.savemat
                save = rept;
            else
                save = [];
            end
            [~,~,~,esttau,Likelihood,MLikelihood]=ZIPFA(X, nf, 0, p.Results,'missing',logical(missing),'rept',save);
            if p.Results.display
                fprintf('\n factor = %d; fold = %d; esttau = %.3g; Likelihood = %.5g; MLikelihood = %.4g \n',nf,rept,esttau,full(Likelihood(end)),full(MLikelihood(end)))
            end
            MLikelihood = MLikelihood(~isinf(MLikelihood));
            Mlike(rept)=MLikelihood(end);
        end
    end
    Allres=[Allres Mlike];
end
end