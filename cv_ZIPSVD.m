function [cvsample,Allres] = cv_ZIPSVD (X, k, fold, iter, Madj, tau, tol ,tolLnlikelihood, maxiter, initialtau, display, initialmode)
if any(isnan(X))
    error('There is Na''s in X !')
end
Allres = [];
[m, n] =size(X);
% X=X(:,sum(X==0)/m<0.9);
% generate assignment
cvsample=[];
for i=1:fold-1
    cvsample=[cvsample i*ones(1,floor(m*n/fold))];
end
cvsample=[cvsample fold*ones(1,m*n-(fold-1)*floor(m*n/fold))];
cvsample=datasample(cvsample,length(cvsample),'Replace',false);
cvsample=reshape(cvsample,m,n);
%%%%
try 
    cvsample=evalin('base','cvsample');
    disp('Use existing cvsample!!!')
end
% 
%%%%%
for nf = k
    Mlike=zeros(fold,1);
    parfor rept=1:fold
        missing=cvsample;
        missing(missing==rept)=0;
        missing(missing~=0)=1;
        [Ufit,Vfit,itr,esttau,Likelihood,MLikelihood]=ZIPSVD(X, nf, iter, Madj, tau, tol ,...
            tolLnlikelihood, maxiter, initialtau, display, initialmode, missing,rept);
        fprintf('\n itr = %.0g; esttau = %.3g; Likelihood = %.4g; MLikelihood = %.4g \n',itr,esttau,full(Likelihood(end)),full(MLikelihood(end)))
        Mlike(rept)=MLikelihood(end);
    end
    Allres=[Allres [nf; mean(Mlike)]];
end
end