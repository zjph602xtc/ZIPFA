function [Ufit,Vfit,itr,finaltau,Likelihood,MLikelihood] = ZIPSVD(X, k, iter, Madj, tau, tol ...
    ,tolLnlikelihood, maxiter, initialtau, display, initialmode, missing,rept)
% X is the matrix to be decomposed
% k is the number of factor
% iter is the max iteration
% Madj: whether adjust for m
% tau is the initial guess for tau
% tolLnlikelihood is the max difference in two iterations
% tol is percentage of max parameter difference in two iterations searching for parameters
% initialtau: 'iteration', 'stable', 'initial'.
% maxiter is the max iterations in two iterations searching fro parameter
% display: whether display process
% initialmode: 'SVD', 'Poisson'. first round fitting mode
% missing: if x missing is not null, Then likelihood is the rest of X !!!! T or F matrix
Xt = X';
[m,n] = size(X);
if any(isnan(X))
    if (~isempty(missing))
        error('Cannot assign ''missing'' when X has missing values')
    end
    indexram = ~isnan(X(:));
    indexram_row = ~isnan(Xt(:));
else
    if (~isempty(missing))
        indexram = logical(missing(:));
        missingt = missing';
        indexram_row = logical(missingt(:));
    else
        indexram = 1:(m*n);
        indexram_row = indexram;
    end
end
clear missingt

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

% interation begins
Ufit = cell(iter,1);
Vfit = cell(iter,1);
Likelihood=[];
MLikelihood=[];

if Madj
    ratio = repmat(sum(X,2,'omitnan')/median(sum(X,2,'omitnan')),n,1);
    ratio_row = reshape(ratio,m,[])';
    ratio_row = ratio_row(:);
end

for itr=1:iter
    fprintf('\n ****************** \n Round %.0f \n ******************\n', itr)
    %%% update U
    if Madj
        mi = ratio_row(indexram_row);
    end
    Voldc = mat2cell(sparse(Vold),n,k);
    Voldc = repmat(Voldc, m, 1);
    dat = [Xt(:) blkdiag(Voldc{:})];
    dat = dat(indexram_row,:);
    if (strcmp(initialmode,'poisson'))
        return %%%
    else
        if itr==1
            uuuu = EMzeropoisson_mat(dat, tau, tol, maxiter, initialtau, Madj, display, [tau reshape(Uold',1,[])], 'off');
        else
            uuuu = EMzeropoisson_mat(dat, tau, tol, maxiter, initialtau, Madj, display, [uuuu(end,1) reshape(Ufit{itr-1}',1,[])], 'off');
        end
    end
%     lnL_mat([Xt(:) blkdiag(Voldc{:})], uuuu, indexram_row, ratio_row, 'off')
    Unew = reshape(uuuu(end,2:end), k,[])';
    
    %%% update V
    if Madj
        mi = ratio(indexram);
    end
    Uoldc = mat2cell(sparse(Uold),m,k);
    Uoldc = repmat(Uoldc, n, 1);
    dat = [X(:) blkdiag(Uoldc{:})];
    dat = dat(indexram,:);
    if (strcmp(initialmode,'poisson'))
        return %%%
    else
        if itr==1
            uuuu = EMzeropoisson_mat(dat, tau, tol, maxiter, initialtau, Madj, display, [tau reshape(Vold',1,[])], 'off');
        else
            uuuu = EMzeropoisson_mat(dat, tau, tol, maxiter, initialtau, Madj, display, [uuuu(end,1) reshape(Vfit{itr-1}',1,[])], 'off');
        end
    end
    Vnew = reshape(uuuu(end,2:end), k,[])';
    
    if Madj
        Likelihood=[Likelihood lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, indexram, ratio, 'off')];
        if ~isempty(missing)
            MLikelihood=[MLikelihood  lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, ~indexram, ratio, 'off',rept,k)];
        end
    else
        Likelihood=[Likelihood lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, indexram, [], 'off')];
        if ~isempty(missing)
            MLikelihood=[MLikelihood  lnL_mat([X(:) blkdiag(Uoldc{:})], uuuu, ~indexram, [], 'off')];
        end
    end
    
    disp(Likelihood)
    disp(MLikelihood)
    %%% next step
    [U,S,V] = svd(Unew * Vnew');
    Uold = U(:,1:k) * S(1:k,1:k);
    Vold = V(:,1:k);
    
    %%% save answer
    Ufit{itr} = Uold;
    Vfit{itr} = Vold;
    
    if (itr~=1)
        fprintf('Max Ln likelihood diff = %.4g %%',full(100*(Likelihood(itr)-Likelihood(itr-1))/abs(Likelihood(itr-1))))
        if (full((Likelihood(itr)-Likelihood(itr-1))/abs(Likelihood(itr-1))))<tolLnlikelihood
            finaltau=uuuu(end,1);
            return
        end
    end
end
finaltau=uuuu(end,1);
end

