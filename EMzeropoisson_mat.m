function [startvar] = EMzeropoisson_mat(data, varargin)
% EMzeropoisson_mat  To fit the zero inflated Poisson regression.
%   fittedbeta = EMzeropoisson_mat([y x])
%
%   fittedbeta = EMzeropoisson_mat([y x], tau, 'display', false, ...)
%
% -data: First y then x.
% -tau (0.1): Initial tau to fit. Will be overwritten by 'initial'.
% -'initial' ([]): Initial [tau beta].
% -'initialtau' ('iteration'): Choose the initial value of tau at the beginning of EM iteration.
%       'stable': estimate tau from fitted beta in last round;
%       'initial': always use the initially assigned tau in 'tau' or 'initial';
%           Use the default tau = 0.1 if 'initial' is empty.
%       'iteration': use fitted tau in last round.
% -'tol' (1e-4): Percentage of l2 norm change of [tau beta].
% -'maxiter' (100): Max iteration.
% -'Madj' (false): Whether adjust for relative library size M.
% -'m' ([]): Relative library size M.
% -'display' (true): Display the fitting procedure.
% -'intercept' (true): Whether the model contains an intercept.
%
% Result contains the fitted results in each row. The last row shows
% the final result. First column is tau, second column is intercept (if the
% model has intercept), other columns are fitted coefficients.
%
%   See also ZIPFA, cv_ZIPFA

% Tianchen Xu


[n, nvar] = size(data);
p = inputParser;
addRequired(p,'data',@ismatrix);
addOptional(p,'tau',0.1,@isscalar);
addParameter(p,'initial',[],@isvector);
addParameter(p,'initialtau','iteration',@(x)any(validatestring(x,{'stable','initial','iteration'})));
addParameter(p,'tol',1e-04,@(x)validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'maxiter',100,@(x)validateattributes(x,{'numeric'},{'scalar','>=',2}));
addParameter(p,'Madj',false,@islogical);
addParameter(p,'m',[],@isvector);
addParameter(p,'display',true,@islogical);
addParameter(p,'intercept',true,@islogical);

if length(dbstack)>1
    p.KeepUnmatched = true;
end
parse(p,data,varargin{:});
if p.Results.Madj 
    validateattributes(p.Results.m,{'numeric'},{'size',[n,1]},'EMzeropoisson_mat','m');
end

x = data(:,2:end);
y = data(:,1);
tau = p.Results.tau;
if p.Results.intercept; intercept = 'on'; else intercept = 'off';end
display = p.Results.display;
% taubackup = p.Results.tau;

if isempty(p.Results.initial)
    if display
        disp('Initializing ...')
    end
    startvar = glmfit(x,y,'poisson','link','log','constant',intercept);
    startvar = [tau; startvar]';
else
    startvar = p.Results.initial;
end

if strcmp(intercept,'on')
    Dx = [ones(n,1) x];
else
    Dx = x;
end

i = 1;
if display
    disp('Start maximizing ...')
end

while (i < p.Results.maxiter && (i==1 || ...
        sqrt(sum((startvar(i,:)-startvar(i-1,:)).^2))/sqrt(sum(startvar(i,:).^2))>p.Results.tol))
    % E step
    Z = zeros(1, n);
    Izero = y==0;
    if p.Results.Madj
        m = p.Results.m;
        Z(Izero) = 1./(1 + exp( startvar(i,1).*(Dx(Izero,:) * startvar(i,2:end)')-exp(m(Izero).*(Dx(Izero,:) * startvar(i,2:end)'))  ));
    else
        Z(Izero) = 1./(1 + exp( startvar(i,1).*(Dx(Izero,:) * startvar(i,2:end)')-exp(Dx(Izero,:) * startvar(i,2:end)')  ));
    end
    
    
    % M step
    options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,...
        'HessianFcn','objective','Display','off');
    if (i ~=1 && strcmp(p.Results.initialtau,'stable'))
        tau = -mean(log( n./ sum(Z==0) - 1) ./ (Dx * startvar(i,2:end)'));
    elseif (strcmp(p.Results.initialtau,'initial'))
        tau = startvar(1, 1);
    elseif (strcmp(p.Results.initialtau,'iteration'))
        tau = startvar(i, 1);
    end
    
    if p.Results.Madj
        startvar = [startvar; fminunc(@lnL, [tau startvar(i,2:end)]', options, Dx,y,Z,m)'];
    else
        startvar = [startvar; fminunc(@lnL, [tau startvar(i,2:end)]', options, Dx,y,Z,1)'];
        %         startvar=[startvar; fsolve(@dlnL,[tau startvar(i,2:end)]', optimoptions('fsolve','display','off',...
        %             'Algorithm', 'levenberg-marquardt','StepTolerance',1e-4,'InitDamping',50),Dx,y,Z)'];
    end
    
    i = i+1;
    if (strcmp('iteration',p.Results.initialtau) && abs(max(diff(startvar(:,1))))>50 && i>3)
        disp('May be divergent tau. Try ''stable'' model.')
        %%continue
    end
    if display
        fprintf('This is %.0f th iteration, Frobenius norm diff = %g. \n',...
            i, sqrt(sum((startvar(i,:)-startvar(i-1,:)).^2))/sqrt(sum(startvar(i,:).^2)))
    end
end
end

function [value, grad, hess] = lnL( beta, Dx,y,Z,m )
tau = beta(1);
beta = beta(2:end);
xb = Dx * beta;
mexb = m.*exp(xb);
etauxb = exp(tau.*xb);
% n = length(y);

value = tau.* Z*xb + sum(log(1+exp(-tau.*xb))) ...
    + (-1+Z) * (y.*xb  -mexb);
% value = -value

if nargout > 1 % gradient required
    %     a = (1 - Z) * spdiags(y-exp(xb),0,n,n) * Dx;
    wz = (1./(etauxb+1))'-Z;
    dldb = (wz.*tau + (1 - Z) .* (y-mexb)') * Dx;
    dldtau = wz * xb;
    grad = -[dldtau dldb]'; % value = -value
    if nargout > 2 % Hessian required
        dldtautau = -sum( xb.^2.*etauxb./(etauxb+1).^2 );
        dldtaudb = (transpose((etauxb+1-tau.*xb.*etauxb)./(etauxb+1).^2)-Z)*Dx;
        dldbdb =  bsxfun(@times, Dx',tau*tau.*transpose(-etauxb./(etauxb+1).^2)...
            +(1-Z).* transpose(-m.*exp(xb)))*Dx;
        hess = -[[dldtautau dldtaudb]; [dldtaudb' dldbdb]];
    end
end

if isinf(value) ||  isnan(value)
    %     sprintf('Warning: objective value is %.5g',value)
    small = (-tau.*xb)<30;
    value = tau.* Z*xb + sum(log(1+exp(-tau.*xb)).*small + (-tau.*xb).*(~small),'omitnan') ...
        + (-1+Z) * (y.*xb - mexb);
end

if any(any(isinf(hess))) ||  any(any(isnan(hess)))
    %     sprintf('Warning: hess contains bad value')
    small = etauxb < exp(30);
    dldtautau = -sum( xb.^2.*etauxb./(etauxb+1).^2.*small ...
        + xb.^2./etauxb.*(~small),'omitnan' );
    v1=(etauxb+1-tau.*xb.*etauxb)./(etauxb+1).^2;
    v2=(1-tau.*xb)./etauxb;
    V=zeros(length(y),1);
    V(small)=v1(small);
    V(~small)=v2(~small);
    dldtaudb = (V'-Z)*Dx;
    R=zeros(length(y),1);
    r1=etauxb./(etauxb+1).^2;
    r2=1./etauxb;
    R(small)=r1(small);
    R(~small)=r2(~small);
    dldbdb =  bsxfun(@times, Dx',-tau*tau.*R'+(1-Z).* transpose(-m.*exp(xb)))*Dx;
    hess = -[[dldtautau dldtaudb]; [dldtaudb' dldbdb]];
end
end