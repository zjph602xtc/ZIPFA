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
        + (-1+Z) * (y.*xb  -mexb);
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