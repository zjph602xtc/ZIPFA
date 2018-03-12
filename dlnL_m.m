function [ res ] = dlnL_m( beta, Dx, y, Z, m)
tau = beta(1);
beta = beta(2:end);
xb = Dx * beta;
n = length(y);
a = (1 - Z) * spdiags(y-m.*exp(xb),0,n,n) * Dx;
b = sum(((1./(exp(tau * xb)+1))'-Z) * xb);
res = [a b];
end

