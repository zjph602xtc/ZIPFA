function [ res ] = dlnL( beta, Dx, y, Z )
tau = beta(1);
beta = beta(2:end);
xb = Dx * beta;
n = length(y);
a = (1 - Z) * spdiags(y-exp(xb),0,n,n) * Dx;
b = ((1./(exp(tau * xb)+1))'-Z) * xb;
res = [a b];
end

