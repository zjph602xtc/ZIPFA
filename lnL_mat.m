function [likelihood] = lnL_mat(data, coefficient, subset, mi, intercept,rept,k)
% note: data and mi are in full data size.
coefficient = coefficient(end,:);
tau = coefficient(1);
beta = coefficient(2:end)';
y = data(:,1);
n = length(y);
if strcmp(intercept,'on')
    Dx = [ones(n,1) data(:,2:end)];
else
    Dx = data(:,2:end);
end
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
% like(indexnotz) = log(1-p(indexnotz)) - mi(indexnotz).*exp(lnlam(indexnotz))...
%     +y(indexnotz).*(log(mi(indexnotz))+lnlam(indexnotz));
if ~isempty(subset)
    like = like(subset);
end
likelihood = sum(like);

%%%%
try
% Alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
% range=[Alphabet(k) '2'];
% xlswrite('mdetail.xlsx',like,rept,range); 
save(['factor' num2str(k) 'fold' num2str(rept)])
end
%%%%
end