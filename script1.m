n = 50000;
x1 = randn(n,1);
x2 = randn(n,1);

lam = exp(x1 - 2*x2 + 1.5);
y = poissrnd(lam, n, 1);

tau = .75;
p = 1./(1+lam.^tau);
mean(p)
Z = binornd(1, p, n ,1);
y(logical(Z)) = 0;
max(y)




EMzeropoisson_mat([y x1 x2], .5, 2e-4, 100,'iteration',false,true, [0 1 1 -5 ],'on')

%%
rng(250)
u = [[2.5*ones(35,1); 2*ones(45,1); 1.7*ones(60,1);0*ones(60,1)] [0.9*ones(60,1); 0*ones(70,1); 00*ones(70,1)] [0*ones(35,1);1.9*ones(125,1); 1.8*ones(40,1)]];
vt = [[0.0*ones(1,30) 0*ones(1,30) 1.5*ones(1,40)]; [0*ones(1,35) 1.5*ones(1,25) 1*ones(1,40)]; [1.3*ones(1,25) 0.8*ones(1,50) 0*ones(1,25)]];
u=normrnd(0,0.08,200,3)+u;
vt=normrnd(0,0.06,3,100)+vt;
a = u * vt;
% HeatMap(a,'Symmetric',false)

lambda = exp(a);
max(max(lambda))
X = poissrnd(lambda,200,100);

tau = 0.28;
P = 1./(1+lambda.^tau);
mean(mean(P))


Z = binornd(1,P,200,100);
X(logical(Z)) = 0;

[Ufit,Vfit,itr,esttau,Likelihood,MLikelihood]=ZIPSVD(X, 3, 10, false, 0.4, 1e-03 ,...
    0.0001, 100, 'iteration', true, 'SVD', []);
[Ufit,Vfit,itr,esttau,Likelihood,MLikelihood]=ZIPSVD(X, 5, 10, true, 0.25, 1e-04 ,...
    0.03, 100, 'iteration', true, 'SVD', []);
subplot(1,2,1)
heatmap(u*vt,[],[],[],'Colormap',hot1,'MinColorValue',-0.2,'MaxColorValue',5.6);
subplot(1,2,2)
heatmap(Ufit{itr}*Vfit{itr}',[],[],[],'Colormap',hot1,'MinColorValue',-0.2,'MaxColorValue',5.6);


X=X(:,sum(X==0)/281<0.8);
[Ufit,Vfit,itr,esttau,Likelihood,MLikelihood]=ZIPSVD(X, 5, 20, true, 0.25, 1e-04 ,...
    0.00001, 100, 'iteration', true, 'SVD', []);

[cvsample,allres] = cv_ZIPSVD (X, 4:6, 8, 10, true, 0.25, 1e-4 ,0.00001, 35,...
    'iteration', true, 'SVD');
system('shutdown -s -t 300')
%%% real data
