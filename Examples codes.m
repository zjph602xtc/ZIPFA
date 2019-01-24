%% Zero inflated Poisson regression
% generate data
n = 50000;
x1 = randn(n,1);
x2 = randn(n,1);
lam = exp(x1 - 2*x2 + 1.5);
y = poissrnd(lam, n, 1);
tau = .75; % tau = 0.75
p = 1./(1+lam.^tau);
mean(p)
Z = binornd(1, p, n ,1);
y(logical(Z)) = 0;
max(y)
clear n lam tau p Z ans

% fit the model 
res = EMzeropoisson_mat([y x1 x2]);
fittedtau = res(end,1);
fittedintercept = res(end,2);
fittedbeta = res(end,3:end);

EMzeropoisson_mat([y x1 x2],0.5,'initialtau','stable','tol',0.01)

%% Zero inflated Poisson factor analysis
% generate data
rng(1)
u = [[0*ones(35,1); 2*ones(45,1); 1.7*ones(60,1);0*ones(60,1)]...
    [1.8*ones(35,1);0.9*ones(45,1); 0*ones(120,1)]...
    [0*ones(35,1);1.7*ones(165,1)]];
vt = [[0.0*ones(1,30) 0*ones(1,30) 1.7*ones(1,40)];...
    [0*ones(1,35) 1.7*ones(1,25) 1*ones(1,40)];...
    [1.7*ones(1,25) 0.9*ones(1,50) 0.9*ones(1,25)]];
u=normrnd(0,0.06,200,3)+u;
vt=normrnd(0,0.05,3,100)+vt;
a = u * vt;

lambda = exp(a);
X = poissrnd(lambda,200,100);

tau = 0.616; % tau = 0.616
P = 1./(1+lambda.^tau);
Z = binornd(1,P,200,100);
X(logical(Z)) = 0;
clear a lambda P u vt Z tau

% fit the model
[Ufit, Vfit, itr, fittedtau, likelihood]= ZIPFA(X,3,'Madj',false);
fittedU = Ufit{itr};
fittedV = Vfit{itr};

% calculate the likelihood of first 2 rows in the predicted model
[~, ~, ~, ~, likelihood, CVlikelihood]= ZIPFA(X,3,'Madj',false,...
    'missing',[false(2,102);true(198,102)]);

%% Cross validation - Zero inflated Poisson factor analysis
% use the simulated data in last section
% use parallel toolbox accelerate
[cvsample,CVlikelihood]=cv_ZIPFA(X,1:6,10,'Madj',false);
errorbar(1:6,mean(CVlikelihood),std(CVlikelihood))
% single thread
% [cvsample,CVlikelihood]=cv_ZIPFA(X,2:4,10,'Madj',false,'parallel',false);
% errorbar(2:4,mean(CVlikelihood),std(CVlikelihood))

% add one more rank = 7, with exising cvsample (do not delete variable cvsample)
[~,CVlikelihood7]=cv_ZIPFA(X,7,10,'Madj',false,'display',false);

% By Tianchen Xu  01/14/2019