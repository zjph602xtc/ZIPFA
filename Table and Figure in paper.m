%% colorbar
hot1=flipud(hot(300));
hot1=hot1(1:225,:);

%% Web Figure 1 (a)
% plot U,V
rng(250)
u = [[0*ones(35,1); 2*ones(45,1); 1.7*ones(60,1);0*ones(60,1)]...
    [1.8*ones(35,1);0.9*ones(45,1); 0*ones(120,1)]...
    [0*ones(35,1);1.7*ones(165,1)]];
vt = [[0.0*ones(1,30) 0*ones(1,30) 1.7*ones(1,40)];...
    [0*ones(1,35) 1.7*ones(1,25) 1*ones(1,40)];...
    [1.7*ones(1,25) 0.9*ones(1,50) 0.9*ones(1,25)]];
u=normrnd(0,0.06,200,3)+u;
vt=normrnd(0,0.05,3,100)+vt;
subplot(1,2,1)
hx=heatmap(u,[],[],[],'Colormap',hot1,'MinColorValue',0,'MaxColorValue',2.5);
title('U');
subplot(1,2,2)
hx=heatmap(vt',[],[],[],'Colormap',hot1,'MinColorValue',0,'MaxColorValue',2.5);
title('V');
colorbar;

%% Web Figure 1 (b)
rng(250)
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
max(max(lambda))
X = poissrnd(lambda,200,100); % save as uvax.mat

hx=heatmap(a,[],[],[],'Colormap',hot1,'MinColorValue',-0.2,'MaxColorValue',6.6);
title('true ln(\lambda)')
colorbar('southoutside')
h=axes(gcf)
htt=histogram(h,a,60,'FaceAlpha',0,'EdgeColor','blue','LineWidth',0.7)
h.Visible='off';
h.XLim=[-0.2,6.6];
h.YLim=[0,2050];

%% Web Figure 1 (c)
color = [[30 129 196]/255; % b
    [237 77 36]/255;  %r
    [0,175,79]/255;   %g
    [237 177 36]/255;  %y
    [182 122 199]/255;   %p
    [22, 171, 192]/255;]; % lb
subplot(2,1,1)
Z = linkage(u,'complete');
d=dendrogram([Z(:,1:2) log(log(Z(:,3)/min(Z(:,3))/0.95*exp(1)))],200,...
    'reorder',1:200,'ColorThreshold',1.7,'checkcrossing',false)
set(d,'color','k')
ax1=gca;
ax1.Visible='off';
subplot(2,1,2)
pt = plot(u,'+');
set(pt(1),'Marker','o','MarkerSize',3,'MarkerEdgeColor',color(1,:),'MarkerFaceColor',color(1,:))
set(pt(2),'Marker','d','MarkerSize',3,'MarkerEdgeColor',color(2,:),'MarkerFaceColor',color(2,:))
set(pt(3),'Marker','s','MarkerSize',3,'MarkerEdgeColor',color(3,:),'MarkerFaceColor',color(3,:))

xlabel('Sample index')
ylabel('Factor value')
l1=legend({'Factor 1','Factor 2','Factor 3'});
outer=get(gca,'OuterPosition');
set(gca,'OuterPosition',[outer(1:3),0.72])
box('off')
f=gcf;
set(f,'Position',[285,326,495,479])
l1.Position=[0.7052    0.5798    0.1529    0.1075];
pnew

%% Web Figure 1 (d)
subplot(2,1,1)
Z = linkage(vt','complete');
d=dendrogram([Z(:,1:2) log(log(Z(:,3)/min(Z(:,3))/0.95*exp(1)))],200,...
    'reorder',1:100,'ColorThreshold',1.7,'checkcrossing',false)
set(d,'color','k')
ax1=gca;
ax1.Visible='off';
subplot(2,1,2)
pt=plot(vt','+');
set(pt(1),'Marker','o','MarkerSize',3,'MarkerEdgeColor',color(1,:),'MarkerFaceColor',color(1,:))
set(pt(2),'Marker','d','MarkerSize',3,'MarkerEdgeColor',color(2,:),'MarkerFaceColor',color(2,:))
set(pt(3),'Marker','s','MarkerSize',3,'MarkerEdgeColor',color(3,:),'MarkerFaceColor',color(3,:))

xlabel('Taxa index')
ylabel('Score value')
l1=legend({'Score 1','Score 2','Score 3'});
outer=get(gca,'OuterPosition');
set(gca,'OuterPosition',[outer(1:3),0.72])
box('off')
f=gcf;
set(f,'Position',[285,326,495,479])
l1.Position=[0.7052    0.1298    0.1529    0.1075];
pnew

%% Web Figure 3
load('alldataX.mat')
[~,CVlikelihood]=cv_ZIPFA(X,1:8,10,'display',false);

%% Figrue 1
load('alldataX.mat')
% X=X(:,sum(X==0)/281<0.8);
m = sum(X,2,'omitnan')/median(sum(X,2,'omitnan'));

Xd = bsxfun(@(x,y)x./y, X, m);
X1 = log(Xd);
X1(isinf(X1)) = NaN;
X1mean = mean(X1, 'omitnan');

p = sum(X==0)/281;
pz=p==0;
plot(X1mean(~pz),p(~pz),'+')
% plot(X1mean(~pz),log(p(~pz)./(1-p(~pz))),'+')
axis([-1 7 -0.02 1.02])

esttau = fitlm(X1mean(~pz),log(p(~pz)./(1-p(~pz))),'Intercept',false);

%% Figure 4 (a)(b)
load('alldataX.mat')
tic
[U,V,itr,esttau] = ZIPFA(X,5)
toc
u = U{itr}; v = V{itr}; % save as 'finaluv.mat';
load('finaluv.mat')
lnlam = u*v';
p = 1./(1+exp(lnlam).^esttau);
N =sum(X,2,'omitnan')/median(sum(X,2,'omitnan'));
p = p+(1-p).*poisspdf(0, bsxfun(@times,exp(lnlam),N));

row = sum((p).^(1/100),1);
[~,I] = sort(row);
pp = p(:,I);

[pp,Ir] = sort(pp,1);
figure(1)
hx=heatmap((1-pp).^2.8,[],[],[],'Colormap',hot1,'MinColorValue',0,'MaxColorValue',.91);
co = colorbar;
co.TickLabels = arrayfun(@(x)num2str(x,'%.3f'),1-[0:0.1:0.9].^(1/2.8),'unif',0);

X=X(:,sum(X==0)/281<0.8);
XX = X(:,I);
for i = 1:297
    XX(:,i) = XX(Ir(:,i),i);
end
figure(2)
hx=heatmap((log1p(XX)),[],[],[],'Colormap',hot1,'MinColorValue',1.5,'MaxColorValue',11.6);
co = colorbar;
co.Ticks = [1.5 2:11]
co.TickLabels = ['0.0' arrayfun(@(x)num2str(x,'%4.1f'),exp(2:11),'unif',0)];

XX(XX<20)=0;
XX(XX~=0)=1;
figure(3)
hx=heatmap(XX,[],[],[],'Colormap',hot1,'MinColorValue',0,'MaxColorValue',1.5);


%% Table 1: generate simulation data
% setting 0: 0 percentage  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    csvwrite(sprintf('g1_%d.csv',i),X);
end

% setting 1  20%  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = 0.616;
    P = 1./(1+lambda.^tau);
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g2_%d.csv',i),X);
end


% setting 1 40% only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = 0.157;
    P = 1./(1+lambda.^tau);
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g3_%d.csv',i),X);
end

% setting 2 20%  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = 0.192;
    P = exp(-lambda.^tau);
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g4_%d.csv',i),X);
end

% setting 2  40%   only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = -0.033;
    P = exp(-lambda.^tau);
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g5_%d.csv',i),X);
end

% setting 3 20%  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = 0.708;
    P = 1-exp(-lambda.^(-tau));
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g6_%d.csv',i),X);
end

% setting 3 40%   only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = 0.272;
    P = 1-exp(-lambda.^(-tau));
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g7_%d.csv',i),X);
end

% setting 4 20%  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = 0.222;
    P = exp(-lambda.*tau);
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g8_%d.csv',i),X);
end

% setting 4  40%  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    tau = 0.0846;
    P = exp(-lambda.*tau);
    Z = binornd(1,P,200,100);
    X(logical(Z)) = 0;
    csvwrite(sprintf('g9_%d.csv',i),X);
end

% seting 5 20%  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    p = unifrnd((0.2-sum(sum(X==0))/20000-0.1),( 0.2-sum(sum(X==0))/20000+0.1),1,100);
    %p=normrnd(0.2-sum(sum(X==0))/20000,0.05,1,100);
    p=repmat(p,200,1);
    Z = binornd(1,abs(p),200,100);
    X(logical(Z)) = 0;
    sum(sum(X==0))/20000;
    csvwrite(sprintf('g10_%d.csv',i),X);
end

% setting 5  40%  only counts inflated zeros
for i=1:200
    rng(i)
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
    max(max(lambda));
    X = poissrnd(lambda,200,100);
    p = unifrnd((0.4-sum(sum(X==0))/20000-0.1),( 0.4-sum(sum(X==0))/20000+0.1),1,100);
    %p=normrnd(0.4-sum(sum(X==0))/20000,0.05,1,100);
    p=repmat(p,200,1);
    Z = binornd(1,abs(p),200,100);
    X(logical(Z)) = 0;
    sum(sum(X==0))/20000;
    csvwrite(sprintf('g11_%d.csv',i),X);
end

%% Table 1: calculate values
% ZIPFA
for g=1:11
    for i=1:200
        X(:,:,i)=csvread(sprintf('g%d_%d.csv',g,i));
    end
    parfor i=1:200
        [Ufit,Vfit,itr,esttau,Likelihood,MLikelihood]=ZIPFA(X(:,:,i), 3, 'Madj',false);
        fita(:,:,i)=Ufit{itr}*Vfit{itr}';
    end
    for i=1:200
        csvwrite(sprintf('g%d_%d_fit.csv',g,i),fita(:,:,i));
    end
end

% log-PCA
for g=10:11
    for i=1:200
        X(:,:,i)=csvread(sprintf('g%d_%d.csv',g,i));
    end
    X(X==0)=0.5;
    parfor i=1:200
        [coeff,score]=pca(log(X(:,:,i)),'centered',false,'NumComponents',3);
        fita(:,:,i)=score*coeff';
    end
    for i=1:200
        csvwrite(sprintf('g%d_%d_fit_log.csv',g,i),fita(:,:,i));
    end
end

% gomms, PSVDOS   see R code 'Simu.R'

%% Table 1: collect results
% ZIPFA/log-PCA
u = [[0*ones(35,1); 2*ones(45,1); 1.7*ones(60,1);0*ones(60,1)]...
    [1.8*ones(35,1);0.9*ones(45,1); 0*ones(120,1)]...
    [0*ones(35,1);1.7*ones(165,1)]];
vt = [[0.0*ones(1,30) 0*ones(1,30) 1.7*ones(1,40)];...
    [0*ones(1,35) 1.7*ones(1,25) 1*ones(1,40)];...
    [1.7*ones(1,25) 0.9*ones(1,50) 0.9*ones(1,25)]];
a = u * vt;
for g=1:11
    for i=1:200
        X(:,:,i)=csvread(sprintf('g%d_%d_fit_log.csv',g,i));
%       X(:,:,i)=csvread(sprintf('g%d_%d_fit.csv',g,i));
        l2(g,i)=sum(mean(abs(X(:,:,i)-a)).^2);
        [taxa(g,i),sub(g,i)]=cluster_fit(X(:,:,i));
        lsum=sum(sum((X(:,:,i)-a).^2));
        dif=sum((X(:,:,i)-a).^2);
        simp(g,i)=sum((dif./lsum).^2);
        bray(g,i)=sum(mean(abs(X(:,:,i)-a)))/2;
    end
end
mean(l2,2,'omitnan')
std(l2,0,2,'omitnan')
mean(bray,2,'omitnan')
std(bray,0,2,'omitnan')
mean(taxa,2,'omitnan')
std(taxa,0,2,'omitnan')
mean(sub,2,'omitnan')
std(sub,0,2,'omitnan')
mean(1000*simp,2,'omitnan')
std(1000*simp,0,2,'omitnan')

% GOMMS
u = [[0*ones(35,1); 2*ones(45,1); 1.7*ones(60,1);0*ones(60,1)]...
    [1.8*ones(35,1);0.9*ones(45,1); 0*ones(120,1)]...
    [0*ones(35,1);1.7*ones(165,1)]];
vt = [[0.0*ones(1,30) 0*ones(1,30) 1.7*ones(1,40)];...
    [0*ones(1,35) 1.7*ones(1,25) 1*ones(1,40)];...
    [1.7*ones(1,25) 0.9*ones(1,50) 0.9*ones(1,25)]];
a = u * vt;
for g=11:11
    for i=1:200
        X(:,:,i)=csvread(sprintf('g%d_%d_fit_gomms.csv',g,i));
        l2(g,i)=sum(mean(abs(X(:,:,i)-a)).^2);
        [taxa(g,i),sub(g,i)]=cluster_fit(X(:,:,i));
        lsum=sum(sum((X(:,:,i)-a).^2));
        dif=sum((X(:,:,i)-a).^2);
        simp(g,i)=sum((dif./lsum).^2);
        bray(g,i)=sum(mean(abs(X(:,:,i)-a)))/2;
    end
end
dd=(l2==l2(10,1));
dd(5,28)=1;
l2(dd)=nan;bray(dd)=nan;taxa(dd)=nan;sub(dd)=nan;simp(dd)=nan;
mean(l2,2,'omitnan')
std(l2,0,2,'omitnan')
mean(bray,2,'omitnan')
std(bray,0,2,'omitnan')
mean(taxa,2,'omitnan')
std(taxa,0,2,'omitnan')
mean(sub,2,'omitnan')
std(sub,0,2,'omitnan')
mean(1000*simp,2,'omitnan')
std(1000*simp,0,2,'omitnan')

% PSVDOS
u = [[0*ones(35,1); 2*ones(45,1); 1.7*ones(60,1);0*ones(60,1)]...
    [1.8*ones(35,1);0.9*ones(45,1); 0*ones(120,1)]...
    [0*ones(35,1);1.7*ones(165,1)]];
vt = [[0.0*ones(1,30) 0*ones(1,30) 1.7*ones(1,40)];...
    [0*ones(1,35) 1.7*ones(1,25) 1*ones(1,40)];...
    [1.7*ones(1,25) 0.9*ones(1,50) 0.9*ones(1,25)]];
a = u * vt;
for g=1:11
    for i=1:200
        X(:,:,i)=csvread(sprintf('g%d_%d_fit_psvdos.csv',g,i));
        l2(g,i)=sum(mean(abs(X(:,:,i)-a)).^2);
        [taxa(g,i),sub(g,i)]=cluster_fit(X(:,:,i));
        lsum=sum(sum((X(:,:,i)-a).^2));
        dif=sum((X(:,:,i)-a).^2);
        simp(g,i)=sum((dif./lsum).^2);
        bray(g,i)=sum(mean(abs(X(:,:,i)-a)))/2;
    end
end
dd=(l2==0);
l2(dd)=nan;bray(dd)=nan;taxa(dd)=nan;sub(dd)=nan;simp(dd)=nan;
mean(l2,2,'omitnan')
std(l2,0,2,'omitnan')
mean(bray,2,'omitnan')
std(bray,0,2,'omitnan')
mean(taxa,2,'omitnan')
std(taxa,0,2,'omitnan')
mean(sub,2,'omitnan')
std(sub,0,2,'omitnan')
mean(1000*simp,2,'omitnan')
std(1000*simp,0,2,'omitnan')

%% Figure 5
load('finaluv.mat')
taxa = readtable('taxa.csv','ReadVariableNames',0);
taxa = taxa(sum(X==0)/281<0.8,:);
taxa = taxa{:,:};
% use results in 'our_result.xlsx'
% fac: loading on factor 2
% sng: whether is clinical meaningful
h=bar(fac);
color = [[30 129 196]/255;
    [237 77 36]/255;
    [0,175,79]/255];

fHand = figure;
aHand = axes('parent', fHand);
hold(aHand, 'on')
h = [];
for i = 1:numel(fac)
    if sng(i)==-1
        c = 3;
    elseif sng(i)==0
        c = 1;
    else
        c = 2;
    end
    h = [h bar(i, fac(i), 'parent', aHand, 'facecolor', color(c,:))];
end

set(h,'edgecolor','none')
xlim([0 231])

xlabel('Taxa species')
ylabel('Absolute loadings on factor 2')

legend(aHand,h([1 3 2]),{'Clinical meaningful taxa with negative loadings' 'Clinical meaningful taxa with positive loadings' 'Other taxa'})
set(gcf,'PaperSize',[11.5000 11])

% permutation test
median(find(sng))
iter = 50000;
sam = [ones(25,iter); zeros(204,iter)];
f = [];
for i = 1:iter
    sam(:,i)=randsample(sam(:,1),229);
    f = [f find(sam(:,i))];
end
hist(mean(f))
sum(mean(f)<mean(find(sng)))

%% Figure 3
% generate simulation data
taun = 1;
tauest = [];
flag = true;
for j = 1:15
    for i=1:200
        rng(i)
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
        max(max(lambda));
        X = poissrnd(lambda,200,100);
        if j~=1
            while (flag)
                Xback = X;
                P = 1./(1+lambda.^taun);
                Z = binornd(1,P,200,100);
                Xback(logical(Z)) = 0;
                if (sum(sum(Xback==0))/20000 - sum(sum(X==0))/20000 - (j-1)*5/100 > 0.005)
                    taun = taun + 0.001;
                elseif (sum(sum(Xback==0))/20000 - sum(sum(X==0))/20000 - (j-1)*5/100 < -0.005)
                    taun = taun - 0.001;
                else
                    flag = false;
                end
                taun
            end
            X(logical(Z)) = 0;
            tauest(j,i) = taun;
            flag = true;
        end
        csvwrite(sprintf('percent_%.2f_%d.csv',(j-1)*5/100,i),X);
    end
end

% ZIPFA
u = [[0*ones(35,1); 2*ones(45,1); 1.7*ones(60,1);0*ones(60,1)]...
    [1.8*ones(35,1);0.9*ones(45,1); 0*ones(120,1)]...
    [0*ones(35,1);1.7*ones(165,1)]];
vt = [[0.0*ones(1,30) 0*ones(1,30) 1.7*ones(1,40)];...
    [0*ones(1,35) 1.7*ones(1,25) 1*ones(1,40)];...
    [1.7*ones(1,25) 0.9*ones(1,50) 0.9*ones(1,25)]];
u=normrnd(0,0.06,200,3)+u;
vt=normrnd(0,0.05,3,100)+vt;
a = u * vt;
for j = 1:15
    for i=1:200
        X(:,:,i)=csvread(sprintf('percent_%.2f_%d.csv',(j-1)*5/100,i));
    end
    parfor i = 1:200
        [Ufit,Vfit,itr,esttau,Likelihood,MLikelihood]=ZIPFA(X(:,:,i), 3, 'Madj', false);
        l2(j,i) =  sum(mean(abs(Ufit{itr}*Vfit{itr}' - a)).^2);
    end
end

% log-PCA
for j=1:15
    for i=1:200
        X(:,:,i)=csvread(sprintf('percent_%.2f_%d.csv',(j-1)*5/100,i));
    end
    X(X==0)=0.5;
    parfor i=1:200
        [coeff,score]=pca(log(X(:,:,i)),'centered',false,'NumComponents',3);
        fita=score*coeff';
        l2p(j,i) =  sum(mean(abs(fita - a)).^2);
    end
end

% Calculate GOMMS and PSVDOS in 'Simu.R'
% GOMMS
for j=1:15
    for i=1:200
        X(:,:,i)=csvread(sprintf('percent_%.2f_%d_fit_gomms_noomit.csv',(j-1)*5/100,i));
        l2g(j,i)=sum(mean(abs(X(:,:,i)-a)).^2);
    end
end
l2g(l2g==0)=nan;
l2g(l2g==l2g(2,1))=nan;
mean(l2g,2,'omitnan')

% PSVDOS
for j=1:15
    for i=1:200
        X(:,:,i)=csvread(sprintf('percent_%.2f_%d_fit_psvdos.csv',(j-1)*5/100,i));
        l2p(j,i)=sum(mean(abs(X(:,:,i)-a)).^2);
    end
end
mean(l2p,2,'omitnan')

% l2back=l2;
% l2back(l2back>800)=nan;
% mean(l2back,2,'omit')

%% Figure 3 Draw the figure
color = [[30 129 196]/255;
    [237 77 36]/255;
    [0,175,79]/255;
    [237 177 36]/255];
x1=[    2.35	2.85	2.15;
    4.240921658	7.86437818	4.438748386;
    4.54043176	13.65672207	5.377224722;
    5.283794914	23.89194611	7.121886851;
    4.53	32.84	7.07;
    6.13096393	61.76567724	14.49894369;
    7.021124924	93.01124534	20.56589872;
    9.023107884	135.25785	28.93038474;
    12.44	171.16	33.93;
    29.1848426	259.5316384	54.28545613;
    48.63327176	332.6732977	73.33458862;
    82.83260885	419.0385476	98.69515293;
    98.53400999	512.4159873	133.6693988;
    107.5079013	616.3426741	183.7525178;
    123.8794309	734.1625621	339.2184045;
    ];
e1=[    0.235203892	0.302464531	0.15;
    0.119432068	0.504770794	0.249428112;
    0.133193154	0.770730042	0.301974446;
    0.202303218	1.158752598	0.435176358;
    0.347342941	1.651221945	1;
    0.339010989	2.136118852	1.462801943;
    0.325386831	2.558363922	1.7084517;
    0.501242975	3.112059861	2.045920038;
    1.510424928	4.838159131	5;
    3.610570779	4.53224287	8.091334857;
    4.505783783	4.558781522	5.452908368;
    15.85675399	5.390863267	6.512214609;
    34.8675838	5.267888533	9.83833453;
    44.51777027	5.976695836	17.08385798;
    46.55567327	6.527096455	57.95944443;
    ];
x2=[4.47 5 6.2 12 17.21 30 32.45199774 45 48.97724955 80 97.06368667 172.7923586];
e2=[0.39 0.5 0.67 0.7 0.82 2 5.50709589 18 9.77480368 45 50.8540091 64.67320132 ];
h1=errorbar(x1,e1);
hold on
h2=errorbar(x2,e2);
set(h1,'linew',1.5)
set(h2,'linew',1.5)
set(h1(1),'color',color(1,:))
set(h1(2),'color',color(2,:))
set(h1(3),'color',color(3,:))
set(h2,'color',color(4,:))
ylim([0,410])
xlim([0 16])
ax=gca;
ax.XTick=[1:2:15];
ax.XTickLabel={'0%' '10%' '20%' '30%' '40%' '50%' '60%' '70%'};
legend([h1(2) h2 h1(3) h1(1)],{'log-SVD' 'GOMMS' 'PSVDOS' 'ZIPFA'},'location','northwest');
xlabel('Inflated zero percentage')
ylabel('L_2 norm')

%% Figure 2
res = [];
for g=1:11
    for i = 1:100
        [cvsample,allres] = cv_ZIPFA(csvread(sprintf('g%d_%d.csv',g,i)), 1:5, 'Madj', false);
        res = [res; allres];
    end
end

color = [[30 129 196]/255; % b
    [237 77 36]/255;  %r
    [0,175,79]/255;   %g
    [237 177 36]/255;  %y
    [182 122 199]/255;   %p
    [22, 171, 192]/255;]; % lb

h1=errorbar(m(:,[1,2,4,6,8,10]),s(:,[1,2,4,6,8,10]))
set(h1,'linew',1)
for i=1:6
    set(h1(i),'color',color(7-i,:))
end
hold on
h2=errorbar(m(:,[3:2:11]),s(:,[3:2:11]),'--')
set(h2,'linew',1)
for i=2:6
    set(h2(i-1),'color',color(7-i,:))
end
axis([1.5,5.2,-7200,-4500])
legend([h1(1:2) h2(1) h1(3) h2(2) h1(4) h2(3) h1(5) h2(4) h1(6) h2(5)],...
    {'0% Inflated Zero','Setting 1 (20%)','Setting 1 (40%)','Setting 2 (20%)','Setting 2 (40%)',...
    'Setting 3 (20%)','Setting 3 (40%)','Setting 4 (20%)','Setting 4 (40%)','Setting 5 (20%)','Setting 5 (40%)'},...
    'location','southeast')
g=gca;
g.XTick=[2:5];
xlabel('Rank')
ylabel('CV likelihood')
pnew

%% Web Figure 2
load('uvax')
% fitted: the matrix to draw.
% fitted=readtable('examplefitted.csv');
% fitted=fitted{:,2:end};

subplot(2,2,4)
hx=heatmap(fitted,[],[],[],'Colormap',hot1,'MinColorValue',-0.2,'MaxColorValue',6.6);
hxax = get(hx,'parent');
hxax.OuterPosition=[0.1111   -0.0959    0.9665  0.9490];
h1=axes(gcf);
Z = linkage(fitted,'complete');
d=dendrogram([Z(:,1:2) log(log(Z(:,3)/min(Z(:,3))/0.95*exp(1)))],200,...
    'reorder',200:-1:1,'ColorThreshold',1.86,'checkcrossing',false,'Orientation','left');
h1.OuterPosition=[ 0.010589568665245  -0.101877672429820   0.347221064488924   0.987378901489331];
h2=axes(gcf);
Z = linkage((fitted)','complete');
d=dendrogram([Z(:,1:2) log(log(Z(:,3)/min(Z(:,3))/0.95*exp(1)))],200,...
    'reorder',1:100,'ColorThreshold',1.7,'checkcrossing',false);
h2.OuterPosition=[  0.058009608785175   0.620218578467571   1.026767330130405   0.358490060776120];
set(gcf,'children',[hxax h1 h2])
httx=axes(gcf)
htt=histogram(httx,fitted,60,'FaceAlpha',0,'EdgeColor','blue','LineWidth',0.7)
httx.OuterPosition=[0.060141843971631  -0.101860605451071   1.021276595744681   0.985618022729559];
httx.Visible='off';
httx.XLim=[-0.2,6.6];
httx.YLim=[0,2050];
h1.Visible='off';
h2.Visible='off';
httx.Visible='off';
set(gcf,'Position',1.0e+02 *[4.882 3.162  5.08 4.456])

l2=sum(mean(abs(fitted-a)).^2);
lsum=sum(sum((fitted-a).^2));
dif=sum((fitted-a).^2);
simp=sum((dif./lsum).^2);
bc=sum(mean(abs(fitted-a)))/2;
annotation('textbox',[ 0.639297 0.604835  0.339007  0.186482],...
    'String',{['L2 Norm: ' num2str(l2,'%.3f')],['Simpson: ' num2str(simp,'%.4f')],['Bray-Curtis: ' num2str(bc,'%.2f')]},...
    'FitBoxToText','on','edgecolor','none','Fontsize',16,'Facealpha',0.65,'Backgroundcolor',[1 1 1]);
pnew

%% Web Table 1
% in 'Simu.R'