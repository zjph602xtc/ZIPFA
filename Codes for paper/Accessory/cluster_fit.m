function [rate_taxa,rate_sub] = cluster_fit(X)
% function to calculate clustering accurarcy

z=linkage(X,'complete');
T = cluster(z,'maxclust',4);

ct = [1 2 3 4];
ct = perms(ct);
ct = [repmat(ct(:,1),1,35),repmat(ct(:,2),1,45),repmat(ct(:,3),1,60),repmat(ct(:,4),1,60)];
rate_taxa=max(sum(bsxfun(@eq,T',ct),2))/200;

z=linkage(X','complete');
T = cluster(z,'maxclust',4);

ct = [1 2 3 4];
ct = perms(ct);
ct = [repmat(ct(:,1),1,25),repmat(ct(:,2),1,10),repmat(ct(:,3),1,25),repmat(ct(:,4),1,40)];
rate_sub=max(sum(bsxfun(@eq,T',ct),2))/100;
end