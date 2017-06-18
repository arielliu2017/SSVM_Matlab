function [Y,Xmean,Xscale] = Gene_Normalize(X)
% Pre-processing of gene features.
[n,p] = size(X);
Xmean = mean(X); X1 = X - repmat(Xmean,n,1);
Xscale = sqrt(sum(X1.*X1,1))/sqrt(n); Y = X1./repmat(Xscale,n,1); % Y = X1;
return