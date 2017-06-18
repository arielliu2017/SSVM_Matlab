function S= get_accuracy(X, labels, IP)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

P = size(X,1);
e = ones(P,1);
c = sign(X*IP.wgt-IP.gamma*e);
a = c.*labels;

K = find(labels == 1);
numOfPos = length(K);
numOfNeg = P - length(K);

pos = 0;
neg = 0;
for i = 1:P
    if(a(i) > 0) && (labels(i) > 0)
        pos = pos + 1;
    elseif (a(i) > 0) && (labels(i)<0)
        neg = neg + 1;
    end
        
end
if(numOfPos == 0)
    S.acc = neg / numOfNeg;
elseif (numOfNeg == 0)
    S.acc = pos/numOfPos;
else
S.acc = 0.5*(pos/numOfPos) + 0.5*(neg / numOfNeg) ;
end
%S.pos = pos;
%S.neg = neg;



