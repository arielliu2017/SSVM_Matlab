function [IPl1] = SSVM_Train(X, labels, C, TOL, method)

%Data X (each point is a row).
%labels a column matrix with the labels of X
%C the margin weight
%TOL error tolerance for interior point method


P = size(X,1);
DIM = size(X,2);

%Variable break down
%x = [w+ w- g+ g- y]

Zdim = zeros(DIM,1);
Zdimp = zeros(DIM,P);
Idim = eye(DIM);
IP = eye(P);
eP = ones(1,P);
edim = ones(1,DIM);

%SET up matrices in A without forming D = diag(d)
M = repmat(labels, 1, size(X,2));
M1 = M.*X;
M2 = labels.*eP';

%A = [M1 -M1 -M2 M2 IP Zdimp' Zdimp'; R2; R3];
A = [M1 -M1 -M2 M2 IP];

b = [eP'];
c = [edim edim 0 0 C*eP]';

size(A);
size(b);
size(c);


method;

switch method
    case {'ReducedKKTxy'}
       disp('Method ReducedKKT xy')    
       IPl1 = LP_IPredxy(-c, -A, -b, TOL);    
    
    case 'KKT'
       disp('Using KKT code LP_IP')     
       IPl1 = LP_IP(c, A, b, TOL);
    
    case 'ReducedKKTx'
       disp('Method is ReducedKKT2')
       IPl1 = LP_IPredx(-c, -A, -b, TOL);    
    
    case 'ReducedKKTy'
       disp('Method is LP_IPredy_fast')
       IPl1 = LP_IPredy_fast(-c, -A, -b, TOL);    

    case 'adaptxy'
       disp('Method is adaptxy')
       IPl1 = LP_IP_adaptxy(-c, -A, -b, TOL);    
   
       
    otherwise
       disp('Unknown method.')
end
 

ww = IPl1.x(1:DIM)-IPl1.x(DIM+1:2*DIM);
gg = IPl1.x(2*DIM+1)- IPl1.x(2*DIM+2);

IPl1.wgt = ww;
IPl1.gamma = gg;
IPl1.C = C;


