%This code solves an Linear Program using a primal dual interior 
%point method of the form
%       min c'x
%sub.       Ax >=b, x>=0.

function [IP] = LP_IP(c, A, b, TOL)

MAX_ITERS = 200;


%This code implements an Interior Point Algorithm for Solving LPs.
%It is an implemtation of the Primal Dual Path Following ALgorithm derived
%in class.

N=size(c,1);
M = size(b,1);

x = ones(N,1);
z = x;
p = ones(M,1);
w = p;
E0 = [x; w; p; z];%initial point
e2 = p;
e3 = z;
delta = 1/10;



for i = 1:MAX_ITERS%i is the iteration of the IP algorithm
    i;
    rho = A*x-w-b;
    sig = A'*p+z-c;
    gamma = z'*x+p'*w; 
    mu = delta*gamma/(M+N);
    
    F = [rho; sig; diag(x)*diag(z)*e3-mu*e3; diag(p)*diag(w)*e2-mu*e2];

    DF1 = [A -eye(M) zeros(M,M) zeros(M,N)];
    DF2 = [zeros(N,N) zeros(N,M) A' eye(N)];
    DF3 = [diag(z) zeros(N,M) zeros(N,M) diag(x)];
    DF4 = [zeros(M,N) diag(p) diag(w) zeros(M,N)];

    DF = [DF1; DF2; DF3; DF4];
    DE = -DF\F;

    theta(i) = 0.9/(max(-DE./E0));

    if theta(i) < 10^(-6)%if the step size is too small stop
        i
        break;
    end
    
    E1 = E0 + theta(i)*DE;
    E1';
    E0= E1;

    x = E1(1:N,1);
    w = E1(N+1:N+M,1);
    p = E1(N+M+1:N+2*M,1);
    z = E1(N+2*M+1:2*N+2*M,1);

  
    IP.xx(:,i) = x;
    IP.ww(:,i) = w;
    IP.pp(:,i) = p;
    IP.zz(:,i) = z;
    
    IP.met1(i) = norm(rho);%measure of primal constraint
    IP.met2(i) = norm(sig);%measure of dual constraint
    IP.met3(i) = norm(gamma);%measure of complementarity
    IP.met4(i) = norm(F);%value of F that should be zero for a solution
   
    IP.amet1(i) = norm(rho,1);%measure of primal constraint
    IP.amet2(i) = norm(sig,1);%measure of dual constraint
    IP.amet3(i) = norm(gamma,1);%measure of complementarity
    IP.amet4(i) = norm(F,1);%value of F that should be zero for a solution
    
    IP.bmet1(i) = norm(rho,inf);%measure of primal constraint
    IP.bmet2(i) = norm(sig,inf);%measure of dual constraint
    IP.bmet3(i) = norm(gamma,inf);%measure of complementarity
    IP.bmet4(i) = norm(F,inf);%value of F that should be zero for a solution
   
    
    
    errb = max([IP.bmet1(i) IP.bmet2(i) IP.bmet3(i) IP.bmet4(i)]);
    erra = max([IP.amet1(i) IP.amet2(i) IP.amet3(i) IP.amet4(i)]);
    err = max([IP.met1(i) IP.met2(i) IP.met3(i) IP.met4(i)]);
    
    
    if err < TOL
        val = c'*x;
        Total_Iters = i;
        exitflag = 1;
        IP.err = err;
        break
    else
        exitflag = 0;
        val = -999;
        Total_Iters = i;
        IP.err = err;
    end
    
end


IP.x = x;
IP.w = w;
IP.p = p;
IP.z =z;
IP.val = val;
IP.iters = Total_Iters;
IP.exitflag = exitflag;
IP.conerr = A*x-b;
IP.TOL = TOL;
IP.MAX_ITERS = MAX_ITERS;

