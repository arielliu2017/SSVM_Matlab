%This code solves an Linear Program using a primal dual interior 
%point method of the form
%       max c'x
%sub.       Ax=<b, x>=0.
%note the different standard form


function [IP] = LP_IP_adaptxy(c, A, b, TOL)

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
delta = 1/10
err = 999;


for i = 1:MAX_ITERS%i is the iteration of the IP algorithm
    i
    
    b1 = b-A*x;
    b2 = c-A'*p;
    
    rho = -w+b1;
    sig = z+b2;
    
    gamma = z'*x+p'*w; 
    mu = delta*gamma/(M+N);
%    diag(x)*diag(z)*e3 can be replaced by
    a1 = x.*z;
    a2 = p.*w;
    
    F = [rho; sig; a1-mu*e3; a2-mu*e2];
    
  a3 = 1./p;
  a4 = 1./x; 

if err > 1%solve for dx
    disp('solving dx')
    a6 = 1./w;
    AtWIY = A'*diag(a6.*p);
    R = AtWIY*A + diag(a4.*z);
%    condRdx(i) = cond(R)
    t2 = b1-mu*a3;
    r = b2+ mu*a4 + AtWIY*t2;
    dx = R\r;
    dp = -diag(a6.*p)*(t2-A*dx);
else
   disp('solving dy')
   a5 = 1./z;
   AXZI = A*diag(x.*a5);
   t1 = b2+mu*a4;
   r = b1-mu*a3-AXZI*t1; 
   R = -diag(a3.*w)-AXZI*A';
%   condRdy(i) = cond(R)
   dp = R\r;
   dx = diag(x.*a5)*(t1-A'*dp);
end    



dz = mu*a4 - z - (a4.*z).*dx;
dw = mu*a3 - w - (a3.*w).*dp; 






mm = max([-dx./x; -dp./p; -dz./z; -dw./w]);
newtheta = 0.9/mm
theta(i)   =  newtheta;

     if theta(i) < 10^(-9)%if the step size is too small stop
         i
         break;
     end
     
     x = x + theta(i)*dx;
     p = p + theta(i)*dp;
     z = z + theta(i)*dz;
     w = w + theta(i)*dw;



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
   
    
    
    errb = max([IP.bmet1(i) IP.bmet2(i) IP.bmet3(i) IP.bmet4(i)])
    erra = max([IP.amet1(i) IP.amet2(i) IP.amet3(i) IP.amet4(i)])
    err = max([IP.met1(i) IP.met2(i) IP.met3(i) IP.met4(i)])
    
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

