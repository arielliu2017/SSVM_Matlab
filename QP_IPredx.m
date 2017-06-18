%solves QP problem 
%                       min x'Qx/2+c'x
%       subject to      Ax>=b,   x>=0

function [IP] = QP_IPredx(Q, c, A, b, TOL)

MAX_ITERS = 200; 


N=size(c,1);
M = size(b,1); 
x = ones(N,1);
z = x;
y = ones(M,1);
w = y;

E0 = [x; w; y; z];%initial point
e2 = y;
e3 = z;
delta = 1/10;

 
for i = 1:MAX_ITERS%i is the iteration of the IP algorithm
     i
     c1 = b-A*x;
     rho = c1+w;
     c2 = -A'*y+c+Q*x
     sig = c2-z;
     gamma = z'*x+y'*w; 
     mu = delta*gamma/(M+N);

%     X = diag(x);
%     Z = diag(z);
%     W = diag(w);
%     Y = diag(y);
     
%     XI = pinv(X);
%     YI = pinv(Y);
     
     xi = 1./x;
     yi = 1./y;
     wi = 1./w;
     %K = -(XI*Z+Q);
     %R = [K A'; A YI*W];
     %r = [sig+z-mu*XI*e3; rho-w+mu*YI*e2];
     
     AtYWI = A'*diag(y.*wi);
     r = c2 - mu*xi - AtYWI*(c1+mu*yi);
     R = -diag(xi.*z)-Q-AtYWI*A;
     
     %%%[U,S,V] = svd(R,0);
     %%%DExy = V*pinv(S)*U'*r;
     
     dx = R\r;
     dy = diag(y.*wi)*(c1+mu*yi-A*dx);
     
     
     %dx = DExy(1:N);
     %dy = DExy(N+1:N+M);
     
  %   dz = XI*(mu*e3-X*Z*e3-Z*dx);
  %   dw = YI*(mu*e2-Y*W*e2-W*dy);
     dz = mu*xi - z - (xi.*z).*dx;
     dw = mu*yi - w - (yi.*w).*dy;


      mm = max([-dx./x; -dy./y; -dz./z; -dw./w]);
      theta(i)   =  0.9/mm;

     if theta(i) < 10^(-9)%if the step size is too small stop
         i
         break;
     end
     
     x = x + theta(i)*dx;
     y = y + theta(i)*dy;
     z = z + theta(i)*dz;
     w = w + theta(i)*dw;
    
    IP.xx(:,i) = x;
    IP.ww(:,i) = w;
    IP.pp(:,i) = y;
    IP.zz(:,i) = z;
     
     
     F = [sig; rho; -x.*z+mu*e3; -y.*w+mu*e2];

     met1(i) = norm(rho);%measure of primal constraint
     met2(i) = norm(sig);%measure of dual constraint
     met3(i) = norm(gamma);%measure of complementarity
     met4(i) = norm(F);%value of F that should be zero for a solution
     
     err = max([met1(i) met2(i) met3(i) met4(i)])
     if err < TOL
         exitflag = 1;
         val = x'*Q*x/2 + c'*x;
         Total_Iters = i
          IP.err = err
         break
     else
         exitflag = 0;%did not converge in maxiters
         val = -999;%return a nonsense value
          IP.err = err
                   Total_Iters = i

     end
end


IP.x = x;
IP.w = w;
IP.p = y;
IP.z =z;
IP.val = val;
IP.iters = Total_Iters;
IP.exitflag = exitflag;
IP.conerr = A*x-b;
IP.TOL = TOL;
IP.MAX_ITERS = MAX_ITERS;

