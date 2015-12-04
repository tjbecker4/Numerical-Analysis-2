function [x,iter,res,actual] = myCG(A,b,tol)

%
%
%
%
%
%
%
%
%
%
%


x = 0;
res(1) = norm(b);
p = res(1);
iter = 0;
reltol = tol*(1+res(1));
actual(1) = res(1);

while iter <= numel(b) && res(iter+1) > reltol
    
    r1 = res(iter+1)'*res(iter+1);
    P = A*p;
    
    alpha = r1/(p'*P);
    x = x + alpha*p;
    res(iter+2) = res(iter+1) - alpha*P;
    beta = (r'*r)/r1;
    p = res(iter+2) + beta*p;
    
    iter = iter + 1;
    actual(iter+2) = norm(b-A*x);
    
end

