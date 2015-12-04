function [x,iter,res,actual] = myCG(A,b,tol)


%
% Tim Becker
% CAAM 554 Homework 4 Problem 5
% 2/7/2012
%
% myCG.m solves Ax=b using the CG version that was posited in class and
% runs for numel(b) iterations
%
% [x,iter] = myCG(A,b,tol)
%
%
% Input:
%
% A         Input matrix
% b         Input vector
% tol       Stopping tolerance
%
% Output:
%
% x         Solution to Ax=b
% iter      Number of iterations
% res       A vector of the residuals calculated
% actual    A vector of the norm of (b-Ax) at each iteration

% Initialization
x = 0;
r = b;
p = r;
iter = 1;
reltol = tol*(1+norm(b));
r1 = r'*r;
res(1) = norm(b);
actual(1) = res(1);

% CG
while iter <= numel(b)
    
    P = A*p;
    
    alpha = r1/(p'*P);
    x = x + alpha*p;
    r = r - alpha*P;
    rr = r'*r;
    beta = rr/r1;
    p = r + beta*p;
    r1 = rr;
    
    % residual updates
    res(iter+1) = sqrt(r1);  
    actual(iter+1) = norm(b-A*x);
    
    
    iter = iter + 1;
    
end

