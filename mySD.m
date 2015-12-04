function [x,iter,res] = mySD(A,b,tol)

%
% Tim Becker
% CAAM 554 Homework 4 Problem 5
% 2/7/2012
%
% mySD.m solves Ax=b using the SD version that was posited in class and
% runs for numel(b) iterations
%
% [x,iter] = mySD(A,b,tol)
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


% Initialization
x = 0;
r = b;
p = r;
iter = 1;
reltol = tol*(1+norm(b));
r1 = r'*r;
res(1) = norm(b);

% SD
while iter <= numel(b) 
    
    P = A*p;
    
    alpha = r1/(p'*P);
    x = x + alpha*p;
    r = r - alpha*P;
    rr = r'*r;
    p = r;
    r1 = rr;
    
    % residual update
    res(iter+1) = sqrt(r1);
    
    iter = iter + 1;
    
end

