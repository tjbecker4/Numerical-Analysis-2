function [x,iter] = myGMRES1(A,b,tol,maxit,relnorm)

%
% Tim Becker
% CAAM 554 Homework 3 Problem 1a
% 1/31/2012
%
% myGMRES1.m solves Ax=b using the GMRES version that was posited in class
%
% [x,iter] = myGMRES1(A,b,tol,maxit,relnorm)
%
%
% Input:
%
% A         Input matrix
% b         Input vector
% tol       Stopping tolerance
% maxit     Maximum number of iterations
% relnorm   Used when called from restarted GMRES to use correct norm
%           for residual
%
% Output:
%
% x         Solution to Ax=b
% iter      Number of iterations

% Variable checking

if ~exist('maxit','var')
    maxit = numel(b);
end

if exist('relnorm','var')
    reltol = tol*(1+relnorm);
else 
    reltol = tol*(1+norm(b));
end

% Initialization

r = norm(b);
Q(:,1) = b/r;
w = 1;
n = 1;


while norm(r) > reltol && n <= maxit
    % ARNOLDI
    v = A*Q(:,n);
    H(:,n) = Q'*v;
    v = v - Q*H(:,n);
    H(n+1,n) = norm(v);
    Q(:,n+1) = v/H(n+1,n);
    
    d = conj(w)*H(1:n,n);
    dd = sqrt((H(n+1,n))^2+d^2);
    
    
    nu = d/dd;
    t = H(n+1,n)/dd;
    w = -t*w;
    w(n+1) = conj(nu);
    r = t*norm(r);
    n = n + 1;
end

iter = n-1;

e = zeros(n,1);
e(1,1)=1;

y = H\e;
x = Q(:,1:n-1)*(norm(b)*y);

