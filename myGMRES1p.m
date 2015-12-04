function [x,iter] = myGMRES1p(A,b,tol,M1,M2,maxit)

%
% Tim Becker
% CAAM 554 Homework 3 Problem 1b
% 1/31/2012
%
% myGMRES1p.m solves Ax=b using the pre-conditioned GMRES version that was posited in class
%
% [x,iter] = myGMRES1p(A,b,tol,M1,M2,maxit)
%
%
% Input:
%
% A         Input matrix
% b         Input vector
% tol       Stopping tolerance
% M1        Part of pre-conditioning matrix M
% M2        Second part of M
% maxit     Maximum number of iterations
%
% Output:
%
% x         Solution to Ax=b
% iter      Number of iterations

% Variable Checking

if ~exist('M1','var') && ~exist('M2','var')
    [x,iter] = mygmres1(A,b,tol);
    return;
end

if ~exist('maxit','var')
    maxit = numel(b);
end

% Initialization

b = M2\(M1\b);
r = norm(b);
Q(:,1) = b/r;
w = 1;
n = 1;
reltol = tol*(1+r);


while norm(r) > reltol && n <= maxit
    % ARNOLDI
    v = M2\(M1\(A*Q(:,n)));
    H(:,n) = Q'*v;
    v = v - Q*H(:,n);
    H(n+1,n) = norm(v);
    Q(:,n+1) = v/H(n+1,n);
    
    f = H(n+1,n);
    d = conj(w)*H(1:n,n);
    dd = sqrt(f^2+d^2);
    
    
    nu = d/dd;
    t = f/dd;
    
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

