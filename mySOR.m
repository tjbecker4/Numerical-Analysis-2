function [x,iter,info] = mySOR(A,b,omega,maxit,tol)
% Tim Becker
% 1/24/13
% CAAM 554 Homework 2, Problem 1.
%
% [x,iter,info] = mySOR(A,b,omega,maxit,tol)
%
% mySOR.m uses the SOR method to solve Ax = b. 
%
% Input Commands
% A          a user supplied matrix.
%
% b          a user supplied vector.
%
% maxit      a positive integer specifying the max number
%            of iterations allowed.
%
% tol        a positive real number (the stopping tolerance).
%
% Output Commands      
% x         The approximate solution to Ax = b.
%
% iter       A vector containing residual norms.
%
% info      A string detailing when the scheme converged, if it did so.


x = zeros(size(b));

Q = diag(diag(A))/omega + tril(A,-1);
Z = A-Q;
reltol = tol*(1 + norm(b));
v1 = b;
res(1) = norm(v1);

if res(1) <= reltol
    info = 'Converged at Iteration 0';
end

for iter = 2:maxit
    x = Q\v1;
    v2 = v1;
    v1 = b - (Z)*x;
    res(iter) = norm(v1-v2);
    if res(iter) <= reltol
       info = ['Converged at Iteration ' num2str(iter)];
       break
    end
end

if i == maxit
    info = 'Maxit reached';
end
    