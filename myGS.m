function [x,info,res] = myGS(A,b,maxit,tol)
% Tim Becker
% 1/17/13
% CAAM 554 Homework 1, Problem 1.
%
% [x,info,res] = myGS(A,b,maxit,tol)
%
% myGS.m uses the Gauss-Seidel method to solve Ax = b. 
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
% info      A string detailing when the scheme converged, if it did so.
%
% res       A vector containing residual norms.

x = zeros(size(b));

Q = tril(A);
U = triu(A,1);
reltol = 1 + norm(b);
r = b - U*x-Q*x;
res(1) = norm(r);

if res(1) <= tol*reltol
    info = 'Converged at Iteration 0';
end

for i = 1:maxit
    x = x + Q\r;
    res(i+1) = norm(r);
    r = b - U*x-Q*x;
    if res(i+1) <= tol*reltol
       info = ['Converged at Iteration ' num2str(i)];
       break
    end
end

if i == maxit
    info = 'Maxit reached';
end
    