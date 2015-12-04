function [x,iter] = myGMRES2(A,b,p,tol,maxit)

%
% Tim Becker
% CAAM 554 Homework 3 Problem 1c
% 1/31/2012
%
% myGMRES2.m solves Ax=b using the restarted GMRES version that was posited in class
%
% [x,iter] = myGMRES1(A,b,p,tol,maxit)
%
%
% Input:
%
% A         Input matrix
% b         Input vector
% p         Amount of iterations before restarting
% tol       Stopping tolerance
% maxit     Maximum number of iterations
%
% Output:
%
% x         Solution to Ax=b
% iter      Number of iterations, both inner and outer

% Variable Checking
if ~exist('maxit','var')
    maxit = numel(b);
end

% Initialization

x = zeros(numel(b),1);
reltol = tol*(1+norm(b));
r = b;

for n = 1:maxit
     
   [dx,inner] = myGMRES1(A,r,tol,p,norm(b));
   
   x = x + dx;
   
   r = b-A*x;
       
   if norm(r) <= reltol,
       break
   end
    
end

iter(1) = n;
iter(2) = inner-1;



