function X = mySolve4X(A,B,tol,maxit)

%
% Tim Becker
% CAAM 554 Homework 5 Problem 1
% 2/21/2012
%
% mysolve4X.m solves the simplified Sylvester Equation AX+XA=B for X 
% using a function handle and then Matlab's preconditioned CG method.
%
% X = mysolve4X(A,B,tol,maxit)
%
%
% Input:
%
% A         Input matrix
% B         Input matrix
% tol       Stopping tolerance
% maxit     Maximum number of iterations
%
% Output:
%
% X         Solution to AX + XA = B



n = size(A,1);

Afun = @(v) reshape((reshape(v,n,n))*A + A*(reshape(v,n,n)),numel(B),1);

b = reshape(B,numel(B),1);

[X,~] = pcg(Afun,b,tol,maxit);

X = reshape(X,n,n);


% Observations: 
% Using the reshape function allows the code to run much faster than the
% Kronecker product that we discussed in class. This took up far too much
% time and memory, even when we reshaphed the matrix into a vector to use
% just matrix-vector products. Thus, using a function handle for some
% vector v and then doing matrix-matrix multiplication after reshaping it,
% and then reshaping back into a vector allows us to use Pre-conditioned
% Conjugate Gradient to solve this generalized Sylvester Equation fairly
% quickly and with high accuracy. The plot show how close the error is to 0
% for the run with n=100, which corroborates how well the code does
% compared to the instructor's code. This code is slightly faster than the
% given program as well, in both the runs for n=100 and n=1000.






