function [f,g] = mylork(x,n,k,A,nrmA2)

%
% Tim Becker
% CAAM 554 Homework 8 Problem 2
% 3/21/2012
%
% mylork.m returns the function value and the gradient of our given function, in vector form.
%
% [f,g] = mylork(x,n,k,A,nrmA2)
%
%
%
% Input:
%
% x        	Vectorized X
% n,k       Dimensions
% A         Input Matrix
% nrmA2     Frobenius norm of A, squared
%
% Output:
%
% f         vectorized function value
% g         vectorized gradient value


X = reshape(x,n,k);

f = ((1/4)*(norm(X*X'-A,'fro'))^2)/nrmA2;

f = f(:);

if nargout > 1
    g = (X*(X'*X)-A*X)/nrmA2;
    g = g(:);
end



% This code will calculate the function value and the gradient value of our
% given function. We use nargout so that we do not calculate the gradient
% if we do not need it. We reshape the vector into a matrix, and then deal
% with matrices, before turning the function and gradient back into vectors
% to output. My code is slower than the professor's, but this could be
% based on a variety of factors, not just in this code. My BFGS might be
% slower based on the rho, c, or alpha values. Additionally, for
% testBFGS2, we have that it does not converge for our first test instance.
% This is due to the fact that our alpha values go down to machine
% precision very quickly, so our step size is very small each time. This
% could be down to the treatment of the initial guess, which is the
% identity matrix in n by k dimensions (then vectorized). 

