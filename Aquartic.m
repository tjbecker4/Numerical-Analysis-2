function [f,g,H] = Aquartic(x,A,a)

%
% Tim Becker
% CAAM 554 Homework 6 Problem 2
% 3/7/2012
%
% Aquartic.m finds the function, gradient, and Hessian of the quartic
% function that is given in the problem.
%
% [f,g,H] = Aquartic(x,A,a)
%
%
% Input:
%
% x         Input vector
% A         Input matrix
% a         Input scalar
%
% Output:
%
% f         Function
% g         Gradient
% H         Hessian

% Initialization
y = A*x;
z = x'*y;
n = numel(x);
e = ones(n,1);


% function
f = z^2/4 + a*(e'*x-n)^2/2;
   
if nargout > 1

    % gradient
    g = z.*y + a*(e'*x-n).*e;
    
end

if nargout > 2

    % Hessian
    H = z.*A + (2.*y)*y' + a.*(e*e');
    
end

