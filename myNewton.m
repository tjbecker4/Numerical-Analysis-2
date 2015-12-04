function [x,iter,nf] = myNewton(func,x0,tol,maxit,varargin)

%
% Tim Becker
% CAAM 554 Homework 6 Problem 2a
% 3/7/2012
%
% myNewton.m used the Newton method to minimize the given convex function.
%
% [x,iter,nf] = myNewton(func,x0,tol,maxit,varargin)
%
%
% Input:
%
% func      Input function
% x0        Beginning vector 
% tol       Stopping tolerance
% maxit     Maximum number of iterations
% varargin  The variables that are taken in and used to evaluate the
%           function
%
% Output:
%
% x         Solution
% iter      Number of iterations
% nf        Number of function evaluations

iter = 0;
x = x0;


[f,g,H] = feval(func,x,varargin{:});
nf = 1;
crit = norm(g)/(1+abs(f));


while crit >= tol && iter < maxit
    
    fprintf('iter:   %3i f = %6.4i, norm(g) = %6.2e, crit = %6.2e\n',iter,f,norm(g),crit)
   
    
    x = x - (H\g);
    
    [f,g,H] = feval(func,x,varargin{:});
    nf = nf + 1;
    
    iter = iter + 1;
    crit = norm(g)/(1+abs(f));
    
    
end


% The Newton method minimizes the function without using line search, but
% converges fairly rapidly due to the quadratic convergence that takes over
% after a certain number of iterations. myNewton.m runs faster than the
% instructor's code due to the optimization of Aquartic.m. Additionally, we
% use H\g rather than computing the inverse, although the other code uses
% that as well (I assume!). For the CG methods, we see that the CG-FR
% method is the slowest of the four, and that the CG-PR+ method converges
% faster than the CG-PR method, which is expected. 