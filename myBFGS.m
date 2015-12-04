function [x,iter,hist] = myBFGS(func,x0,tol,maxit,varargin)


%
% Tim Becker
% CAAM 554 Homework 8 Problem 1
% 3/21/2012
%
% Quasi-Newton method using inverse BFGS update 
% (if you do backtracking, make sure you check y'*s > 0)
%
% Usage: 
%   [x,iter,hist] = myBFGS(func,x0,tol,maxit,varargin);
%
% INPUT:
%       func - matlab function for f and g evaluations
%         x0 - starting point (column vector)
%        tol - stopping tolerance for gradient norm
%      maxit - maximum number of iterations
%
% OUTPUT:
%          x - computed approximate solution
%       iter - number of iterations taken
%      hist(1,:) - vector of function values at all iterations
%      hist(2,:) - vector of gradient norms  at all iterations

% Initialization
x = x0;
a0 = 1;
c = 10^(-4);
rho = 0.5;
n = numel(x);
hist = zeros(2,maxit);

[f,g] = feval(func,x,varargin{:});

p = -g;

crit = norm(g);

fprintf('  k       f_k       norm(g_k)      alpha_k\n');
fprintf('--------------------------------------------\n');



for iter=0:maxit
    
   if crit < tol*(1+abs(f))
       break
   end
   
   a = a0;
   x1 = x+a*p;
   [f1,~] = feval(func,x1,varargin{:});
    
   % Backtracking
   while f1 > f+c*a*g'*p    
       a = rho*a;
       x1 = x + a*p;
       [f1,~] = feval(func,x1,varargin{:});
   end
   
   x1 = x + a*p;
   [f1,g1] = feval(func,x1,varargin{:});
   
   s1 = a*p;
   y1 = g1 - g;
   
   % BB step for initial H
   if iter == 0
       H = ((s1'*s1)/(s1'*y1))*eye(n);
   end
    
   % Update H
   if s1'*y1 > 0
       eta = 1/(s1'*y1);
       Hy = H*y1;
       Hys = Hy*s1';
       H = H - eta*Hys-eta*Hys'+((eta^2*(y1'*Hy)+eta)*s1)*s1';
   end
   
   % Update
   g = g1;
   p = -H*g;
   x = x1;
   f = f1;
   
   
   hist(1:2,iter+1) = [f, norm(g)];
   crit = norm(g);
   

   % Print out values every 10 iterations
   if mod(iter,10) == 0
       fprintf(' %3i  %11.4e  %11.4e   %9.2e\n',iter,f,norm(g),a);
   end
    
end

% Print out final iteration
if iter ~= maxit
    fprintf(' %3i  %11.4e  %11.4e   %9.2e\n',iter,f,norm(g),a);
    fprintf('Converged!\n');
end



% In this problem, we notice that we achieve quadratic convergence after a
% certain number of iterations. This is what we expect, since we have a
% Quasi-Newton method. We have to begin with a BB step for our initial H,
% so we take our first step in the direction of the negative gradient and
% then find our first H. This way, we have a multiple of the identity to
% begin with as H and this helps the convergence of the method.
% Additionally, we increase the speed of the method by calculating our new
% H in a clever way.
