function [x,iter] = MyGN(func,x,tol,maxiter,prt,varargin)

%
% Tim Becker
% CAAM 554 Homework 7 Problem 9
% 3/14/2012
%
% MyGN.m uses the Gauss-Newton method to solve a set of nonlinear least
% squares problems
%
% [x,iter] = MyGN(func,x,tol,maxiter,prt,varargin)
%
%
%
% Input:
%
% func      Input function
% x         Beginning vector 
% tol       Stopping tolerance
% maxiter   Maximum number of iterations
% prt       switch for if printing info is on or off
% varargin  The variables that are taken in and used to evaluate the
%           function
%
% Output:
%
% x         Solution
% iter      Number of iterations


% Initialization
iter = 0;
alpha0 = .75;
c = 10^(-4);
rho = .11;

[r,J] = feval(func,x,varargin{:});
g = J'*r;
f = r'*r/2;
p = pinv(J)*r;


crit = norm(g);

while crit >= tol && iter < maxiter
    
    alpha = alpha0;
    
    [rk,Jk] = feval(func,x-alpha*p,varargin{:});
    fk = rk'*rk/2;
    
   
    % Line search
    while fk > f+c*alpha*g'*p
        
        alpha = rho*alpha;
        [rk,Jk] = feval(func,x-alpha*p,varargin{:});
        fk = rk'*rk/2;
       
    end
    
    % Update
    x = x - alpha*p;
    g = Jk'*rk;
    f = fk;  
    p = (Jk'*Jk)\g;
    
    iter = iter + 1;
    crit = norm(g); 
    
    % Print out values
    if prt                                     
         fprintf('iter: %2i  norm(F) = %7.3e\n',iter,crit);
    end
    
end



