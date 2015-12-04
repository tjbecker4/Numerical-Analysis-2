function [x,iter,nf] = myCGD(func,method,x0,tol,maxit,varargin)

%
% Tim Becker
% CAAM 554 Homework 6 Problem 1
% 3/7/2012
%
% myCGD.m uses a variety of CG methods to minimize the given function.
%
% [x,iter,nf] = myCGD(func,method,x0,tol,maxit,varargin)
%
%
%
% Input:
%
% func      Input function
% method    Integer detailing which gradient method should be used
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


% Initialization
alpha = 1;
c = 10^(-4);
rho = .1;
iter = 0;
x = x0;
[f,g] = feval(func,x,varargin{:});
p = -g;

nf = 1;
crit = norm(g)/(1 + abs(f));


while crit >= tol && iter < maxit
    
    if mod(iter,50) == 0
      fprintf('iter:   %3i f = %6.4i, norm(g) = %6.2e, crit = %6.2e\n',iter,f,norm(g),crit)
    end
    
   
    x1 = x+alpha*p;
    [fk,gk] = feval(func,x1,varargin{:});
    nf = nf + 1;
    
    % Backtracking
    while fk > f+c*alpha*g'*p
        
        alpha = rho*alpha;
        x1 = x + alpha*p;
        [fk,gk] = feval(func,x1,varargin{:});
        nf = nf + 1;
        
    end
    x1 = x + alpha*p;

    % BB Step
    s = x1 - x;
    y = gk - g;
    alpha = (s'*s)/abs(s'*y);
    
   
    
    % Check method
    if method == 1
        beta = 0;
    elseif method == 2
        beta = gk'*(gk-g)/(g'*g);
    elseif method == 3
        beta = max(0,gk'*(gk-g)/(g'*g));
    elseif method == 4
        beta = (gk'*gk)/(g'*g);
    end
  
   % Check if descent direction
    if gk'*(-gk+beta*p) >= 0 || ~exist('method','var')
        beta = 0;
    end
    
    f = fk;
    g = gk;
    x = x1;
    
    p = -g + beta*p;
    
    iter = iter + 1;
    crit = norm(g)/(1 + abs(f));

end

if mod(iter,50) ~= 0
     fprintf('iter:   %3i f = %6.4i, norm(g) = %6.2e, crit = %6.2e\n',iter,f,norm(g),norm(g)/(1 + abs(f)))
end
    

% We notice in problem 1 that myCGD.m uses many fewer function evaluations
% than the given yzCGD.m. In both Deblur and Denoise, myCGD is faster due
% to this. The use of line-search and the BB step in the code acts to help
% expedite the convergence. Both methods converge to the same f value, so
% the differences come from the cpu time. The images in Deblur and Denoise
% are quite similar and cannot be distinguished from each other by the
% naked eye. This is due to their similarities in relative error.
% (Discussion of the methods is on the myNewton.m page with other
% discussion from question 2)

