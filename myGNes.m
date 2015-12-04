function [X,iter] = myGNes(A,k,tol,maxit)


%
% Tim Becker
% CAAM 554 Homework 8 Problem 3
% 3/21/2012
%
% myGNes.m uses a variation on Gauss-Newton to solve the NLS problem
%
% [X,iter] = myGNes(A,k,tol,maxit)
%
%
%
% Input:
%
% A         Input Matrix
% k         Dimension of lower rank
% tol       Stopping tolerance
% maxit     Maximum number of iterations
%
% Output:
%
% X         Solution
% iter      Number of iterations


[n,~] = size(A);
x = rand(n,k);
x = x(:);
a = 1;
c = 10^(-4);
rho = 0.1;


[f,g,p] = mlork(x,n,k,A);

crit = norm(g)/(1+abs(f));

for iter =1:maxit
    if crit < tol
        break 
    end
    
    x1 = x - a*p;
%     f1 = mlork(x1,n,k,A);
%     
%     %Backtracking
%     while f1 > f-c*a*g'*p    
%         a = rho*a;
%         x1 = x - a*p;
%         f1 = mlork(x1,n,k,A);
%     end
   
    x = x1;
    
    [f,g,p] = mlork(x,n,k,A);

    crit = norm(g)/(1+abs(f));

end

X = reshape(x,n,k);
end

%
% mlork.m does the same as mylork.m, but also calculates the search
% direction, and assumes nrmA2=1.
%
% [f,g,p] = mlork(x,n,k,A)
%
%
%
% Input:
%
% x        	Vectorized X
% n,k       Dimensions
% A         Input Matrix
%
% Output:
%
% f         vectorized function value
% g         vectorized gradient value
% p         search direction

function [f,g,p] = mlork(x,n,k,A)


X = reshape(x,n,k);

f = ((1/4)*norm(X*X'-A,'fro')^2);

if nargout > 1
    
    AX = A*X;
    XX = X'*X;
    g = (X*(XX)-AX);
    g = g(:);
end

if nargout > 2

    p = X-AX/(XX);
    p = p(:);
end

end

% For this problem, we see that we are able to quickly recreate both panda
% pictures with a small proportion of the original information we were
% given. We see that EIGs is slower than our Gauss-Newton estimation. Here,
% I use the estimated Gauss-Newton method as well as a rewritten low rank
% code based on mylork.m. This one finds the search direction as well.
% Using the line search decreases the speed of this algorithm for our
% pandas. If we just used alpha=1 for the entire algorithm, then we would
% have very fast convergence. Additionally, it depends on our random matrix
% at the beginning. Sometimes there are instances in which this code and
% the professor's run much slower than EIGs, just based on the specific
% attributes of the random matrix.


