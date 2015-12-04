function [x,iter] = myCGN(A,b,tol)

%
% Tim Becker
% CAAM 554 Homework 4 Problem 1
% 2/7/2012
%
% myCGN.m solves Ax=b using the CG version that was posited in class and
% applied to the normal equations.
%
% [x,iter] = myGMRES1(A,b,tol)
%
%
% Input:
%
% A         Input matrix
% b         Input vector
% tol       Stopping tolerance
%
% Output:
%
% x         Solution to Ax=b
% iter      Number of iterations

% Initialization

b = A'*b;

x = 0;
r = b;
p = r;
iter = 0;
reltol = tol*(1+norm(b));
r1 = r'*r;

% CGN
while iter <= numel(b) && sqrt(r1) >= reltol
    
    P = A'*(A*p);
    
    alpha = r1/(p'*P);
    x = x + alpha*p;
    r = r - alpha*P;
    rr = r'*r;
    beta = rr/r1;
    p = r + beta*p;
    r1 = rr;
    
    iter = iter + 1;
    
end

% EXPLANATION
%
% myCGN.m converges in fewer iterations and a smaller relative residual
% than Matlab PCG since I used the relative tolerance of tol*(1+norm(b)),
% as opposed to tol*norm(b). As such, I generated output for each case to
% show this. My code is faster than Matlab's due to the optimized nature of
% the code that we were given in class, in which we do not store anything,
% as well as calculate the inner products only one time each. Lastly, for
% the n<=1000 cases, the other solvers failed. This is most likely down to
% the fact that the matrix A was non-symmetric. We used the normal
% equations so that the CG method would work, but the other methods did not
% work. I do not know the inner workings of these methods, but I imagine
% the failure was due to the lack of ability to deal with such a matrix.

