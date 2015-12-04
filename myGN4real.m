function [X,iter] = myGN4real(A,k,tol,maxit,quiet)

%
% Tim Becker
% CAAM 554 Homework 9 Problem 1
% 4/2/2012
%
% myGN4real.m uses real Gauss Newton without line search to solve the NLS
% problem
%
% [X,iter] = myGN4real(A,k,tol,maxit)
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
X = rand(n,k);
X = X/(norm(X,'fro'));
nrmA = norm(A,'fro')^2;

for iter = 1:maxit
    XX = X'*X;
    Y = X/XX;
    AY = A*Y;
    P = AY - 0.5*X*(Y'*(AY)) - 0.5*X;
    X = X + P;
    f = objective(X,A,nrmA,XX);
    relchg = norm(P,'fro')/norm(X,'fro');
    
    if ~quiet 
        fprintf('iter %3i: f = %.4e, relcgh = %.2e\n',iter,f,relchg);
    end
    
    if relchg < tol
        break;
    end
    
end

end


function [f] = objective(X,A,nrmA,XX)

% objective.m gets the value of the objective function
%
% [f] = objective(X,A,nrmA,XX)
%
%
% Input:
%
% X         Input n by k matrix
% A         Input Matrix
% nrmA      Frobenius norm of A, squared.
% XX        X'*X
%
% Output:
%
% f         function value
%


f = 0.5*trace((XX)^2) - trace(X'*A*X) + 0.5*nrmA;


end




% My code ran faster than EIGS in general, but slower than the instructor's
% code. As k increased, the codes, predictably, ran slower, but the
% performance of my code in contrast with EIGS became stronger as k
% increased. The manner in which the objective value was evaluated, as well
% as the optimization of P and storing values gave my code its speed. The
% profiler showed that there was still the most time spent in computing P,
% which is to be expected even after optimization, but also that
% calculating XX takes up a good amount of the cpu time. The computational
% cost profile is attached at the back of this assignment. The function
% value and the relative change were fairly similar to the instructor's
% printouts, so it seems that the objective evaluations and the other
% computations are similar to the instructor's but his code did something
% tricky to make it even faster! One option that I could have done was
% change norm(X,'fro') to X(:)'*X(:), which might speed up slightly, but
% should not account for the differences in our codes. Regardless, EIGS is
% much slower, as can be seen in the output and the plot. Lastly, the
% difference in residuals for some of the runs is down to the randomness of
% the initial matrix, and the fact that our tolerance is just 2e-02. For 
% most of the runs, the residuals were around the same value, but a few had 
% quite different residuals than EIGS (both for mine and the instructor's 
% code).

