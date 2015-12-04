%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This is a front-end script for solving "standard" LP:
%
%                min   c'*x
%                st  A*x = b
%                     x >= 0
%
%  by primal-dual interior point method.  It needs
%  problem data A, b and c.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 5e-8; [m, n] = size(A);

tcpu = cputime;
% call your function pdipm.m
[x,y,z,iter] = pdipm(A,b,c,tol);
% call instructor's Zpdipm.p
%[x,y,z,iter] = Zpdipm(A,b,c,tol);
tcpu = cputime - tcpu; 

fprintf('\n  Prob. name = %s\n', name);
fprintf('  Prob. size [m, n] = [%g %g]\n',m,n);
fprintf('  No. of iterations used: %i\n',iter);
fprintf('  Primal obj. value = %8.4e\n',c'*x);
fprintf('  Total CPU = %g sec.\n', tcpu);

%  check feasibility and optimality
tol = 1e-7;
if abs(b'*y - c'*x)/(1 + abs(b'*y)) < tol && ...        % duality gap
   all(x > -tol) && norm(b-A*x)/(1+norm(b)) < tol && ...% primal feasible
   all(z > -tol) && norm(A'*y+z-c)/(1+norm(c)) < tol    % dual feasible
   fprintf('  Accuracy up to %g\n',tol);
else
   fprintf('  Accuracy less than %g\n',tol);
end
