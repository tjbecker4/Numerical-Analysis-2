% Homework 4 Problem 5 (Trefethen 38.6) Script
% Runs both myCG and mySD on a given matrix and plots their residuals, as
% well as a CG estimate on the same plot. The plot is output in the normal
% fashion, as well as a semilog plot.

% Initialization of matrix
d = 1:100;
A = diag(d) + diag(ones(99,1),-1) + diag(ones(99,1),1);
b = ones(100,1);
tol = 1.e-8;

% Calling each function
[x,iter,res,actual] = myCG(A,b,tol);

[y,iter1,residual] = mySD(A,b,tol);

n = 1:101;
Y = 2*((sqrt(cond(A))-1).^n)./((sqrt(cond(A))+1).^n);

% Plot 1
figure(1)
plot(n,res,'b');
hold on
plot(n,actual,'y');
hold on
plot(n,residual,'g');
hold on
plot(n,Y,'r');
legend('CG residual','CG actual', 'SD residual', 'CG estimate','location','northeast');
title('Residuals for SD and CG Plot');
xlabel('Iterations');
ylabel('Residual Values');

% Plot 2--semilog
figure(2)
semilogy(n,res,'b');
hold on
semilogy(n,actual,'y');
hold on
semilogy(n,residual,'g');
hold on
semilogy(n,Y,'r');
legend('CG residual','CG actual', 'SD residual', 'CG estimate','location','northeast');
title('Residuals for SD and CG Semilog Plot');
xlabel('Iterations');
ylabel('Residual Values');


% Comments
%
% We see that the actual residual as well as our calculated residual for
% the CG method are almost identical. The only difference can be seen at
% the end of the semilog plot. Additionally, both far outperform the SD
% method, which is what we expected to see. Additionally, our upper bound
% for our CG, which we called our CG estimate here, was only slightly
% better performing than the SD method, so our method converged much faster
% than this estimate.