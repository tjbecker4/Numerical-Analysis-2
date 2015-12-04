%clear; close all

M0 = delsq(numgrid('S',102));
n = size(M0,1);

tol = 2e-2;
maxit = 20;
fprintf('n = %i, tol = %.e\n',n,tol)

d = 1./randperm(n)'; d = (d).^.5;
A = M0*spdiags(d,0,n,n)*M0';
nrmA = norm(A,'fro');
A = A / nrmA;
r = @(U,D) max(sqrt(sum((A*U-U*D).^2)));

Ks = input('k = ');
if isempty(Ks); Ks = 100; end
nk = numel(Ks);
t1 = zeros(nk,1);
t2 = zeros(nk,1);

for j = 1:nk
    
k = Ks(j);
fprintf('\nk = %i\n',k)

% EIGS
opts.tol = tol;
tic; [U1,D1] = eigs(A,k,'LA',opts); 
[d1,p] = sort(diag(D1),'descend');
U1 = U1(:,p); D1 = diag(d1);
%X1 = U1*sqrt(D1);
t1(j) = toc; 
fprintf('EIGs ... '); toc

% GN4real
kw = round(1.1*k); quiet = false;
% edit the following line to run your code
tic; [X2,iter2] = myGN4real(A,kw,tol,maxit,quiet);
[U2,~] = qr(X2,0); d20 = d(1:k);
[V,D2] = eig(U2'*A*U2);
[d2,p] = sort(diag(D2),'descend');
U2 = U2*V(:,p); d2 = d2(1:k);
U2 = U2(:,1:k); D2 = diag(d2);
t2(j) = toc; 
fprintf('GN   ... '); toc

if 0 < 1
    fprintf('\nEIGs    : residual = %.4e, trace = %.8e\n',...
        r(U1,D1),sum(d1));
    fprintf('iter %3i: residual = %.4e, trace = %.8e\n',iter2,...
        r(U2,D2),sum(d2));
end

end

h = plot(Ks,t1,'ro:',Ks,t2,'bs--');
set(h,'linewidth',2,'markersize',6)
legend('EIGS','GN','Location','best')
title(sprintf('Solution Time (n = %i)',n),'fontsize',16), 
xlabel('Eigenspace Dimension k','fontsize',16)
set(gca,'fontsize',16)
grid on; shg