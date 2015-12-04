%close all

% test panda1 or panda2
load Panda1;
%load Panda2;

M0 = Panda;
[m,n] = size(M0); 
nrm0 = norm(M0,'fro')^2;
f = @(M).25*norm(M-M0,'fro')^2/nrm0;

% target rank: e.g., 32 to 128
k = input('k = ');
tol = 1e-2;
maxit = 100;
fprintf('m = %i, n = %i, k = %i, tol = %.e\n\n',...
    m,n,k,tol)

A = M0*M0';
nrmA = normest(A,1e-3);
A = A / nrmA;
opts.tol = tol;

% Dimension reduction by EIGS
tic;
[U1,~] = eigs(A,k,'LM',opts); 
V1 = (U1'*M0)';
M1 = U1*V1';
fprintf(' EIGs  ... '); toc

% Dimension reduction by yzGNes
tic; [X,iter2] = yzGNes(A,k,tol,maxit);
[U2,~] = qr(X,0); V2 = (U2'*M0)'; 
M2 = U2*V2';
fprintf('yzGNes ... '); toc

subplot(221); imshow(M0,[]); 
xlabel(['Original Rank ' num2str(min(m,n))])
subplot(222); imshow(M1,[]);
xlabel(['SVDs Rank ' num2str(k)])
subplot(223); imshow(M2,[]);
xlabel(['yzGNes Rank ' num2str(k)])
shg

% Dimension reduction by myGNes
if exist('myGNes','file')
    tic; [X,iter3] = myGNes(A,k,tol,maxit);
    [U3,~] = qr(X,0); V3 = (U3'*M0)';
    M3 = U3*V3';
    fprintf('myGNes ... '); toc
    subplot(224); imshow(M3,[]);
    xlabel(['myGNes Rank ' num2str(k)])
end

% more info
if 0 < -1
    fprintf('\nEIGs   : f = %.4e\n',f(M1));
    fprintf('iter %2i: f = %.4e\n',iter2,f(M2));
    if exist('M3','var')
        fprintf('iter %2i: f = %.4e\n\n',iter3,f(M3));
    end
end
