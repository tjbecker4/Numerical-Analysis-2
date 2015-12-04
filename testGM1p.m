% testGM1p: GMRES with pre-conditioner

close; clear

%n = 500; 
n = input('n = ');
A0 = spdiags((n-1:-1:1)',-1,n,n);
A0(:,n) = (n:-1:1)'; 
e = ones(n,1);
b = zeros(n,1);
b(1) = 1;

fullscreen = get(0,'ScreenSize');
figure('Position',[0 -50 .8*fullscreen(3) .8*fullscreen(4)])

for s = 1:4 % s-loop
    
fprintf('\n\t ------- s = %g -------\n',s);
D = spdiags(e,s+2,n,n); %default is s+2
A = A0 + D - D';

% pre-conditioners
dtol = 1e-3;
[L,U] = ilu(A,struct('type','ilutp','droptol',dtol));
fprintf('\tnnz(A) = %i,   nnz(L+U) = %i\n',nnz(A),nnz(L)+nnz(U))

% plot eigenvalues
eval(['subplot(24' int2str(s) ');']);
eA = eig(full(A)); 
plot(real(eA),imag(eA),'b.'); axis square;
title('eig(A)','FontSize',14);
eval(['subplot(24' int2str(s+4) ');']);
eB = eig(full(U\(L\A))); 
plot(real(eB),imag(eB),'r*'); axis square;
title('eig(inv(M)*A)','FontSize',14);

bnrm = norm(b); tol = 1e-12;
if exist('myGMRES1p','file') && 1;
    fprintf('myGMRES1p:\n');
    t0 = cputime;
    [x,iter] = myGMRES1p(A,b,tol);
    cput = cputime - t0;
    rres = norm(A*x - b);
    fprintf('          Without pre-conditioner\n');
    fprintf(' Iter %4i: rel_res = %6.2e  cpu = %6.2e\n',...
        iter,rres,cput);
    t0 = cputime;
    [x,iter] = myGMRES1p(A,b,tol,L,U);
    cput = cputime - t0;
    rres = norm(A*x - b);
    fprintf('          With    pre-conditioner\n');
    fprintf(' Iter %4i: rel_res = %6.2e  cpu = %6.2e\n',...
        iter,rres,cput);
end

if exist('yzGMRES1p','file');
    fprintf('yzGMRES1p:\n');
    t0 = cputime;
    [x,iter] = yzGMRES1p(A,b,tol);
    cput = cputime - t0;
    rres = norm(A*x - b);
    fprintf('          Without pre-conditioner\n');
    fprintf(' Iter %4i: rel_res = %6.2e  cpu = %6.2e\n',...
        iter,rres,cput);
    t0 = cputime;
    [x,iter] = yzGMRES1p(A,b,tol,L,U);
    cput = cputime - t0;
    rres = norm(A*x - b);
    fprintf('          With    pre-conditioner\n');
    fprintf(' Iter %4i: rel_res = %6.2e  cpu = %6.2e\n',...
        iter,rres,cput);
end

end % s-loop