% test myCGN and other matlab solvers

n = input('n = ');
if isempty(n); n = 500; end

t = 1-1.e-6; p = 10;
A = mat4gmres(n,'mat3',p,t);
b = zeros(n,1); 
b(1) = 1;

if 1; % A = HAH for H = I-2uu'
    u = sparse(n,1);
    u(1:10:n) =  1;
    u(1:20:n) = -1;
    u = u/norm(u); Au = A*u;
    A = A - 2*u*(u'*A) - 2*Au*u' + 4*(u'*Au)*(u*u');
    b = b - 2*(u'*b)*u;
end

if n <= 1000; 
    condA = cond(full(A));
    Lam1 = eig(full(A));
    Lam2 = eig(full(A'*A));
    subplot(121); plot(Lam1,'.'); title('eig(A)');
    subplot(122); plot(Lam2,'.'); title('eig(A''*A)');
    drawnow;
else
    condA = condest(A);
end;
fprintf('cond(A) = %9.2e\n',condA);

tol = 1.e-8; 

if exist('myCGN','file');
fprintf('\nMyCGN:      n = %i\n',n);
t0 = cputime;
[x1,iter] = myCGN(A,b,tol);
cput = cputime - t0;
rres = norm(A*x1 - b)/norm(b);
fprintf(' Iter %4i: rel_res = %6.2e  cpu = %6.2e\n',...
        iter,rres,cput);
end;

% Matlab PCG
Afun = @(x) A'*(A*x);
fprintf('\nMatlab PCG: n = %i\n',n);
t0 = cputime;
[x2,flag,relres,iter] = pcg(Afun,A'*b,tol);
cput = cputime - t0;
rres = norm(A*x2 - b)/norm(b);
fprintf(' Iter %4i: rel_res = %6.2e  cpu = %6.2e\n\n',...
        iter,rres,cput);

% test other Matlab solvers
if n <= 1000
    %warning off;
    [~,f1ag] = qmr(A,b,tol,n);
    if flag > 0; disp('QMR       failed'); end
    [~,flag] = cgs(A,b,tol,n);
    if flag > 0; disp('CGS       failed'); end
    [~,flag] = bicg(A,b,tol,n);
    if flag > 0; disp('BiCG      failed'); end
    [~,flag] = bicgstab(A,b,tol,n);
    if flag > 0; disp('BiCGStab  failed'); end
    [~,flag] = gmres(A,b,10,tol,min(n,500));
    if flag > 0; disp(['GMRES(' int2str(p) ') failed']); end
end
