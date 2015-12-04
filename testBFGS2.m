% CAAM 454/554 - test script for BFGS

func = @yzlork;
if exist('mylork','file')
    func = @mylork;
end
n = 200; k = 5;

%n = input('n = ');
%k = input('k = ');
tol = 1.e-9;
maxit = 300;

T = randn(n,2*k);
A = T*diag(1./(1:2*k))*T';
nrmA2 = norm(A,'fro')^2;

[V,D] = eigs(A,k);
Xs = V*sqrt(D); 
fs = func(Xs(:),n,k,A,nrmA2);

fprintf('n = %i, k = %i, fs = %.8e\n',n,k,fs)

x01 = eye(n,k);
x01 = x01(:);
x02 = randn(n*k,1);

Solvers(1,:) = 'yzBFGS';
Solvers(2,:) = 'myBFGS';

for j = 1:2
    solver = Solvers(j,:);
    if exist(solver,'file')
        t0 = cputime;
        disp('Initial Guess 1 ...')
        [x1,iter1,hist] = ...
            eval([solver '(func,x01,tol,maxit,n,k,A,nrmA2)']);
        f1 = func(x1,n,k,A,nrmA2);
        fh1 = hist(1,:); gh1 = hist(2,:);
        fprintf('f1 - fs = %.8e\n\n',f1-fs)
        disp('Initial Guess 2 ...')
        [x2,iter2,hist] = ...
            eval([solver '(func,x02,tol,maxit,n,k,A,nrmA2)']);
        f2 = func(x2,n,k,A,nrmA2);
        fprintf('f2 - fs = %.8e\n\n',f2-fs)
        fh2 = hist(1,:); gh2 = hist(2,:);
        fprintf([solver ': total CPU = %8.4e\n\n'],cputime-t0);
        figure(j);
        subplot(221); semilogy(1:length(fh1),fh1,'b-*');
        xlabel('Iteration vs. f'); ylabel('func')
        title([solver ', x0 = eye']);
        grid on; set(gca,'fontsize',14)
        subplot(223); semilogy(1:length(gh1),gh1,'b-*');
        xlabel('Iteration vs. g'); ylabel('||g||_2')
        title([solver ', x0 = eye']);
        grid on; set(gca,'fontsize',14)
        subplot(222); semilogy(1:length(fh2),fh2,'r-*');
        xlabel('Iteration vs. f'); ylabel('func')
        title([solver ', x0 = randn']);
        grid on; set(gca,'fontsize',14)
        subplot(224); semilogy(1:length(gh2),gh2,'r-*');
        xlabel('Iteration vs. g'); ylabel('||g||_2')
        title([solver ', x0 = randn']);
        grid on; set(gca,'fontsize',14)
    end
end