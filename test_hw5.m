% test_hw5

n = input('n = ');
% No. of var = n^2
fprintf('\tNumber of variables: %i\n',n^2)

% exact solution is Z
[~,~,Z] = peaks(n);
% generate A p.s.d and B
A = gallery('minij',n);
B = A*Z + Z*A;

ferr = @(x,y) norm(x-y)/norm(y);
tol = 1e-8; maxit = n^2;

% call yzSolve4X
if exist('yzSolve4X','file')
    tic; X = yzSolve4X(A,B,tol,maxit); t = toc;
    fprintf('yzSolve4X: RelErr = %6.2e time = %f\n',...
        ferr(X(:),Z(:)),t)
end

% call mySolve4X
if exist('mySolve4X','file')
    tic; X = mySolve4X(A,B,tol,maxit); t = toc;
    fprintf('mySolve4X: RelErr = %6.2e time = %f\n',...
        ferr(X(:),Z(:)),t)
end

if n <= 150
    surf(X-Z); 
    h = title('Error Plot');
    set(h,'fontsize',16)
    set(gca,'fontsize',16)
end
