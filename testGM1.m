% test script for GMRES
clear; close all;

%m = 40; 
m = input('m = ');
a = 2; b = 8; c = 10;

A = fdm4pde1(m,a,b,c);
n = size(A,1);
fprintf('n = %g\n',n);

h = 1/(m+1);
[X,Y] = meshgrid(0:h:1,0:h:1);
F = X .* cos(2*pi*((1-X).^2+(1-Y).^2));
F = 300*F(2:end-1,2:end-1);
f = F(:);

tol = 1e-6;
       
if exist('myGMRES1','file')
    t0 = cputime;
    [u1,iter] = myGMRES1(A,f,tol);
    cput = cputime - t0;
    res = norm(A*u1 - f);
    fprintf('\nmyGMRES1 code:\n');
    fprintf('Iter %i: res = %6.2e  cpu = %6.2e\n',...
        iter,res,cput);
end

if exist('yzGMRES1','file');
    t0 = cputime;
    [u2,iter] = yzGMRES1(A,f,tol);
    cput = cputime - t0;
    res = norm(A*u2 - f);
    fprintf('\nyzGMRES1 code:\n');
    fprintf('Iter %i: res = %6.2e  cpu = %6.2e\n',...
        iter,res,cput);
end

% Matlab gmres
if 1;
    t0 = cputime;
    [u3,~,~,iter] = gmres(A,f,n,tol,n);
    cput = cputime - t0;
    res = norm(A*u3 - f);
    fprintf('\nMatlab gmres code:\n');
    fprintf('Iter %i: res = %6.2e  cpu = %6.2e\n\n',...
        iter(2),res,cput);
end

str1 = []; str2 = [];
fprintf('\nRelative Difference from matlab gmres:\n');
if exist('u1','var');
    rdf = norm(u1-u3)/norm(u3);
    fprintf('My solution: %g\n',rdf);
    U = reshape(u1,m,m); 
    str1 = 'My solution:  ';
end

if exist('u2','var');
    rdf = norm(u2-u3)/norm(u3);
    fprintf('YZ solution: %g\n',rdf);
    U = reshape(u2,m,m); 
    str2 = 'YZ solution:  ';
end

str = ['m = ' int2str(m) ',   (a,b,c) ='];
str = [str ' (' num2str(a)];
str = [str ', ' num2str(b)];
str = [str ', ' num2str(c) ')'];

figure(1); surfc(U);
title ([str1 str2 str],'fontsize',18);