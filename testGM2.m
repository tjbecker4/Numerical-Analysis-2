% test script for GMRES(p)
clear;

m = input('mesh size m = ');
a = 2; b = 8; c = 10;
A = fdm4pde1(m,a,b,c);
n = size(A,1);

h = 1/(m+1);
[X,Y] = meshgrid(0:h:1,0:h:1);
F = X .* cos(2*pi*((1-X).^2+(1-Y).^2));
F = 300*F(2:end-1,2:end-1);
f = F(:); fnrm = norm(f);

tol = 1e-6;
prange = [10 20];

for p = prange; % p-loop

fprintf('\n----------------------------------------\n');
fprintf(' *** n = %g, tol = %g, p = %g ***\n',n,tol,p);
fprintf('----------------------------------------\n');

if exist('myGMRES2','file')
    t0 = cputime;
    [u1,it1] = myGMRES2(A,f,p,tol,n);
    cput = cputime - t0;
    res = norm(A*u1 - f)/fnrm;
    iter1 = p*(it1(1)-1)+it1(2);
    fprintf('myGMRES2: ');
    fprintf('Iter %i rres = %6.2e cpu = %6.2e\n',...
        iter1,res,cput);
end

if exist('yzGMRES2','file');
    t0 = cputime;
    [u2,it2] = yzGMRES2(A,f,p,tol,n);
    cput = cputime - t0;
    res = norm(A*u2 - f)/fnrm;
    iter2 = p*(it2(1)-1)+it2(2);
    fprintf('yzGMRES2: ');
    fprintf('Iter %i rres = %6.2e cpu = %6.2e\n',...
        iter2,res,cput);
end

% Matlab gmres: turned off 
if 0 < 1; 
    t0 = cputime;
    [u3,tmp1,tmp2,it3] = gmres(A,f,p,tol,n);
    cput = cputime - t0;
    res = norm(A*u3 - f)/fnrm;
    iter3 = p*(it3(1)-1)+it3(2);
    fprintf('  Matlab: ');
    fprintf('Iter %i rres = %6.2e cpu = %6.2e\n',...
        iter3,res,cput);
end
    fprintf('\n');
end; % p-loop
