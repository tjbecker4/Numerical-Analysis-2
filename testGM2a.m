% testGM2a: test script for CAAM 454/554 HW 4
% experiments on the relations between the 
% restart number p and cputime

mrange = 50:10:80;
for km = 1:length(mrange); 

m = mrange(km);
fprintf('--- m = %i ---\n',m);
prstr = 'p = %i\t OutIter %4i: res = %6.2e  cpu = %6.2e\n';

a = 2; b = 8; c = 10;
A = fdm4pde1(m,a,b,c);
n = size(A,1);

h = 1/(m+1);
[X,Y] = meshgrid(0:h:1,0:h:1);
F = X .* cos(2*pi*((1-X).^2+(1-Y).^2));
F = 300*F(2:end-1,2:end-1);
f = F(:);

tol = 1e-6;
prange = [2 6 10 20:20:80];
T = zeros(1,length(prange));

for kp = 1:length(prange);
    p = prange(kp);
    t0 = cputime;
    if exist('myGMRES2','file');
        [u,iter] = myGMRES2(A,f,p,tol,n);
    else
        [u,iter] = yzGMRES2(A,f,p,tol,n);
    end;
    T(kp) = cputime - t0;
    res = norm(A*u - f);
    fprintf(prstr,p,iter(1),res,T(kp));
end;

subplot(['22' int2str(km)]);
plot(prange,T);
title(['Mesh size m = ' int2str(m)]);
xlabel('Restarting Value p');
ylabel('CPU time');

end;
