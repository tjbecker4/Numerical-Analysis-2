%          Test script for Jacobi and GS
% First define: method = 'Jacobi' or: method = 'GS' .
% Also define: myname = 'xxxxx xxxxx' , when you are 
% ready to print pictures.
% ------------------------------------------------------

if ~exist('method','var'); method = 'Jacobi'; end

fprintf('\n====== CAAM554/454 Jacobi/GS Tests ======\n')
format1 = 'err = %8.4f  cputime = %5.2f  \n'; 
format2 = ' r1 = %8.4f       r2 = %5.2f  \n'; 
fprintf(['\nMethod is ' method '\n\n'])

N = input(' N = ');
A = gallery('poisson',N);
n = size(A,1);
b = A*ones(n,1);
tol = 5.e-4;
fprintf('\n----- Poisson Matrix, size = %g -----\n',n)  

% ------ run my code ------

   fprintf('Mycode: ')
   t = cputime;
   [myx, myinfo, myres] = feval(['my' method],A,b,10*n,tol);
   mycput = cputime - t;
   myerr = norm(myx-1)/sqrt(n);
   fprintf(format1,myerr,mycput);

% ------ run instructor's code ------
   fprintf('YZcode: ')
   t = cputime;
   [zx, yzinfo, zres] = feval(['z' method '07'],A,b,10*n,tol);
   zcput = cputime - t;
   zerr = norm(zx-1)/sqrt(n);
   fprintf(format1,zerr,zcput);

% ------ calculate ratios ------
   fprintf('ratios: ')
   r1 = myerr/zerr;
   r2 = mycput/zcput;
   fprintf(format2,r1,r2);
   disp(['Myinfo: ' myinfo]);
   disp(['YZinfo: ' yzinfo]);
   
   str1 = ['ERR ratio = ' num2str(r1)]; 
   str2 = ['CPU ratio = ' num2str(r2)]; 
   if exist('myname','var')
      str3 = [myname ', ' date];
   else
      str3 = ['Test Date: ' date];
   end

% ------ plotting ------
m = length(myres);
figure(1)
semilogy(1:m, myres, '.r');
xlabel('Iteration');  
ylabel('Residual')
title([method ' Method Residual Profile'])
text(m/2,8,str1,'fontsize',14); 
text(m/2,4,str2,'fontsize',14); 
text(m/2,2,str3,'fontsize',14); 
set(gca,'fontsize',16)

figure(2); surfc(reshape(myx,[N N]))
title(['My Solution for ' method])
set(gca,'fontsize',14)
fprintf('\n')
