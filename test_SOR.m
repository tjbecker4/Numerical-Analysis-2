function test_SOR(matrix)
% ------------------------------------------------------
%              test function for SOR
% ------------------------------------------------------

if ~nargin; matrix = 1; end;

fprintf('\n  ====== CAAM 454/554 SOR Testing ======\n')
format1 = 'Omega = %5.2f  Iter = %5i  rerr = %8.4e\n'; 
format2a = ': Total Cputime = %g\n\n'; 

N = input(' N = ');
Range = input(' omega range = ');
do_plot = numel(Range) > 2;

switch matrix
    case 1; A = gallery('poisson',N);
    case 2; A = gallery('wathen',N,N);
    otherwise; error('matrix must be 1 or 2');
end
n = size(A,1);
b = A*ones(n,1);
tol = 1e-4;

Str = [' Poisson '; ' Wathen  '];
fprintf(['\n\t---' Str(matrix,:) 'Matrix, n = %i ---\n'],n)

fprintf('\nRun zSOR:\n')  
solver = 'zSOR';
format2 = [solver format2a];
if do_plot; figure(1); end
run_solver;

if exist('mySOR','file')
    fprintf('\nRun mySOR:\n')  
    solver = 'mySOR';
    format2 = [solver format2a]; 
    if do_plot; figure(2); end
    run_solver;
end


%%%%%%%%%%%%%%%%%%%
function run_solver

IT = zeros(length(Range),1); 
k = 1; RE = IT; t0 = tic;
for omega = Range
   [x, iter] = feval(solver,A,b,omega,50*n,tol/n);
   rerr = norm(x-1)/sqrt(n);
   IT(k) = iter; RE(k) = rerr;
   fprintf(format1,omega,iter,rerr);
   k = k + 1;
end
fprintf(format2,toc(t0));

% ------ plotting ------
if do_plot
    plot(Range, IT, Range, IT, '*');
    title([solver ': Omega vs. Iterarion'])
    ylabel('Iteration');  xlabel('Omega')
    [mint, at] = min(IT); bestO = Range(at);
    str = ['Min Iter = ' num2str(mint) ];
    str = [str ' at omega = ' num2str(bestO)];
    text(Range(1), (max(IT)+min(IT))/2, str,'fontsize',14);
    set(gca,'fontsize',14)
end

end %run_solver

end %main
