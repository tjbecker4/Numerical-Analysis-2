function [x,y,z,iter] = pdipm(A,b,c,tol)

%
% Tim Becker
% CAAM 554 Homework 10 Problem 1
% 4/11/2013
%
% pdipm.m uses the interior point primal dual method for solving the linear
% program
%
% [x,y,z,iter] = pdipm(A,b,c,tol)
%
% INPUT:
%
% A      constraint coefficient matrix
% b      constraint right-hand side vector
% c      objective coefficient vector
% tol    tolerance
%
% OUTPUT:
%
% x      final primal solution
% y      final dual solution
% z      final dual slacks
% iter   iteration counter

% Initialization
[m,n] = size(A);
I = speye(m);

fprintf(' --- [m,n] = [%3i %3i] ---\n',m,n);

p = symamd(abs(A)*abs(A'));

x = ones(n,1)*n;
z = ones(n,1)*n;
y = zeros(m,1);
dy = zeros(m,1);

nrmc = 1+norm(c);
nrmb = 1+norm(b);

primal = norm(A*x-b)/nrmb;
dual = norm(A'*y+z-c)/nrmc;
gap = abs(b'*y-c'*x)/(1+abs(b'*y));

rd = c-A'*y-z;
rp = b-A*x;
    

for iter = 1:100
    
    % Print out each iteration
    fprintf('iter %3i: [primal dual gap] = [%.2e %.2e %.2e]\n',iter,primal,dual,gap);
    
    mu = (x'*z)/n;
    theta = min(0.2,10*mu);
    
    rc = theta*mu-x.*z;
    
    zinv = 1./z;
    t = x.*zinv;
    q = rc.*zinv;
    d = min(t,10^15);
    B = A*sparse(1:n,1:n,d)*A';
    %B = (bsxfun(@times,A,d'))*A';
    rhs = A*(d.*rd-q)+rp;
    
    % Find dx,dy,dz
    
    %R = cholinc(B(p,p) + sparse(1:m,1:m,1.e-10),'inf');
    R = chol(B(p,p) + sparse(1:m,1:m,1.e-10));
    dy(p) = R\(R'\rhs(p));
    dz = rd - A'*dy;
    dx = q-d.*dz;
    
    % Generate step length
    alphap = -1/min([dx./x;-1]);
    alphad = -1/min([dz.*zinv;-1]);
    
    gamma = max(0.995,1-10*mu);
    
    % Update
    w = min(1,gamma*alphad);
    x = x + min(1,gamma*alphap)*dx;
    y = y + w*dy;
    z = z + w*dz;
    
    rd = c-A'*y-z;
    rp = b-A*x;
    
    by = b'*y;
    primal = norm(rp)/nrmb;
    dual = norm(rd)/nrmc;
    gap = abs(by-c'*x)/(1+abs(by));
    
    crit = max([primal,dual,gap]);
    
    % Check if converged
    if crit < tol
        break
    end
    
end



