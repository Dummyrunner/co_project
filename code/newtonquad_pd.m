function [x_new, lambda_new, nu_new] = newtonquad_pd(Q, c, Aineq, bineq, beq, Aeq, ls_alpha, ls_beta, x,lambda, mu_barrier)
%UNTITLED Summary of this function goes here
%   f(x) = 0.5*x'Qx + c'x

% initialize dimensions
n = size(Q,1);
m = size(Aineq,1);
p = size(Aeq,1);

% split stacked r into x, lambda, nu %remove?
% x = r(1:n);
% lambda = r((n+1):(n+m));
% nu = r((n+m+1):(n+m+p));

% Define Matrices for KKT-equality M_kkt*deltar = b_kkt
 M_kkt =  [Q                        Aineq'                  Aeq';...
            -diag(lambda)*Aineq     -diag(Aineq*x-beq)      zeros(m,p);
            Aeq                     zeros(p,m)              zeros(m,p) ];

b_kkt = -rfromxln(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier);

% Solve KKT-equality to get search direction
deltaxln = M_kkt\b_kkt;


deltax = deltaxln(1:n);
deltalambda = deltaxln((n+1):(n+m));
deltanu = deltaxln((n+m+1):(n+m+p));

% Perform Backtracking linesearch for determining suitable step size

% compute maximal step size smax
lambdaquot = -lambda./deltalambda;
lambdaquot = lambdaquot + (deltalambda >= 0).*lambdaquot;
smax = min([1 ; lambdaquot]);

% Backtracking-Linesearch
% Typical Parameter choices:
% alpha in [0.01,0.1]; beta in [0.3,0.8]
s = smax; %S4 0.99 smax??
found = 0;

while found == 0
    s = s*ls_beta;
    x_new = x + s.*deltax;
    lambda_new = lambda + s.*deltalambda;
    nu_new = nu + s.*deltanu;
    
    %TODO compute r
    if Aineq*x_new < 0 && norm(rfromxln(x_new,lambda_new,nu_new,Q,c,Aineq,Aeq,beq,mu_barrier)) <= (1-ls_alpha*s)*norm(rfromxln(x,lambda,nu,Q,c,Aineq,Aeq,beq,mu_barrier))
        found = 1;
    end
end
end


% -----------------------------------------------------------------------
function r = rfromxln(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier)
          r = ...
              [Q*x + c + Aineq'*lambda + Aeq'*nu;...
              -diag(lambda)*(Aineq*x - bineq) - mu_barrier.*eye(m,1);...
              Aeq*x - beq];
end