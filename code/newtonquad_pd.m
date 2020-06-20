function [r_new] = newtonquad_pd(Q, c, Aineq, bineq, beq, Aeq, arm_beta, arm_gamma, r,mu_barrier, tol)
%UNTITLED Summary of this function goes here
%   f(x) = 0.5*x'Qx + c'x

% initialize dimensions
n = size(Q,1);
m = size(Aineq,1);
p = size(Aeq,1);

% split stacked r to x, lambda, nu
x = r(1:n);
lambda = r((n+1):(n+m));
nu = r((n+m+1):(n+m+p));

% Define Matrices for KKT-equality M_kkt*deltar = b_kkt
 M_kkt =  [Q                        Aineq'                  Aeq';...
            -diag(lambda)*Aineq     -diag(Aineq*x-beq)      zeros(m,p);
            Aeq                     zeros(p,m)              zeros(m,p) ];

b_kkt = [Q*x + c + Aineq'*lambda + Aeq'*nu;...
          -diag(lambda)*(Aineq*x - bineq) - mu_barrier.*eye(m,1);...
          Aeq*x - beq];

r_new = [1;1];

end

