function [x_new, lambda_new, nu_new] = newtonquad_pd(Q, c, Aineq, bineq, Aeq, beq, ls_alpha, ls_beta, x,lambda, nu, mu_barrier)
%NEWTONQUAD_PD Computes search direction for pd-ip-algorithm via Newton's method and step size
%via backtracking line search
% Based on a given primal-dual point x,lambda,nu, this functions returns new points
% x_new, lambda_new, nu_new with smaller kkt-residual.
% First, a search direction is determined by applying newton's method to
% the nonlinear equation system r = 0 with r the residual of the
% kkt-conditions, second a suitable step-size is determined via a
% backtracking linesearch
% -------------------------------------------------------------------------
% - mu_barrier is the current weight on the barrier function
% - for other input arguments, see comments in ipquad_pd.m
% -------------------------------------------------------------------------
% Created: 24.06.20, Daniel Bergmann
%--------------------------------------------------------------------------

% initialize dimensions
n = size(Q,1);
m = size(Aineq,1);
p = size(Aeq,1);

% split stacked r into x, lambda, nu %remove?
% x = r(1:n);
% lambda = r((n+1):(n+m));
% nu = r((n+m+1):(n+m+p));

% Define Matrices for KKT-equality M_kkt*deltar = b_kkt
if ~isempty(Aeq)
    M_kkt =  [Q                        Aineq'                  Aeq';...
        -diag(lambda)*Aineq     -diag(Aineq*x-bineq)      zeros(m,p);
        Aeq                     zeros(p,m)              zeros(p,p) ];
else
    M_kkt =  [Q                        Aineq';
        -diag(lambda)*Aineq     -diag(Aineq*x-bineq)];
end

b_kkt = -res_kkt(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier);

% Solve KKT-equality to get search direction
deltaxln = M_kkt\b_kkt;


deltax = deltaxln(1:n);
deltalambda = deltaxln((n+1):(n+m));
deltanu = deltaxln((n+m+1):(n+m+p));

% Perform Backtracking linesearch for determining suitable step size

% compute maximal step size smax
lambdaquot = -lambda./deltalambda;
% eliminate zeros and entries with deltalambda greater or equal zero
lambdaquot = lambdaquot - (deltalambda >= 0).*lambdaquot;
lambdaquot(lambdaquot == 0) = [];
smax = min([1 ; lambdaquot]);

% Backtracking-Linesearch
% Typical Parameter choices:
% alpha in [0.01,0.1]; beta in [0.3,0.8]
s = 0.99*smax; %S4 0.99 smax??
found = 0;
while found == 0
    s = s*ls_beta;
    x_new = x + s.*deltax;
    lambda_new = lambda + s.*deltalambda;
    nu_new = nu + s.*deltanu;
    r_new = res_kkt(x_new,lambda_new,nu_new,Q,c,Aineq,bineq,Aeq,beq,mu_barrier);
    r_old = res_kkt(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier);
    if all(Aineq*x_new - bineq < 0) && (norm(r_new) <= (1-ls_alpha*s)*norm(r_old))
        found = 1;
    end
end
disp(['Backtracking search yields step size s = ',num2str(s)])
end

