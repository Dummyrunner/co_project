function [x,fval,lambda,nu,eta] = ipquad_pd(Q,c,Aineq,bineq,Aeq,beq,x0,lambda0,nu0,gamma,eps_feas,eps_opt,ls_alpha,ls_beta)
%IPQUAD_PD Quadratic optimization via primal-dual-interior-point-method.
% Convex Quadr. function f(x) = (1/2)x'Qx + c'x with linear equality and
% inequality constraints.
% -------------------------------------------------------------------------
% Input Arguments:
% - Q,c define objective function
%   dim(Q) = nxn, dim c = nx1
% - Aineq, bineq define inequality constraints Aineq*x <= bineq
%   dim(Aineq) = mxn, dim(bineq) = mx1
% - Aeq, beq define inequality constraints Aeq*x == beq
%   dim(Aineq) = pxn, dim(beq) = px1
% - x0 init. value for the primal problem. lambda0, nu0. labda0 >= 0.
%   Aineq*x0 <= bineq.
%   initial values for the dual problem
% - gamma is a reduction factor for reducing the barrier weight mu_barrier
%   in each iteration. gamma in (0,1)
% - eps_feas > 0 specifies the tolerance for the 2norms of the primal and the
%   dual residual
% - eps_opt > 0 specifies a tolerance on the surrogate duality gap
% - ls_alpha, ls_beta are parameters for the backtracking linesearch,
%   performed in each iteration. Typical choices: ls_alpha in [0.01,0.1].
%   ls_beta in [0.3,0.8]
%--------------------------------------------------------------------------
% Created: 24.06.20, Daniel Bergmann
%--------------------------------------------------------------------------


x = x0;
lambda = lambda0;
nu = nu0;

% initialize dimensions
n = size(Q,1);
m = size(Aineq,1);
p = size(Aeq,1);

% TODO init dual variables if not passed to the algo. HOW?
found = 0;
count = 0;

while found == 0
    % compute surrogate duality gap
    eta = -(Aineq*x - bineq)'*lambda;
    mu_barrier = gamma*eta/m;
    % Compute KKT residual vector
    r_mu = res_kkt(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier);
    r_dual = r_mu(1:m);
%     r_cent = r_mu((m+1):(m+p));
    r_pri = r_mu((m+p+1):(m+p+n));
    
    if norm(r_pri) <= eps_feas && norm(r_dual) <= eps_feas && eta <= eps_opt
        found=1;
        break;
    else
        % update barrier weighting parmeter mu_barrier
        mu_barrier = gamma*eta/m;
        
        % compute search vector
        [x, lambda, nu] = newtonquad_pd(Q, c, Aineq, bineq, Aeq, beq, ls_alpha, ls_beta, x,lambda, nu, mu_barrier);
    end
    count = count+1;
    
    disp(['Iteration No. ',num2str(count)])
end


% compute value of objective at x
fval =  0.5.*x'*Q*x + c'*x;
end

