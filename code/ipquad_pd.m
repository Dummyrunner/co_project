function [x,fval,lambda,nu,eta] = ipquad_pd(Q,c,Aineq,bineq,Aeq,beq,x0,lambda0,nu0,gamma,eps_feas,eps_opt,ls_alpha,ls_beta)
%UNTITLED Summary of this function goes here
% Aineq*x0 - bineq < 0,  lambda > 0, eps_feas,eps_opt > 0, gamma in (0,1)
count = 0;
x = x0;
lambda = lambda0;
nu = nu0;

n = size(Q,1);
m = size(Aineq,1);
p = size(Aeq,1);

% TODO init dual variables if not passed to the algo. HOW?
found = 0;

while found == false
    % compute surrogate duality gap
    eta = -(Aineq*x - bineq)'*lambda;
    mu_barrier = gamma*eta/m;
    % Compute KKT residual vector
    r_mu = rfromxln(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier);
    r_dual = r_mu(1:m);
    r_cent = r_mu((m+1):(m+p));
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
end


% compute value of objective at x
f = @(x) x'*Q*x + c'*x;
fval = f(x);
end

