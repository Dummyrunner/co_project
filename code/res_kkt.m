function r = res_kkt(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier)
%RES_KKT compute the current KKT-residual of a quadratic
%optimization problem inside a primal-dual interior point methad.
% For the quadr. convex opt. problem defined by Q,c,Aineq,bineq,Aeq,beq
% (for details see comments on ipquad_pd.m) and a current primal-dual point
% (x,lambda,nu), this function computes the current KKT-residual.
% 
% --------------------------------------
% Input Arguments: for detailed explanation see comments on ipquad_pd.m and
% newtonquad_pd.m.
% --------------------------------------
% Created: 24.06.20, Daniel Bergmann
% --------------------------------------


m = size(Aineq,1);
if ~isempty(Aeq)
    % if there are equality constraints
    r = ...
        [Q*x + c + Aineq'*lambda + Aeq'*nu;...
        -diag(lambda)*(Aineq*x - bineq) - mu_barrier.*ones(m,1);...
        Aeq*x - beq];
else
    % If there are no equality constraints
    r = [Q*x + c + Aineq'*lambda;...
        -diag(lambda)*(Aineq*x - bineq) - mu_barrier.*ones(m,1)];   
end

end