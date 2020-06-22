function r = rfromxln(x,lambda,nu,Q,c,Aineq,bineq,Aeq,beq,mu_barrier,m)
%RFROMXLN bla
% blabla
r = ...
    [Q*x + c + Aineq'*lambda + Aeq'*nu;...
    -diag(lambda)*(Aineq*x - bineq) - mu_barrier.*ones(m,1);...
    Aeq*x - beq];
end