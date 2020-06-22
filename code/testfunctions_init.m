% Initialization of Example Optimization Problem


% for NewtonOpt (primal only)
% q = [ 1; 3];
% r = 0.5;
% f = @(x) x'*P*x + q'*x + r;
% gradf = @(x) P*x + q;
% Hessf = @(x) P;
% arm_beta = 0.5;
% arm_gamma = 0.01;

Aineq = [0 1]; bineq = 10;
Aeq = [1 -1]; beq = 0;

Q = diag([1 2]);
c = [1;2];
% x0 = [5;5];

%  m,p = 1
%  n = 2


ls_alpha = 0.05;
ls_beta = 0.5;
r = [5;5;1;1];
% mu_barrier = 0.5;
% l = newtonquad_pd(Q, c, Aineq, bineq, Aeq, beq, ls_alpha, ls_beta,[1; 1],1,1,0.9);

gamma = .1;
eps_feas = 1e-4;
eps_opt = 1e-4;

[x,fval,lambda,nu,eta] = ipquad_pd(Q,c,Aineq,bineq,Aeq,beq,[1; 1],1,1,gamma,eps_feas,eps_opt,ls_alpha,ls_beta);