% Initialization of Example Optimization Problem

% Matrices defining constraints
Aineq = [0 1]; bineq = 10;
Aeq = [1 -1]; beq = 0;

% Matrices defining objective
Q = diag([1 2]);
c = [1;2];

%  m,p = 1
%  n = 2

% linesearch parameters
ls_alpha = 0.05;
ls_beta = 0.5;

% decreasement of mu and tolerances
gamma = .1;
eps_feas = 1e-4;
eps_opt = 1e-4;

% initial primal dual point
x0 = [1;1];
lambda0 = 1;
nu0 = 1;

[x,fval,lambda,nu,eta] = ipquad_pd(Q,c,Aineq,bineq,Aeq,beq,x0,lambda0,nu0,gamma,eps_feas,eps_opt,ls_alpha,ls_beta)

%-------------Apply quadprog() to double check-------------------
% [quadprog_x,quadprog_fval] = quadprog(Q,c)

