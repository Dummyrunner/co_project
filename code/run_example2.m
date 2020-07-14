% Initialization of Example Optimization Problem

% Matrices defining constraints
Aineq = [-2 -3; -1 0; 0 -1]; bineq = [-4;0;0];
Aeq = []; beq = [];

% Matrices defining objective
Q = 2.*[3 1; 1 1];
c = [1; 6];

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
lambda0 = [0.2;1; 1];
nu0 = 0;

[x,fval,lambda,nu,eta] = ipquad_pd(Q,c,Aineq,bineq,Aeq,beq,x0,lambda0,nu0,gamma,eps_feas,eps_opt,ls_alpha,ls_beta)

%-------------Apply quadprog() to double check-------------------
% [quadprog_x,quadprog_fval] = quadprog(Q,c)