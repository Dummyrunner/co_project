delta = 0.1;
T = 0.2;
xnormbound = 1;
unormbound = 1;
N = T/delta;

Ac = [-1 -1; 0 0.5];
Bc = [0;1];
Cc = eye(2);

m = size(Bc, 2);
n = size(Ac, 2);

sysc = ss(Ac,Bc,Cc,0);
Q = diag([0.5 0.5]);
R = 1;

sysd = c2d(sysc, delta, 'zoh');

[Ad,Bd,~,~] = ssdata(sysd);


% linesearch parameters
ls_alpha = 0.05;
ls_beta = 0.5;

% initial predicted signals
u0 = repmat(zeros(m,1),N,1);
x0 = repmat([0.6; 0.8],N+1,1);
xtilde_init = [x0; u0];

% state at initial time t
xt = x0(1:n);

[H,c,Aineq,bineq,Aeq,beq] = ...
              mpcZTC2quadprog(Ad,Bd,Q,R,delta,N,xt,xnormbound,unormbound);
          
lambda0 = arrayfun(@(x) -1/x, Aineq*xtilde_init - bineq);
p = size(Aeq,1);
nu0 = zeros(p,1);
% decreasement of mu and tolerances
gamma = .1;
eps_feas = 1e-3;
eps_opt = 1e-3;

tic
res = quadprog(H,c,Aineq, bineq, Aeq, beq)
% ipquad_pd(H,c,Aineq,bineq,Aeq,beq,xtilde_init,lambda0,nu0,gamma,eps_feas,eps_opt,ls_alpha,ls_beta)
toc