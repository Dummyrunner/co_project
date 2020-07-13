delta = 0.01;
T = 4;
xt = [0.5;0.5];
xnormbound = 1;
unormbound  =1;
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


       
          
% uinit = zeros(size(Bc,2));
% xtilde_init = [xt; uinit];


% linesearch parameters
ls_alpha = 0.05;
ls_beta = 0.5;

u0              = repmat(zeros(m,1),N,1);
x0              = repmat([0.6; 0.8],N+1,1);
xtilde_init = [u0;x0];

xt = x0(1:n);

[H,c,Aineq,bineq,Aeq,beq] = ...
              mpcZTC2quadprog(Ad,Bd,Q,R,delta,N,xt,xnormbound,unormbound);
          
lambda0 = arrayfun(@(x) -1/x, Aineq*[xtilde_init] - bineq);
nu0 = 0;
% decreasement of mu and tolerances
gamma = .1;
eps_feas = 1e-4;
eps_opt = 1e-4;


ipquad_pd(H,c,Aineq,bineq,Aeq,beq,xtilde_init,lambda0,nu0,gamma,eps_feas,eps_opt,ls_alpha,ls_beta)