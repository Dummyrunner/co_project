delta = 0.01;
T = 4;
xt = [0.5;0.5];
xnormbound = 1;
unormbound  =1;
N = T/delta;

Ac = [-1 -1; 0 0.5];
Bc = [0;1];
Cc = eye(2);

sysc = ss(Ac,Bc,Cc,0);
Q = diag([0.5 0.5]);
R = 1;

sysd = c2d(sysc, delta, 'zoh');

[Ad,Bd,~,~] = ssdata(sysd);

[Q,c,Aineq,bineq,Aeq,beq] = mpcZTC2quadprog(Ad,Bd,Q,R,delta,N,xt,xnormbound,unormbound);
