% fpor NewtonOpt (primal only)
% q = [ 1; 3];
% r = 0.5;
% f = @(x) x'*P*x + q'*x + r;
% gradf = @(x) P*x + q;
% Hessf = @(x) P;

Aineq = [0 1]; bineq = 10;
Aeq = [1 -1]; beq = 0;

Q = diag([1 2]);
c = [1;2];
x0 = [5;5];

arm_beta = 0.5;
arm_gamma = 0.01;
r = [5;5;1;1];
mu_barrier = 0.5;
l = newtonquad_pd(Q,c,Aineq, bineq, beq, Aeq, arm_beta, arm_gamma,r,0.9, 1e-5);