function [Q,c,Aineq,bineq,Aeq,beq] = mpcZTC2quadprog(Ad,Bd,Q,R,delta,N,xinit,xnormbound,unormbound)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   reprasing ztc mpc at time t, horizon N samples
n = size(Ad,1);
m = size(Bd,2);

% Build Matrix H
Qs = repmat({Q},N+1,1);
Rs = repmat({R},N,1);
H = zeros(N*(n+m)+n);
Qblk = blkdiag(Qs{:});
Rblk = blkdiag(Rs{:});
H = delta.*blkdiag(Qblk,zeros(n),Rblk);

% express box constraints for pred. state and input in inequality constr.
Aineq = [xnormbound.*eye(n) zeros(n,m);...
         -xnormbound.*eye(n) zeros(n,m);...
         zeros(m,n) unormbound.*eye(m);...
         zeros(m,n) -unormbound.*eye(m)];
bineq = ones(2*(n+m));

% express system dynamics and init. cond. in equality constr.
Aeqx = kron(eye(N+1),eye(n));  Aeqx(1:n,1:n) = -1.*Aeqx(1:n,1:n);
Aeqx = Aeqx + kron(diag(ones(N,1),-1),Ad);
lastrow = zeros(1,N+1); lastrow(end) = 1;
Aeqx = [Aeqx; kron(lastrow,eye(n))];

Aequ = kron(eye(N),Bd);
Aequ = [zeros(N,N*m); Aequ];
beq = [xinit; zeros(N*m,1)];
end

