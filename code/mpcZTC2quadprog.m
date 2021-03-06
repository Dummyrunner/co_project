function [H,c,Aineq,bineq,Aeq,beq] = mpcZTC2quadprog(Ad,Bd,Q,R,delta,N,xt,xnormbound,unormbound)
%MPCZTC2QUADPROG Translates a MPC problem with zero terminal constraints
%into a quadratic program
%   Model predictive control optimization problem of a time-discrete
%   system determined by Ad, Bd with sample time delta and
%   prediction horizon N, with current state xt and bounds on the maximum
%   norm of state x and input u. Q,R define the objective that has to be
%   minimized

n = size(Ad,1);
m = size(Bd,2);

% Build Matrix H
Qs = repmat({Q},N,1);
Rs = repmat({R},N,1);
Qblk = blkdiag(Qs{:});
Rblk = blkdiag(Rs{:});
H = delta.*blkdiag(Qblk,zeros(n),Rblk);

% express box constraints for pred. state and input in inequality constr.
Aineq = [xnormbound.*eye(n*(N+1)) zeros(n*(N+1),m*N);...
         -xnormbound.*eye(n*(N+1)) zeros(n*(N+1),m*N);...
         zeros(m*N,n*(N+1)) unormbound.*eye(m*N);...
         zeros(m*N,n*(N+1)) -unormbound.*eye(m*N)];
bineq = ones(2*(n*(N+1) + m*N),1);

% express system dynamics and init. cond. in equality constr.
Aeqx = kron(eye(N+1),eye(n));
Aeqx(n+1:end,n+1:end) = -1.*Aeqx(n+1:end,n+1:end);
Aeqx = Aeqx + kron(diag(ones(N,1),-1),Ad);
lastrow = zeros(1,N+1); lastrow(end) = 1;
Aeqx = [Aeqx; kron(lastrow,eye(n))];


Aequ = kron(eye(N),Bd);
Aequ = [zeros(n,N*m); Aequ; zeros(n,N*m)];

beq = [xt; zeros((N+1)*n,1)];

Aeq = [Aeqx Aequ];

c = zeros(N*(m+n)+n,1);
end

