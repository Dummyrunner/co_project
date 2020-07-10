function [Q,c,Aineq,Aeq] = mpcZTC2quadprog(Ad,Bd,Q,R,delta,N,xinit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   reprasing ztc mpc at time t, horizon N samples
n = size(Ad,1);
m = size(Bd,2);

Qs = repmat({Q},N+1,1);
Rs = repmat({R},N,1);
H = zeros(N*(n+m)+n);
Qblk = blkdiag(Qs{:});
Rblk = blkdiag(Rs{:});
H = blkdiag(Qblk,zeros(n),Rblk)
end

