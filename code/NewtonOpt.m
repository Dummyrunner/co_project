function [outputArg1,outputArg2] = NewtonOpt(f, gradf, Hessf, x0, arm_beta, arm_gamma, Aeq, tol)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%  Example for Armijo Parameters: arm_beta = 0.5, arm_gamma = 0.01
n = length(x0);
p = size(Aeq,1);
lambda = 1;
x = x0;

while lambda^2/2 > tol
    Hessfinv = inv(Hessf(x));
    % For newton step, solve linear equation A_kkt * [x_newton; nu ] = b_kkt
    % resulting from KKT conditions
    % for equality constrained Problem
    A_kkt = [Hessf(x) Aeq';...
             Aeq        zeros(p)];
    b_kkt = [-gradf(x); zeros(p,1)];
    x_newton = A_kkt\b_kkt;
    x_newton = x_newton(1:n);
    
   % Compute suitable step-size via backtracking search and Amijo-rule
   s = 1/arm_beta;
   while f(x + s.*x_newton) > f(x) + arm_gamma*s*(gradf(x)'*x_newton)
      s = arm_beta*s;
   end
   
   % Update of examined point x
   x = x + s.*x_newton;
end
end
