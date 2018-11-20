function [x, iter, diff_norm] = SOR(A,b,x0,x_exact,omega)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SOR Solver
%
% Author: Aaron J Reynolds, for MTH 551: Numerical Linear Algebra
% -------------------------------------------------------------------------
% This function uses the SOR method to calculate the solution to 
% a linear system of the form Ax=b. It returns the numerical solution, the 
% number of iterations required to reach it, and the inf-norm of the 
% difference between the exact and numerical solutions
% 
% INPUTS
% A,b: of Ax = b, x0: initial guess for x of Ax=b, x_exact: exact solution
% omega: relaxation factor of SOR method 
%
% OUTPUTS
% x: numerical solution, iter: number of iterations required,
% diff_norm: inf-norm of difference between x and x_exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize parameters and vectors
max_iterations = 5000; % maximum number of iterations if no sol'n found
tolerance = 10^(-12); % convergence criteria
N = length(x0); % number of columns in matrix A
x_old = x0; x_new = zeros(N,1);  % intialize solutions vectors
%% SOR iterative solver
tic
for iter = 1:max_iterations
    for i = 1:N
        sum = 0;
        for j = 1:i-1
            sum = sum + A(i,j)*x_new(j); % use new information
        end
        for j = i:N
           sum = sum + A(i,j)*x_old(j); % use old information
        end
        sum = (b(i) - sum)/A(i,i);
        x_new(i) = x_old(i) + omega*(sum); 
    end
    diff_norm(iter) = norm(x_exact-x_new, inf);
    if diff_norm(iter) < tolerance, break; end
    x_old = x_new;
end
x = x_new;
%% Print statements
if omega == 1
    fprintf('--Gauss Seidel solve--\n')
else
    fprintf('--SOR solve, omega = %g--\n', omega)
end
if iter == max_iterations 
    fprintf('MAXIMUM NUMBER OF ITERATIONS REACHED ');
end
toc; fprintf('Iterations: %g \n',iter);
fprintf('Absolute error: %g \n', diff_norm(iter));
fprintf('---\n\n')
end
