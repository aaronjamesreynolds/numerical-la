function [x, iter, diff_norm] = Jacobi(A,b,x0,x_exact)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Jacobi Iterative Solver
%
% Author: Aaron J Reynolds, for MTH 551: Numerical Linear Algebra
% -------------------------------------------------------------------------
% This function uses the Jacobi method to calculate the solution to 
% a linear system of the form Ax=b. It returns the numerical solution, the 
% number of iterations required to reach it, and the inf-norm of the 
% difference between the exact and numerical solutions
% 
% INPUTS
% A,b: of Ax = b, x0: initial guess for x of Ax=b, x_exact: exact solution
% 
% OUTPUTS
% x: numerical solution, iter: number of iterations required,
% diff_norm: inf-norm of difference between x and x_exact at each iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize parameters and vectors
max_iterations = 5000; % maximum number of iterations if no sol'n found
tolerance = 10^(-12); % convergence criteria
N = length(x0); % number of columns in matrix A
x_old = x0; x_new = zeros(N,1);  % intialize solutions vectors
%% Jacobi iterative solver
tic
for iter = 1:max_iterations
    for i = 1:N
        sum = 0;
       for j = 1:N
          if i~=j
              sum = sum + A(i,j) * x_old(j);
          end
       end
       x_new(i) = (b(i) - sum)/A(i,i); 
    end
    diff_norm(iter) = norm(x_exact-x_new,inf);
    if diff_norm(iter) < tolerance, break; end
    x_old = x_new;
end
x = x_new;
%% Print statements 
if iter == max_iterations 
    fprintf('MAXIMUM NUMBER OF ITERATIONS REACHED ');
end
fprintf('--Jacobi solve--\n'); toc;
fprintf('Iterations: %g \n',iter);
fprintf('Absolute error: %g \n', diff_norm(iter));
fprintf('---\n\n')
end
