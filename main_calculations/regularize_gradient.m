


function [f_reg, f_prime_reg, f_prime_reg_interpolated, norm_cost, x_grad_mesh] = regularize_gradient(f_in, x_mesh, alpha)


% Determine data length
n = length(f_in);
% Make column (vector)
f_in = reshape(f_in, [n, 1]);
% % Make string
% x_steps = reshape(x_steps, [1, n]);
% Shift the left border to 0
f_shift = f_in(1);
f_in = f_in - f_shift;
% Normalize std to 1
f_scale = std(f_in);
f_in = f_in / f_scale;

% Calculate mesh steps
x_steps = x_mesh(2:n) - x_mesh(1:n-1);
x_steps = reshape(x_steps, [1, n - 1]);
x_steps = [0, x_steps];
x_step_mean = mean(x_steps);
x_grad_mesh = (x_mesh(2:end) + x_mesh(1:end-1)) / 2;
x_mesh_grad_2 = (x_grad_mesh(2:end) + x_grad_mesh(1:end-1)) / 2;
% Rescale alpha
alpha = alpha * x_step_mean .^ 4;
w_0L = 2 * (x_mesh(n) - x_mesh(1));


%% Simple one-sided small finite difference matrix, (n-2) x (n-1)
D = zeros(n - 2, n - 1);
for i = 1:(n-2)
    step = (x_steps(i + 1) + x_steps(i + 2)) / 2;
    D(i, i:i+1) = [-1, 1] / step;
end;


%% Integration matrix A, n x (n - 1)
% The first row is zeros
A = zeros(n, n-1);
for i = 2:n
    A(i, :) = A(i-1, :);
    A(i, i-1) = x_steps(i);
end;


%% Boundary conditions matrix (subtract overall trend)
B = zeros(n-2, n-1);
B(:, 1) = - 1 / w_0L;
B(:, n - 1) = 1 / w_0L;


%% Make sparse
% D = sparse(D);
% A = sparse(A);
% B = sparse(B);


%% Optimize
Lagr = (D - B)' * (D - B);
% Hessian
H = A' * A + alpha .* Lagr ;
% Regularize the gradient
f_prime_reg = H \ (A' * f_in);
% Regularized function
f_reg = A * f_prime_reg;


%% Interpolate / extrapolate the gradient to the original bins centers
f_prime_reg_interpolated = interp1(x_grad_mesh, f_prime_reg, x_mesh, 'spline', 'extrap');


%% Calculate cost
cost_1 = 1/2 .* (f_reg - f_in) .^ 2;
cost_1 = trapz(x_mesh, cost_1);
cost_2 = 1/2 .* alpha .* ((D - B) * f_prime_reg) .^ 2;
cost_2 = trapz(x_mesh_grad_2, cost_2);
norm_cost = cost_1 + cost_2;


%% Restore 
% Restore scale
f_prime_reg = f_prime_reg .* f_scale;
f_prime_reg_interpolated = f_prime_reg_interpolated .* f_scale;
f_reg = f_reg .* f_scale;
% Restore shift
f_reg = f_reg + f_shift;




1;












