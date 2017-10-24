%% Calculate regularized graident based on a lattice-defined input function.
% The parameter alpha defines smoothing degree.
%
% It is assumed some of the input values may be NaN. In this cases, the values are omitted for the calculations, but the results are interpolated to the meshes
% expected if there were no NaN values
%
% 


function [f_reg, f_prime_reg, f_prime_reg_interpolated, norm_cost, x_grad_mesh] = regularize_gradient(f_in, x_mesh, alpha)

%% Initialize
% Make column (vector)
f_in = reshape(f_in, [length(f_in), 1]);



%% Extract non-Nan points

non_nan_indices = find(~isnan(f_in));
% Extract non NaN positions in the x mesh
x_mesh_full = x_mesh;
x_mesh = x_mesh(~isnan(f_in));

% Modify f_in
f_in_full = f_in;
f_in = f_in(~isnan(f_in));

% Determine data length
n_full = length(f_in_full);
n = length(f_in);



%% Calculate mesh
x_steps = x_mesh(2:n) - x_mesh(1:n-1);
x_steps = reshape(x_steps, [1, n - 1]);

% Add 0 in front of the steps
x_steps = [0, x_steps];
x_step_mean = mean(x_steps);

% First derivative mesh (in-between the original points)
x_grad_mesh = (x_mesh(2:end) + x_mesh(1:end-1)) / 2;

% Second derivative mesh (in-between first gradient points)
x_mesh_grad_2 = (x_grad_mesh(2:end) + x_grad_mesh(1:end-1)) / 2;

% Create the gradient mesh for the full x mesh
x_grad_mesh_full = (x_mesh_full(2:end) + x_mesh_full(1:end-1)) / 2;



%% Pre-process data
% Shift the left f value to 0
f_shift = f_in(1);
f_in = f_in - f_shift;

% Normalize std to 1
f_scale = std(f_in);
f_in = f_in / f_scale;

% Rescale alpha
alpha = alpha * x_step_mean .^ 4;
w_0L = 2 * (x_mesh(n) - x_mesh(1));



%% Calculate

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



%% Calculate gradients
Lagr = (D - B)' * (D - B);
% Hessian
H = A' * A + alpha .* Lagr ;
% Regularize the gradient
f_prime_reg = H \ (A' * f_in);
% Regularized function
f_reg = A * f_prime_reg;



%% Interpolate / extrapolate
% Interpolate the gradient to the original bins centers
f_prime_reg_interpolated = interp1(x_grad_mesh, f_prime_reg, x_mesh_full, 'spline', 'extrap');

% Interpolate the gradient to the full gradient mesh
f_prime_reg_full = interp1(x_grad_mesh, f_prime_reg, x_grad_mesh_full, 'spline', 'extrap');

% Interpolate the regularized function to the original mesh
f_reg_full = interp1(x_mesh, f_reg, x_mesh_full, 'spline', 'extrap');


%% Calculate cost
cost_1 = 1/2 .* (f_reg - f_in) .^ 2;
cost_1 = trapz(x_mesh, cost_1);
cost_2 = 1/2 .* alpha .* ((D - B) * f_prime_reg) .^ 2;
cost_2 = trapz(x_mesh_grad_2, cost_2);
norm_cost = cost_1 + cost_2;


%% Restore 
% Restore scale

% f'
f_prime_reg = f_prime_reg .* f_scale;
f_prime_reg_full = f_prime_reg_full .* f_scale;
f_prime_reg_interpolated = f_prime_reg_interpolated .* f_scale;

% f
f_reg = f_reg .* f_scale;
f_reg_full = f_reg_full .* f_scale;

% Restore shift
f_reg = f_reg + f_shift;
f_reg_full = f_reg_full + f_shift;

% % % Restore NaN values in f_reg
% % f_reg_new = zeros(1, length(x_mesh_full)) * NaN;
% % f_reg_new(non_nan_indices) = f_reg;
% % f_reg = f_reg_new;
% % % for i = 1:length(x_mesh_orig)
% % % 	f_reg_new(i) = f_reg(i);
% % % end;


%% Final assignments to return *regularized values* on the *original mesh*
f_reg = f_reg_full;
f_prime_reg = f_prime_reg_full;
x_grad_mesh = x_grad_mesh_full;


1;












