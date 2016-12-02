set (0, 'DefaultAxesFontSize', 20);


%% Constants
m = 1;
kBT = 1;


L = 40;
x_min = -L/2;
x_max = L/2;
% % L = x_max - x_min;
% bl_periodic = true;
% T = 100000;
t_step = 0.01;
N = 1e7;
T = t_step * N;
update_progress_every = 1e4;

% Alpha func
alpha_func = @(x) 2.75 + 2.25 * sin(2*pi*x/L);
alpha_prime_func = @(x) 2.25 * 2*pi/L * cos(2*pi*x/L);

% f func
f_func = 0; % = @(x) 0*x;

% Choosing the working mode
% str_mode = 'periodic';
str_mode = 'inf_walls';

% % Initializing mode-specific parameters
% if strcmp(str_mode, 'periodic')
%     bl_periodic = true;
% elseif strcmp(str_mode, 'inf_walls')
%     



%% Initialization
rng('shuffle');
lambda_array = [0, 0.5, 1];
lambda_names_array = {'Itô', 'Stratonovich', 'Isothermal'}; 
lambda_count = length(lambda_array);
x_lambda = zeros(lambda_count, N+1);
v_lambda = zeros(lambda_count, N+1);
t_mesh = (0:N) * t_step;

% Choosing the first point randomly from the interval
x_0 = x_min + rand * L;
% The initial velocity is taken randomly from Maxwell's distribution
v_0 = kBT/2/m * randn;

x_lambda(:, 1) = x_0;
v_lambda(:, 1) = v_0;


%% Using Verlet method to make iterations.
% We have to choose the convention for noise also that we wish to choose.
% For the case of constant diffusivity there will be no difference
for i = 1:N
    q = randn;  % Using the same random white noise for all constructs
    for l_ind = 1:lambda_count
        lambda = lambda_array(l_ind);

        %% Extracting current x and v
        x_i = x_lambda(l_ind, i);
        v_i = v_lambda(l_ind, i);
        alpha_i = alpha_func(x_i);

        %% Creating noise W
        dW = sqrt(t_step) * q;
        
        %% Calculating c
        c = 1/(1 + alpha_i * t_step / (2 * m));
        
        %% Calculating increments
        
        dx = c * t_step * (v_i + (f_func + lambda * kBT / m^2 * alpha_prime_func(x_i)) * t_step/2 ...
            + sqrt(2 * kBT * alpha_i) / (2 * m) * dW);
                
        v_next = 2 * dx / t_step - v_i;
        x_next = x_i + dx;


        %% Taking into account the BCs
        if strcmp(str_mode, 'periodic')
            if x_next > x_max
                x_next = x_next - L;
            elseif x_next < x_min
                x_next = x_next + L;
            end;
        elseif strcmp(str_mode, 'inf_walls')
            if x_next > x_max
                x_next = 2 * x_max - x_next;
            elseif x_next < x_min
                x_next = 2 * x_min - x_next;
            end;
        else
            disp('Error! Wrong mode selected!');
        end;
        
        %% Saving
        x_lambda(l_ind, i+1) = x_next;
        v_lambda(l_ind, i+1) = v_next;
    
    % Iterating lambdas
    end;
    
    %% Printing out the simulation progress
    if mod (i, update_progress_every) == 0
        fprintf('Simulation progress: %.1f %%\n', i/N*100);
    end;
    
    % Iterating
end;




%% Plotting
plot_trajectory;
plot_velocity;
plot_diffusivity;
plot_density;








