set (0, 'DefaultAxesFontSize', 20);

%% All globals are loaded here to be initialized. Other files use only globals they need
global bl_simulation_succeeded;
% global bl_use_adaptive_mesh;
global kBT;
global L;
% global min_points_in_bin;
global N;
% global n_j;
global t_step;
% global SIGMA_MIN;
% global dx_Mean;
% global V;
% global V_j;
% global ABS_TOLERANCE;
% global REL_TOLERANCE;
global lambda_array;
global lambda_count;
global selected_x_over_L;
global max_D_case_number;
global D_case_number;
global max_f_case_number;
global f_case_number;
% global fD_marginalized_steps;
% global one_column_figure_width_in;
% global two_columns_figure_width_in;
% global mac_screen_dpi;
% global font_size;
% global x_bins_centers_saved;
% global dx_mean_in_bins_saved;


%% Constants
load_constants;


%% Choose the regime of D and f
% D_case_number = 1;  % Parabolic diffusivity
% D_case_number = 2;  % Discontinuity in D 
% D_case_number = 3;  % Exponentially decreasing diffusivity
% D_case_number = 4;

% f_case_number = 1;  % f0 = 0
% f_case_number = 2;  % f0 = 1
% f_case_number = 3; 
% f_case_number = 4;  % Linearly growing force from f = -1 to f = 1
% f_case_number = 5;  % Linearly decreasing force from f = 1 to f = -1



% % Initializing mode-specific parameters
% if strcmp(str_mode, 'periodic')
%     bl_periodic = true;
% elseif strcmp(str_mode, 'inf_walls')
%     



%% Initialization
t_step_internal = t_step / internal_steps_number;
N_internal = N * internal_steps_number;
rng('shuffle');
x_lambda = zeros(lambda_count, N+1);
dx_lambda = zeros(lambda_count, N);
% v_lambda = zeros(lambda_count, N+1);
t_mesh = (0:N) * t_step;

% % % Choosing the first point randomly from the interval
% % x_0 = x_min + rand * L;
% Choosing the first point manually
x_0 = selected_x_over_L * L;
% The initial velocity is taken randomly from Maxwell's distribution
% v_0 = kBT/2/m * randn;

x_lambda(:, 1) = x_0;
% v_lambda(:, 1) = v_0;


%% Using Verlet method to make iterations.
% We have to choose the convention for noise also that we wish to choose.
% For the case of constant diffusivity there will be no difference
for i = 1:N
    q = randn(1, internal_steps_number);  % Using the same random white noise for all constructs
    for l_ind = 1:lambda_count
        lambda = lambda_array(l_ind);

        %% Extracting current x and v
        x_i = x_lambda(l_ind, i);
%         b_prime_b_i = D_prime_func(x_i);
        for m = 1:internal_steps_number
            %% Calculate f and D
            f_i = f_func(f_case_number, x_i, L);
            [D_i, b_prime_b_i] = D_func(D_case_number, x_i, L);
            a_lambda_i = f_i * D_i / kBT;
            b_i = sqrt(2 * D_i);
            
            %% Creating the noise increment W_n
            dW = sqrt(t_step_internal) * q(m);  
            
            %% Calculating increments
            dx = a_lambda_i * t_step_internal ...
                + b_i * dW ...
                + 1/2 * b_prime_b_i * (dW^2 + (2 * lambda - 1) * t_step_internal);

            x_next = x_i + dx;
                        
            %% Taking into account the BCs
            switch bc_type
                case ENUM_BC_PERIODIC
                    if x_next > x_max
                        x_next = x_next - L;
                    elseif x_next < x_min
                        x_next = x_next + L;
                    end;
                case ENUM_BC_INF_WALLS 
                    if x_next > x_max
                        x_next = 2 * x_max - x_next;
                    elseif x_next < x_min
                        x_next = 2 * x_min - x_next;
                    end;
                otherwise
                    disp('Error! Wrong mode selected!');
            end;
            % Save the starting point for the next internal round
            x_i = x_next;   
        
        end;
        
        %% Saving
        x_lambda(l_ind, i+1) = x_next;
        dx_lambda(l_ind, i) = x_lambda(l_ind, i+1) - x_lambda(l_ind, i);
    
    % Iterating lambdas
    end;
    
    %% Printing out the simulation progress
    if mod (i, update_progress_every) == 0
        fprintf('D case: %i/%i. f case: %i/%i. Simulation progress: %.3f %%\n',...
            D_case_number, max_D_case_number, f_case_number, max_f_case_number, i/N*100);
    end;
    
    % Iterating
end;


%% Extracting jump lengths and calculating mean
% dx_lambda = x_lambda(:, 2:end) - x_lambda(:, 1:end-1);
mean_jump_length = mean(abs(dx_lambda), 2);


%% Plotting

% bin_and_infer_bin_distributions;

%% Check if the simulation succeeded and the bin we are interested in is not the last bin
bl_simulation_succeeded = true;


% plot_diffusivity;
% plot_force;
% plot_trajectory;
% plot_density;
% 
% plot_D_posterior_in_selected_bins;
% plot_D_grad_posterior_in_selected_bins;
% plot_fD_posterior_in_selected_bins;

% plot_fD_posterior_with_lambda_in_selected_bins;


%% Saving data


if bl_save_data
    fprintf('Saving trajectory only...\n');
    for l_ind = 1:lambda_count
        filename = sprintf('D_%i_f_%i_lambda_%.2f_trajectory.csv', D_case_number, f_case_number, lambda_array(l_ind));
        output_full_path = strcat(output_trajectories_folder, filename);
        output_data = [x_lambda(l_ind, 1:end-1)', dx_lambda(l_ind, :)'];
        dlmwrite(output_full_path, output_data, 'delimiter', CSV_DELIMITER);
    %     save(output_full_path, 'f_case_number', 'D_case_number', 'L', 'N', 'kBT', 't_step', ...
    %     'lambda_array', 't_mesh', 'x_lambda', 'dx_lambda', 'x_bins_centers_saved', ...
    %     'x_bins_borders_saved', 'x_bins_widths_saved', 'x_bins_number_saved', 'dx_mean_in_bins_saved', 'dx_std_in_bins_saved', 'n_j', 'dx_Mean', ...
    %     'V_j', 'V', 'fD_min', 'fD_max', 'fD_step', 'fD_mesh', 'fD_pdf_plot_data', 'pdf_norm', ...
    %     'elements_in_bins_count_saved', 'bl_use_adaptive_mesh');
    end;
    fprintf('Trajectory saved successfully!\n');
end;








