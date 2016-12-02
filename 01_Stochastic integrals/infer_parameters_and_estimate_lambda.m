

%% Constants
bin_width_multiplier = 10;
max_bins_number = 1000;

%% Initialization
D_map_inferred = cell(1, lambda_count);
b_map_inferred = cell(1, lambda_count);
a_map_inferred = cell(1, lambda_count);
f_Pb2_map_inferred = cell(1, lambda_count);
f_lambda_map_inferred = cell(1, lambda_count);
x_bins_centers_saved = cell(1, lambda_count);


%% First binning my trajectory
% I will use the same mesh for the equilibrium distribution to facilitate
% the calculations
for l_ind = 1:lambda_count
    x_bins_number = round(L / mean_jump_length(l_ind) / bin_width_multiplier);
    if x_bins_number > max_bins_number
        x_bins_number = max_bins_number;
    end;
    x_bin_width = L / x_bins_number;
    x_bins_borders = x_min + (0:x_bins_number) * x_bin_width;
    x_bins_centers = (x_bins_borders(1:end-1) + x_bins_borders(2:end))/2;


    %% Binning the trajectory 
    % Calculating the number of elements
    bin_occupation_numbers_raw = ceil((x_lambda(l_ind, :) - x_min) / x_bin_width);
    elements_in_bins_count = histcounts(bin_occupation_numbers_raw, 0.5 + (0:x_bins_number));
    % Initializing. I think using cells may be too slow. To check
    binned_jumps = cell(1, x_bins_number);
    for bin = 1:x_bins_number
        binned_jumps{bin} = zeros(2, elements_in_bins_count(bin));
    end;
    % Binning
    collected_points = zeros(1, x_bins_number);
    for i = 1:N
        bin = bin_occupation_numbers_raw(i);
        collected_points(bin) = collected_points(bin) + 1;
        binned_jumps{bin}(:, collected_points(bin)) = [x_lambda(l_ind, i), dx_lambda(l_ind, i)];
    end;
    % The data is binned in the form (x, dx) for the specified lambda


    %% Now I can infer diffusivity/force parameters in each bin

    D_map = zeros(1, x_bins_number);
    a_map = zeros(1, x_bins_number);
    for bin = 1:x_bins_number
        % Making the maximum likelihood estimates
        f_over_alpha_ML = mean(binned_jumps{bin}(2, :))/t_step;
        D_ML = var(binned_jumps{bin}(2, :)) / 2 /t_step;
        f_ML = f_over_alpha_ML * kBT / D_ML;
        % Saving
        D_map(bin) = D_ML;
        a_map(bin) = f_ML;
    end;
    
    % Calculating b from D
    b_map = sqrt(2 * D_map);
    
    D_map_inferred{l_ind} = D_map;
    b_map_inferred{l_ind} = b_map;
    a_map_inferred{l_ind} = a_map;
    x_bins_centers_saved{l_ind} = x_bins_centers;
    
    %% Now working with the equilibrium distribution
    % Binning normalizing the sum to 1
    P_equilibrium = histcounts(x_lambda(l_ind, :), x_bins_borders, 'Normalization', 'probability');
    
    %% Calculating different forces
    % P * b^2
    Pb2_log = log(P_equilibrium .* (b_map_inferred{l_ind} .^2));
   
    
    
    % Constructing a spline interpolation
    Pb2_log_spline = spline(x_bins_centers, Pb2_log);
    % Numerically differentiating with splines
    Pb2_log_derivative = smooth_derivative(Pb2_log_spline);
    
    % Plot
    figure(5);
    clf;
    plot(x_bins_centers, Pb2_log, 'LineWidth', 2);
    hold on;
    plot(x_bins_centers, ppval(Pb2_log_derivative, x_bins_centers), 'LineWidth', 2);
    legend('ln(Pb^2)', '\partial_x ln(Pb^2)', 'Location', 'best');
    
    % Plotting to check the calculations
% %     deriv_values = ppval(Pb2_log_derivative, x_mesh_temp);
% %     plot(x_mesh_temp, deriv_values, 'b', 'LineWidth', 2);
     
    
    % Now calculating the force component that can be inferred from Pb^2
    f_Pb2_map = kBT * ppval(Pb2_log_derivative, x_bins_centers);
    
    % Repeating the same calculations for the spurious component of the
    % force
    log_b = log(b_map_inferred{l_ind});
    grad_log_b = smooth_derivative(spline(x_bins_centers, log_b));
    % Sampling the curve at bins centers
    f_lambda_map = kBT * 2 * ppval(grad_log_b, x_bins_centers);
    
% %     %% Plotting the two forces
% %     figure(3);
% %     clf;
% %     plot(x_bins_centers, f_Pb2_map, 'b', 'LineWidth', 2);
% %     hold on;
% %     plot(x_bins_centers, f_Pb2_map - f_lambda_map, 'r', 'LineWidth', 2);
% %     legend('\lambda=0', '\lambda=1', 'Location', 'best');
    
    % Saving force values
    f_lambda_map_inferred{l_ind} = f_lambda_map;
    f_Pb2_map_inferred{l_ind} = f_Pb2_map;
    
    
    1;
end;


%% Plotting
% plot_diffusivity_and_force;





