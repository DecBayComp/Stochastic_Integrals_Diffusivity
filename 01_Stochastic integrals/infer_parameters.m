

%% Constants
bin_width_multiplier = 10;
max_bins_number = 1000;

%% Initialization
D_map_inferred = cell(1, lambda_count);
f_map_inferred = cell(1, lambda_count);
x_bins_centers_saved = cell(1, lambda_count);


%% First binning my trajectory
for l_ind = 1:lambda_count
    x_bins_number = round(L / mean_jump_length(l_ind) / bin_width_multiplier);
    if x_bins_number > max_bins_number
        x_bins_number = max_bins_number;
    end;
    x_bin_width = L / x_bins_number;
    x_bins_borders = x_min + (0:x_bins_number) * x_bin_width;
    x_bins_centers = (x_bins_borders(1:end-1) + x_bins_borders(2:end))/2;


    %% Binning the data (not sure it's the most efficient algorithm)
    % Calculating the number of elements
    bin_numbers_raw = ceil((x_lambda(l_ind, :) - x_min) / x_bin_width);
    elements_in_bins_count = histcounts(bin_numbers_raw, 0.5 + (0:x_bins_number));
    % Initializing. I think using cells may be too slow. To check
    binned_jumps = cell(1, x_bins_number);
    for bin = 1:x_bins_number
        binned_jumps{bin} = zeros(2, elements_in_bins_count(bin));
    end;
    % Binning
    collected_points = zeros(1, x_bins_number);
    for i = 1:N
        bin = bin_numbers_raw(i);
        collected_points(bin) = collected_points(bin) + 1;
        binned_jumps{bin}(:, collected_points(bin)) = [x_lambda(l_ind, i), dx_lambda(l_ind, i)];
    end;
    % The data is binned in the form (x, dx) for the specified lambda


    %% Now I can infer diffusivity/force parameters in each zone

    D_map = zeros(1, x_bins_number);
    f_map = zeros(1, x_bins_number);
    for bin = 1:x_bins_number
        % Making the maximum likelihood estimates
        f_over_alpha_ML = mean(binned_jumps{bin}(2, :))/t_step;
        D_ML = var(binned_jumps{bin}(2, :)) / 2 /t_step;
        f_ML = f_over_alpha_ML * kBT / D_ML;
        % Saving
        D_map(bin) = D_ML;
        f_map(bin) = f_ML;
    end;
    
    D_map_inferred{l_ind} = D_map;
    f_map_inferred{l_ind} = f_map;
    x_bins_centers_saved{l_ind} = x_bins_centers;
end;


%% Plotting
plot_diffusivity_and_force;





