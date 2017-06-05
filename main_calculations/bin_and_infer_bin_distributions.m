

%% Globals
global bl_use_adaptive_mesh;
global dx_Mean;
global elements_in_bins_count_saved;
% global min_points_in_bin;
global n_j;
global V;
global V_j;
global x_bins_centers_saved;
global x_bins_borders_saved;
global x_bins_widths_saved;
global x_bins_number_saved;
global dx_mean_in_bins_saved;
global dx_std_in_bins_saved;


%% Constants
load_constants;
initial_bin_width_multiplier = 1/8;



%% Initialize
x_bins_centers_saved = cell(1, lambda_count);
dx_mean_in_bins_saved = cell(1, lambda_count);
x_bins_widths_saved = cell(1, lambda_count);
dx_2_mean_in_bins_saved = cell(1, lambda_count);
dx_std_in_bins_saved = cell(1, lambda_count);
dx_Mean = zeros(1, lambda_count);
x_bins_number_saved = zeros(1, lambda_count);
V = zeros(1, lambda_count);
V_j = cell(1, lambda_count);
n_j = cell(1, lambda_count);
x_bin_width = zeros(1, lambda_count);
elements_in_bins_count = cell(1, lambda_count);


%% First binning my trajectory
% I will use the same mesh for the equilibrium distribution to facilitate
% the calculations
% If an empty bin is detected after binning, increase the bin width. Proceed until converged
for l_ind = 1:lambda_count
        
    if bl_use_adaptive_mesh
        fprintf('Binning progres (adaptive mesh): %i/%i\n', l_ind, lambda_count);
        [x_bins_borders, x_bins_centers, x_bins_number, x_bins_widths,...
            elements_in_bins_count, binned_jumps] = select_bins_adaptive_mesh(x_lambda(l_ind, :), dx_fwd_lambda(l_ind, :), min_points_in_bin_short_trajectories);
    else   
        fprintf('Binning progres (regular mesh): %i/%i\n', l_ind, lambda_count);
        [x_bins_borders, x_bins_centers, x_bins_number, x_bins_widths,...
            elements_in_bins_count, binned_jumps] = select_bins_equal_width(x_lambda(l_ind, :), dx_fwd_lambda(l_ind, :), min_points_in_bin_short_trajectories);
    end;
    % The data is binned in the form (x, dx) for the specified lambda
    
    
    %% Now I can infer diffusivity/force distributions in each bin
    dx_mean_in_bins = zeros(1, x_bins_number);
    dx_2_mean_in_bins = zeros(1, x_bins_number);
    dx_std_in_bins = zeros(1, x_bins_number);
    for bin = 1:x_bins_number
        dx_mean_in_bins(bin) = mean(binned_jumps{bin}(2, :));
        dx_2_mean_in_bins(bin) = mean(binned_jumps{bin}(2, :).^2);
        dx_std_in_bins(bin) = std(binned_jumps{bin}(2, :));
    end;
    V_in_bins = dx_2_mean_in_bins - dx_mean_in_bins.^2;
    % And the number of jumps per bin is stored in elements_in_bins_count

    
    %% Saving data for current lambda
    % Bin parameters
    elements_in_bins_count_saved{l_ind} = elements_in_bins_count;
    x_bins_centers_saved{l_ind} = x_bins_centers;
    x_bins_borders_saved{l_ind} = x_bins_borders;
    x_bins_widths_saved{l_ind} = x_bins_widths;
    x_bins_number_saved(l_ind) = x_bins_number;
    % Trajectory parameters
    dx_mean_in_bins_saved{l_ind} = dx_mean_in_bins;
    dx_2_mean_in_bins_saved{l_ind} = dx_2_mean_in_bins;
    dx_std_in_bins_saved{l_ind} = dx_std_in_bins;
    V_j{l_ind} = V_in_bins;
    n_j{l_ind} = elements_in_bins_count;
    dx_Mean(l_ind) = mean(dx_fwd_lambda(l_ind, :));
    V(l_ind) = mean(dx_fwd_lambda(l_ind, :).^2) - dx_Mean(l_ind)^2;
    
end;

l_ind = 1;
fprintf('Progres: Finished\n\n');
fprintf('Lambda = %.2f. Total points: %i\nBins number: %i\nMin points in bin: %i\nMax points in bin: %i\n',...
    lambda_array(l_ind), N, x_bins_number_saved(l_ind), min(min(elements_in_bins_count_saved{l_ind})), max(max(elements_in_bins_count_saved{l_ind})));











