


%% Constants
load_constants;
x_lim_vec = [0, 0.41];
trials = 10;


l_ind = 3;

% Load bin structure
x_bins_number = data_struct.x_bins_number(l_ind);
x_bins_centers = data_struct.x_bins_centers_saved{l_ind};

% Initialize
D_KS_Ito = cell(1, lambda_count);
D_KS_Ito_mean = cell(1, lambda_count);
D_KS_Stratonovich = cell(1, lambda_count);
D_KS_Stratonovich_mean = cell(1, lambda_count);
D_KS_divine = cell(1, lambda_count);
D_KS_divine_mean = cell(1, lambda_count);
D_KS_marginalized = cell(1, lambda_count);
D_KS_marginalized_mean = cell(1, lambda_count);
for l_ind1 = 1:lambda_count
    D_KS_Ito{l_ind1} = zeros(trials, x_bins_number);
    D_KS_Stratonovich{l_ind1} = zeros(trials, x_bins_number);
    D_KS_divine{l_ind1} = zeros(trials, x_bins_number);
    D_KS_marginalized{l_ind1} = zeros(trials, x_bins_number);
end;

% Calculate exact fD
fD_exact = f_func(selected_f_case, x_bins_centers * L, L)...
    .* D_func(selected_D_case, x_bins_centers * L, L);

% Collect the MAP force distribution in all bins over trials
MAP_fD_data_Ito = zeros(trials, x_bins_number);
MAP_fD_data_Stratonovich = zeros(trials, x_bins_number);
MAP_fD_data_divine = zeros(trials, x_bins_number);
MAP_fD_data_marginalized = zeros(trials, x_bins_number);
for trial = 1:trials
    MAP_fD_data_Ito(trial, :) = trials_data{trial}.MAP_fD_Ito{l_ind}(1, :);
    MAP_fD_data_Stratonovich(trial, :) = trials_data{trial}.MAP_fwd_fD_Stratonovich{l_ind}(1, :);
    MAP_fD_data_divine(trial, :) = trials_data{trial}.MAP_fwd_fD_divine{l_ind}(1, :);
    MAP_fD_data_marginalized(trial, :) = trials_data{trial}.MAP_fwd_fD_marginalized{l_ind}(1, :);
end;


% Select bin
Ito_results = zeros(trials, x_bins_number);
Stratonovich_results = zeros(trials, x_bins_number);
divine_results = zeros(trials, x_bins_number);
marginalized_results = zeros(trials, x_bins_number);
% count_completed = 0;
% count_total = trials * x_bins_number;
fprintf('Starting K-S distance calculations.\n');
tic; 
parfor trial = 1:trials
%     fprintf('Calculating K-S distance. Lambda: %i/%i. Trial: %i/%i\n', l_ind, lambda_count, trial, trials);
    for bin = 1:x_bins_number
        fprintf('Calculating K-S distance. Lambda: %i/%i. Trial: %i/%i. Bin: %i/%i.\n', l_ind, lambda_count, trial, trials, bin, x_bins_number);
        %% Ito distribution
        emp_distrib = fD_exact(bin) - MAP_fD_data_Ito(:, bin);
        % Get posterior
        cont_distrib_func = @(f) bin_fD_posterior_func (trials_data{trial}, l_ind, bin, f + MAP_fD_data_Ito(trial, bin), 'forward'); 
        % Calculate distance
        Ito_results(trial, bin) = calculate_KS_distance(emp_distrib, cont_distrib_func);     
        
        %% Stratonovich
        emp_distrib = fD_exact(bin) - MAP_fD_data_Stratonovich(:, bin);
        % Get posterior
        cont_distrib_func = @(f) bin_fD_simple_Stratonovich_posterior_func(trials_data{trial}, ...
                l_ind, bin,  f + MAP_fD_data_Stratonovich(trial, bin), trials_data{trial}.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
        % Calculate distance
        Stratonovich_results(trial, bin) = calculate_KS_distance(emp_distrib, cont_distrib_func);     
        
        %% Divine
        emp_distrib = fD_exact(bin) - MAP_fD_data_divine(:, bin);
        % Get posterior
        cont_distrib_func = @(f) bin_fD_divine_inference_posterior_func(trials_data{trial}, ...
            l_ind, bin, f + MAP_fD_data_divine(trial, bin), trials_data{trial}.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
        % Calculate distance
        divine_results(trial, bin) = calculate_KS_distance(emp_distrib, cont_distrib_func); 
        
        %% Marginalized
        emp_distrib = fD_exact(bin) - MAP_fD_data_marginalized(:, bin);
        % Get posterior
        cont_distrib_func = @(f) bin_fD_lambda_marginalized_posterior_func(trials_data{trial}, ...
                l_ind, bin, f + MAP_fD_data_marginalized(trial, bin), trials_data{trial}.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
        % Calculate distance
        marginalized_results(trial, bin) = calculate_KS_distance(emp_distrib, cont_distrib_func); 
%         count_completed = count_completed + 1;
    end
end;
D_KS_Ito{l_ind} = Ito_results;
D_KS_Stratonovich{l_ind} = Stratonovich_results;
D_KS_divine{l_ind} = divine_results;
D_KS_marginalized{l_ind} = marginalized_results;

fprintf('Calcualtions completed in %.2f s\n', toc);

% Calculate mean
D_KS_Ito_mean{l_ind} = mean(D_KS_Ito{l_ind}, 1);
D_KS_Stratonovich_mean{l_ind} = mean(D_KS_Stratonovich{l_ind}, 1);
D_KS_divine_mean{l_ind} = mean(D_KS_divine{l_ind}, 1);
D_KS_marginalized_mean{l_ind} = mean(D_KS_marginalized{l_ind}, 1);


%% Plot
figure(7);
clf;
hold on;
% Divine
plot(x_bins_centers, D_KS_divine_mean{l_ind}, strcat('-', markers_list{enum_divine}), 'color', color_sequence(enum_divine, :),...
    'markers', marker_size, 'LineWidth', line_width);
% Ito
plot(x_bins_centers, D_KS_Ito_mean{l_ind}, strcat('-', markers_list{enum_Ito}), 'color', color_sequence(enum_Ito, :),...
    'markers', marker_size, 'LineWidth', line_width);
% Stratonovich
plot(x_bins_centers, D_KS_Stratonovich_mean{l_ind}, strcat('-', markers_list{enum_Stratonovich}), 'color', color_sequence(enum_Stratonovich, :),...
    'markers', marker_size, 'LineWidth', line_width);
% Marginalized
plot(x_bins_centers, D_KS_marginalized_mean{l_ind}, strcat('-', markers_list{enum_marginalized}), 'color', color_sequence(enum_marginalized, :),...
    'markers', marker_size, 'LineWidth', line_width);
% Adjust
xlim(x_lim_vec);
ylim([0, 1]);
xlabel('x');
ylabel('K-S distance');
title(sprintf('$\\lambda^* = %.2f$', lambda_array(l_ind)), 'interpreter', 'latex');
str_legend = {'Orcl', 'Ito', 'Str', 'Mar'};
legend(str_legend, 'location', 'southwest');

















% 
% figure(8);
% clf;
% hold on;
% emp_distrib = sort(emp_distrib);
% plot(emp_distrib, 1:trials);
% plot(emp_distrib, cum_cont_distrib);
% 
% 



% % Adjust
% title(sprintf('Posterior for $x\\approx%.2f$, $\\lambda^* = %.1f$', selected_x_over_L, lambda_array(l_ind)), 'interpreter', 'latex');
% % xlim(x_lim_vec);
% % ylim(y_lim_vec);
% xlabel('fD', 'interpreter', 'latex');
% ylabel('PDF', 'interpreter', 'latex');
% % % Legend
% % str_legend_local = {'Orcl', 'Ito', 'Str', 'Mar'};
% % legend(str_legend_local, 'location', 'northwest', 'interpreter', 'latex', 'fontsize', font_size);
% % 
% % 
% 





% % % %% === PLOT ===
% % % %% Plot force posterior
% % % figure(7);
% % % clf;
% % % hold on;


% % % % Detect the longest interval and centrally increase its length
% % % fD_min = min(MAP_fD_data_Ito_bin );
% % % fD_max = max(MAP_fD_data_Ito_bin);
% % % fD_interval = fD_max - fD_min;
% % % fD_min = fD_min - fD_interval * (factor-1)/2;
% % % fD_max = fD_max + fD_interval * (factor-1)/2;
% % % % fD_min = -0.4;
% % % % fD_max = 0.4;
% % % % Prepare mesh
% % % fD_step = (fD_max - fD_min) / (fD_steps - 1);
% % % fD_mesh = fD_min:fD_step:fD_max;
% % % % Add special points
% % % fD_mesh = sort([fD_mesh, fD_values]);
% % % % Calculate data to plot
% % % size1 = lambda_count;
% % % size2 = length(fD_mesh);
% % % fD_pdf_plot_data_divine = zeros(size1, size2);
% % % fD_pdf_plot_data_Ito = zeros(size1, size2);
% % % fD_pdf_plot_data_Stratonovich = zeros(size1, size2);
% % % fD_pdf_plot_data_marginalized = zeros(size1, size2);
% % % % fD posterior data
% % % trial = trial_selected;
% % % for l_ind = 1:lambda_count
% % %     bin = selected_bins_indices(l_ind);
% % %     % Divine
% % %     fD_pdf_plot_data_divine(l_ind, :) = bin_fD_divine_inference_posterior_func(trials_data{trial}, ...
% % %         l_ind, bin, fD_mesh, trials_data{trial}.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
% % %     % Ito
% % %     fD_pdf_plot_data_Ito(l_ind, :) = bin_fD_posterior_func (trials_data{trial}, l_ind, bin, fD_mesh, 'forward');    
% % %     % Stratonovich
% % %     fD_pdf_plot_data_Stratonovich(l_ind, :) = bin_fD_simple_Stratonovich_posterior_func(trials_data{trial}, ...
% % %                 l_ind, bin, fD_mesh, trials_data{trial}.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
% % %     % Marginalized
% % %     fD_pdf_plot_data_marginalized(l_ind, :) = bin_fD_lambda_marginalized_posterior_func(trials_data{trial}, ...
% % %                 l_ind, bin, fD_mesh, trials_data{trial}.MAP_fwd_D_grad_regular_interp{l_ind}(bin), 'forward');
% % % end;
% % % % Plot
% % % % % Divine
% % % % plot(fD_mesh, fD_pdf_plot_data_divine(l_ind, :), strcat('-', markers_list{1}), 'LineWidth', line_width/2);
% % % 
% % % % Ito
% % % plot(fD_mesh - MAP_fD_data_Ito_bin(trial), fD_pdf_plot_data_Ito(l_ind, :), strcat('-', markers_list{2}), 'LineWidth', line_width/2);
% % % % plot(fD_mesh, fD_pdf_plot_data_Ito(l_ind, :), strcat('-', markers_list{2}), 'LineWidth', line_width/2);
% % % 
% % % % % Stratonovich
% % % % plot(fD_mesh, fD_pdf_plot_data_Stratonovich(l_ind, :), strcat('-', markers_list{3}), 'LineWidth', line_width/2);
% % % % % Marginalized
% % % % plot(fD_mesh, fD_pdf_plot_data_marginalized(l_ind, :), '-', 'LineWidth', line_width);
% % % 
% % % % % True value
% % % % y_lim_vec = [0, max([fD_pdf_plot_data_marginalized(l_ind, :), fD_pdf_plot_data_Ito(l_ind, :)]) * 1.05];
% % % % plot(fD_values(2) .* [1, 1], y_lim_vec, '--k', 'LineWidth', line_width);
% % % 
% % % % ylim(y_lim_vec);    % I don't know why it changes scale without it
% % % 
% % % 
% % % %% Plot MAP histogram
% % % MAP_fD_bias_Ito_bin = fD_values(2) - MAP_fD_data_Ito_bin;
% % % % MAP_fD_bias_Ito_bin = MAP_fD_data_Ito_bin;
% % % h_hist = histogram(MAP_fD_bias_Ito_bin, 'Normalization', 'pdf');
% % % uistack(h_hist, 'bottom');




