set (0, 'DefaultAxesFontSize', 20);



%% Constants
load_constants;
% bl_batch_start = 1;


%% Initialize
fD_min = [];
fD_max = [];
fD_step = [];
fD_mesh = [];
fD_pdf_plot_data = [];
pdf_norm = [];


% for D_case_number = 1:max_D_case_number
%     for f_case_number = 1:max_f_case_number
D_case_number = selected_D_case;
f_case_number = selected_f_c5ase;
        %% Initialize variables
        x_lambda = zeros(lambda_count, N);
        dx_lambda = zeros(lambda_count, N);
        for l_ind = 1:lambda_count
            fprintf('Processing case D: %i/%i, f: %i/%i, lambda: %i/%i\n',...
                D_case_number, max_D_case_number, f_case_number, max_f_case_number, l_ind, lambda_count);
            %% Load CSV file with x and dx data
            filename = sprintf('D_%i_f_%i_lambda_%.2f_trajectory.csv', D_case_number, f_case_number, lambda_array(l_ind));
            output_full_path = strcat(output_trajectories_folder, filename);
            
%             output_data = [x_lambda(l_ind, 1:end-1)', dx_lambda(l_ind, :)'];
            input_data = dlmread(output_full_path, CSV_DELIMITER);
            % Separate data into x and dx
            x_lambda(l_ind, 1:N) = input_data(1:N, 1);
            dx_lambda(l_ind, 1:N) = input_data(1:N, 2);
        
        end;
        input_data = [];
        %% Bin data
        bin_and_infer_bin_distributions;
        
        %% Calculate margnalized force
        fD_step = [];
        fD_mesh = [];
        fD_pdf_plot_data = [];
%         calculate_fD_posterior_with_lambda_in_selected_bins;
        
        %% Calculate missing variables
        t_mesh = (0:N) * t_step;
       
        %% Save results
        filename = sprintf('D_%i_f_%i_data.mat', D_case_number, f_case_number);
        output_full_path = strcat(output_data_folder, filename);

        fprintf('Saving processed data...\n');
        save(output_full_path, 'f_case_number', 'D_case_number', 'L', 'N', 'kBT', 't_step', ...
        'lambda_array', 't_mesh', 'x_lambda', 'dx_lambda', 'x_bins_centers_saved', ...
        'x_bins_borders_saved', 'x_bins_widths_saved', 'x_bins_number_saved', 'dx_mean_in_bins_saved', 'dx_std_in_bins_saved', 'n_j', 'dx_Mean', ...
        'V_j', 'V', 'fD_min', 'fD_max', 'fD_step', 'fD_mesh', 'fD_pdf_plot_data', 'pdf_norm', ...
        'elements_in_bins_count_saved', 'bl_use_adaptive_mesh');
        fprintf('Processed data saved successfully!\n');
        
        
        1;
%     end;
% end;









