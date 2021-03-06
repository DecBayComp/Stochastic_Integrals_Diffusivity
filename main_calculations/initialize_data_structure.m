

function data_struct = initialize_data_structure(bins_count, fine_mesh_steps_count, conventions_count, lambda_types_count)

%% Comments
% - bins_count is a vector of length lambda_count

%% Constants
% Bin level
s1 = 4;


%%
data_struct = struct();

% Lambda level
% data_struct.elements_in_bins_count = cell(1, lambda_count);
data_struct.dx_mean_all_bins_all_trials = 0;
data_struct.dx_mean_all_bins = 0;
data_struct.dx_mean_in_bins = zeros(1, bins_count);
data_struct.mean_jump_length_bins = zeros(1, bins_count);
data_struct.mean_jump_bins_all_trials = zeros(1, bins_count);
data_struct.mean_log_K_L_avg_neg_x = 0;
data_struct.mean_log_K_L_avg_pos_x = 0;
% data_struct.dx_bck_mean_in_bins_saved = cell(1, lambda_count);
% data_struct.points_in_bins = cell(1, max_bins_number);
% data_struct.MAP_fD_Hanggi = cell(1, lambda_count);
% data_struct.MAP_fD_Hanggi_mean = cell(1, lambda_count);
data_struct.MAP_b = zeros(bins_count, s1);
data_struct.MAP_D = zeros(bins_count, s1);
data_struct.MAP_b_regular = zeros(1, bins_count);
data_struct.MAP_bb_prime_regular = zeros(1, bins_count - 1);
data_struct.MAP_bb_prime_regular_interp = zeros(1, bins_count);
data_struct.MAP_bb_prime_regular_interp_mean = 0;
% data_struct.MAP_fwd_fD_divine = cell(1, lambda_count);
% data_struct.MAP_fwd_fD_divine_mean = cell(1, lambda_count);
data_struct.MAP_b_mean = 0;
data_struct.MAP_D_mean = 0;
data_struct.b_theor_data = zeros(bins_count, 3);							% At bins centers, b and its two first derivatives (??)
% data_struct.D_theor_data = zeros(bins_count, 3);                          % At bins centers, D and its two first derivatives
data_struct.D_theor_fine_data = zeros(1, fine_mesh_steps_count);
data_struct.b_theor_fine_data = zeros(1, fine_mesh_steps_count);
data_struct.bb_prime_theor_fine_data = zeros(1, fine_mesh_steps_count);
data_struct.a_theor_data = zeros(1, bins_count);                          % At bins centers
data_struct.fD_theor_fine_data = zeros(1, fine_mesh_steps_count);           % On a fine mesh
data_struct.MAP_a = zeros(bins_count, conventions_count, s1);
data_struct.MAP_a_mean = 0;
data_struct.bl_empty_bin = zeros(1, bins_count);
% data_struct.MAP_fD_Ito_widths = cell(1, lambda_count);
% data_struct.MAP_fwd_fD_marginalized = cell(1, lambda_count);
% data_struct.MAP_fwd_fD_marginalized_mean = cell(1, lambda_count);
% data_struct.MAP_fwd_fD_Stratonovich = cell(1, lambda_count);
% data_struct.MAP_fwd_fD_Stratonovich_mean = cell(1, lambda_count);
% data_struct.MAP_fD_mean_sigma2_Ito = cell(1, lambda_count);
% data_struct.MAP_fD_mean_sigma2_Stratonovich = cell(1, lambda_count);
% data_struct.MAP_fD_mean_sigma2_marginalized = cell(1, lambda_count);
% data_struct.MAP_fD_mean_sigma2_divine = cell(1, lambda_count);
% data_struct.fD_theor = cell(1, lambda_count);
data_struct.lambda = -1;
data_struct.n_j = zeros(1, bins_count);
data_struct.n_j_mean = zeros(1, bins_count);
data_struct.trials_MAP_b = 0;
data_struct.trials_MAP_a = 0;
data_struct.trials_MAP_bb_prime_regular_interp = 0;
data_struct.trials_b_KS_distance = 0;
data_struct.trials_b_KS_distance_mean = zeros(lambda_types_count, bins_count);
data_struct.trials_log_K_L = 0;
data_struct.trials_log_K_G = 0;
data_struct.UR_b = 0;
data_struct.UR_a = 0;
data_struct.UR_a_bin_mean = 0;
data_struct.UR_a_bin_max = 0;
% data_struct.UR_D = cell(1, lambda_count);
% data_struct.UR_fwd_divine = cell(1, lambda_count);
% data_struct.UR_Ito = cell(1, lambda_count);
% data_struct.UR_fwd_marginalized = cell(1, lambda_count);
% data_struct.UR_fwd_Stratonovich = cell(1, lambda_count);
% data_struct.UR_Hanggi = cell(1, lambda_count);
data_struct.V_all_bins = 0;
data_struct.V_j = zeros(1, bins_count);
% data_struct.V_bck_j = cell(1, lambda_count);
data_struct.x_bins_centers = zeros(1, bins_count);
data_struct.x_bins_widths = zeros(1, bins_count);
data_struct.x_bins_borders = zeros(2, bins_count);
data_struct.x_fine_mesh = zeros(1, fine_mesh_steps_count);
data_struct.x_bins_number = bins_count;
data_struct.x_grad_mesh = zeros(1, bins_count - 1);


% % Bin level
% s1 = 4;
% for l_ind = 1:lambda_count
%     data_struct.dx_fwd_mean_in_bins_saved{l_ind} = zeros(1, bins_count(l_ind));
%     data_struct.dx_bck_mean_in_bins_saved{l_ind} = zeros(1, bins_count(l_ind));
% %     data_struct.MAP_bck_fD_Hanggi{l_ind} = zeros(3, bins_count(l_ind));
%     data_struct.MAP_D{l_ind} = zeros(s1, bins_count(l_ind));
%     data_struct.MAP_fwd_fD_divine{l_ind} = zeros(s1, bins_count(l_ind));
%     data_struct.MAP_fD_Ito{l_ind} = zeros(s1, bins_count(l_ind));
%     data_struct.MAP_fwd_fD_marginalized{l_ind} = zeros(s1, bins_count(l_ind));
%     data_struct.MAP_fwd_fD_Stratonovich{l_ind} = zeros(s1, bins_count(l_ind));
%     data_struct.MAP_fD_mean_sigma2_Ito{l_ind} = zeros(2, bins_count(l_ind));
%     data_struct.MAP_fD_mean_sigma2_Stratonovich{l_ind} = zeros(2, bins_count(l_ind));
%     data_struct.MAP_fD_mean_sigma2_marginalized{l_ind} = zeros(2, bins_count(l_ind));
%     data_struct.MAP_fD_mean_sigma2_divine{l_ind} = zeros(2, bins_count(l_ind));
%     % data_struct.fD_theor{l_ind} = zeros(1, bins_count(l_ind));
%     data_struct.n_j{l_ind} = zeros(1, bins_count(l_ind));
%     data_struct.V_fwd_j{l_ind} = zeros(1, bins_count(l_ind));
%     data_struct.V_bck_j{l_ind} = zeros(1, bins_count(l_ind));
% end;




