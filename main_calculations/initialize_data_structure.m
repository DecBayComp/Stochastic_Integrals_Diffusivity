

function data_struct = initialize_data_structure(lambda_count, bins_count, fine_mesh_steps_count)

%% Comments
% - bins_count is a vector of length lambda_count


%%
data_struct = struct();

% Lambda level
% data_struct.elements_in_bins_count = cell(1, lambda_count);
data_struct.dx_Mean = zeros(1, lambda_count);
data_struct.dx_fwd_mean_in_bins_saved = cell(1, lambda_count);
data_struct.dx_bck_mean_in_bins_saved = cell(1, lambda_count);
% data_struct.points_in_bins = cell(1, max_bins_number);
data_struct.MAP_fD_Hanggi = cell(1, lambda_count);
data_struct.MAP_fD_Hanggi_mean = cell(1, lambda_count);
data_struct.MAP_fwd_D = cell(1, lambda_count);
data_struct.MAP_fwd_D_regular = cell(1, lambda_count);
data_struct.MAP_fwd_D_grad_regular = cell(1, lambda_count);
data_struct.MAP_fwd_D_grad_regular_interp = cell(1, lambda_count);
data_struct.MAP_fwd_fD_divine = cell(1, lambda_count);
data_struct.MAP_fwd_fD_divine_mean = cell(1, lambda_count);
data_struct.MAP_D_mean = cell(1, lambda_count);
data_struct.D_theor_data = cell(1, lambda_count);                          % At bins centers
data_struct.D_theor_fine_data = zeros(1, fine_mesh_steps_count);
data_struct.D_grad_theor_fine_data = zeros(1, fine_mesh_steps_count);
data_struct.fD_theor_data = cell(1, lambda_count);                          % At bins centers
data_struct.fD_theor_fine_data = zeros(1, fine_mesh_steps_count);           % On a fine mesh
data_struct.MAP_fD_Ito = cell(1, lambda_count);
data_struct.MAP_fD_Ito_mean = cell(1, lambda_count);
data_struct.MAP_fD_Ito_widths = cell(1, lambda_count);
data_struct.MAP_fwd_fD_marginalized = cell(1, lambda_count);
data_struct.MAP_fwd_fD_marginalized_mean = cell(1, lambda_count);
data_struct.MAP_fwd_fD_Stratonovich = cell(1, lambda_count);
data_struct.MAP_fwd_fD_Stratonovich_mean = cell(1, lambda_count);
data_struct.MAP_fD_mean_sigma2_Ito = cell(1, lambda_count);
data_struct.MAP_fD_mean_sigma2_Stratonovich = cell(1, lambda_count);
data_struct.MAP_fD_mean_sigma2_marginalized = cell(1, lambda_count);
data_struct.MAP_fD_mean_sigma2_divine = cell(1, lambda_count);
% data_struct.fD_theor = cell(1, lambda_count);
data_struct.n_j = cell(1, lambda_count);
data_struct.UR_D = cell(1, lambda_count);
data_struct.UR_fwd_divine = cell(1, lambda_count);
data_struct.UR_Ito = cell(1, lambda_count);
data_struct.UR_fwd_marginalized = cell(1, lambda_count);
data_struct.UR_fwd_Stratonovich = cell(1, lambda_count);
data_struct.UR_Hanggi = cell(1, lambda_count);
data_struct.V = zeros(1, lambda_count);
data_struct.V_fwd_j = cell(1, lambda_count);
data_struct.V_bck_j = cell(1, lambda_count);
data_struct.x_bins_centers_saved = cell(1, lambda_count);
data_struct.x_bins_widths_saved = cell(1, lambda_count);
data_struct.x_fine_mesh = zeros(1, fine_mesh_steps_count);
data_struct.x_bins_number = bins_count;
data_struct.x_grad_mesh = cell(1, lambda_count);


% Bin level
s1 = 4;
for l_ind = 1:lambda_count
    data_struct.dx_fwd_mean_in_bins_saved{l_ind} = zeros(1, bins_count(l_ind));
    data_struct.dx_bck_mean_in_bins_saved{l_ind} = zeros(1, bins_count(l_ind));
%     data_struct.MAP_bck_fD_Hanggi{l_ind} = zeros(3, bins_count(l_ind));
    data_struct.MAP_D{l_ind} = zeros(s1, bins_count(l_ind));
    data_struct.MAP_fwd_fD_divine{l_ind} = zeros(s1, bins_count(l_ind));
    data_struct.MAP_fD_Ito{l_ind} = zeros(s1, bins_count(l_ind));
    data_struct.D_theor_data{l_ind} = zeros(3, bins_count(l_ind));
    data_struct.fD_theor_data{l_ind} = zeros(1, bins_count(l_ind));
    data_struct.MAP_fwd_fD_marginalized{l_ind} = zeros(s1, bins_count(l_ind));
    data_struct.MAP_fwd_fD_Stratonovich{l_ind} = zeros(s1, bins_count(l_ind));
    data_struct.MAP_fD_mean_sigma2_Ito{l_ind} = zeros(2, bins_count(l_ind));
    data_struct.MAP_fD_mean_sigma2_Stratonovich{l_ind} = zeros(2, bins_count(l_ind));
    data_struct.MAP_fD_mean_sigma2_marginalized{l_ind} = zeros(2, bins_count(l_ind));
    data_struct.MAP_fD_mean_sigma2_divine{l_ind} = zeros(2, bins_count(l_ind));
    % data_struct.fD_theor{l_ind} = zeros(1, bins_count(l_ind));
    data_struct.n_j{l_ind} = zeros(1, bins_count(l_ind));
    data_struct.V_fwd_j{l_ind} = zeros(1, bins_count(l_ind));
    data_struct.V_bck_j{l_ind} = zeros(1, bins_count(l_ind));
    data_struct.x_bins_centers_saved{l_ind} = zeros(1, bins_count(l_ind));
    data_struct.x_bins_widths_saved{l_ind} = zeros(1, bins_count(l_ind));
end;




