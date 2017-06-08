

function plot_article_all(data_struct, trials_data)


%% Globals
load_constants;


%% Constants
bl_save_figures = true;


%% Initialize
fig_count = 0;
% str_legend = {'Orcl', 'Ito', 'Str', 'Mar'};


%% Plot diffusivity posterior in one bin, diffusivity profile and the fail rate
fig_count = fig_count + 1; 
plot_article_D(data_struct, trials_data, fig_count, bl_save_figures);
% Plot D bias as a function of D gradient ???


%% === Plot fD profile in one bin for each inference type ===
fig_count = fig_count + 1; 
plot_article_fD_profile_in_bin(data_struct, trials_data, fig_count, bl_save_figures);


%% === Plot mean fD profile and fail rate for each inference type ===
fig_count = fig_count + 1; 
plot_article_mean_fD_profile_and_fail_rate(data_struct, fig_count, bl_save_figures);


%% === Plot fD bias as a function fo the gradient ===
plot_article_fD_bias(data_struct, fig_count, bl_save_figures);


%% === Plot fD error cones ===
fig_count = plot_article_fD_error_cones(data_struct, fig_count, bl_save_figures);


%% === Print out fail rates ===
print_fail_rates(data_struct);




















