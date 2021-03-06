

function plot_article_all(stat_struct, trials_data)


%% Globals
load_constants;


%% Constants
bl_save_figures = true;


%% Initialize
fig_count = 0;



%% Plot simulated diffusivity and force
fig_count = fig_count + 1; 
plot_simulated_diffusivity_and_force(fig_count, bl_save_figures);

% % % 
% % % 
% % % 
%% Plot point density and bin locations
fig_count = fig_count + 1; 
plot_article_point_density(stat_struct, fig_count, bl_save_figures);
% % % 
% % % 
% % % 

%% Plot diffusivity posterior in one bin, diffusivity profile and the fail rate
fig_count = fig_count + 1; 
plot_article_b(stat_struct, trials_data, fig_count, bl_save_figures);

% % % 
% % % 
% % % 
% % % % % % %% Compare MAP b distribution with the b posterior from one trial
% % % % % % fig_count = fig_count + 1; 
% % % % % % plot_article_b_MAP_vs_posterior(data_struct, trials_data, fig_count, bl_save_figures);
% % % 
% % % 
% % % 
% % % % % % %% === Plot drift profile in one bin for each inference type ===
% % % % % % fig_count = fig_count + 1; 
% % % % % % plot_article_a_profile_in_bin(data_struct, trials_data, fig_count, bl_save_figures);
% % % 
% % % 
% % % 
% % % %% === Plot mean fD profile and fail rate for each inference type ===
% % % fig_count = fig_count + 1;
% % % plot_article_mean_a_profile_and_fail_rate(data_struct, trials_data, data_struct.bl_force, fig_count, bl_save_figures);
% % % 
% % % 
% % % 
% % % %% === Plot global Bayes factor ===
% % % fig_count = fig_count + 1; 
% % % plot_article_global_bayes_factor(data_struct, data_struct.bl_force, fig_count, bl_save_figures);
% % % 
% % % 
% % % 

%% === Plot local Bayes factor ===
fig_count = fig_count + 1; 
plot_article_local_bayes_factor(stat_struct, fig_count, bl_save_figures);

% % % 
% % % 
% % % 
% % % % % % %% === Plot fD bias as a function fo the gradient ===
% % % % % % fig_count = fig_count + 1; 
% % % % % % plot_article_a_bias(data_struct, fig_count, bl_save_figures);
% % % 
% % % 
% % % %% === Print out different important numbers ===
% % % print_numbers(data_struct);
% % % 
% % % 
% % % 
% % % %% === Print Global Bayes factor ===
% % % disp(data_struct.mean_log_K_G);

















