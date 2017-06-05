


function print_fail_rates(data_struct)

%% Constants
load_constants;
w = 10.0;


%% Initialize
x_left = (1/2 + 2*1)/w;
x_right = (1/2 + 2*2)/w;

%% Calculate
% The means are weighted with the size of bins
AFR_Ito_mean = zeros(lambda_count, 1);
% AFR_Ito_max = zeros(lambda_count, 1);
AFR_Stratonovich_mean = zeros(lambda_count, 1);
AFR_Hanggi_mean = zeros(lambda_count, 1);
% AFR_Stratonovich_max = zeros(lambda_count, 1);
AFR_marginalized_mean = zeros(lambda_count, 1);
% AFR_marginalized_max = zeros(lambda_count, 1);
AFR_divine_mean = zeros(lambda_count, 1);
% AFR_divine_max = zeros(lambda_count, 1);
UR_Ito_mean = zeros(lambda_count, 1);
UR_Ito_max = zeros(lambda_count, 1);
UR_Stratonovich_mean = zeros(lambda_count, 1);
UR_Stratonovich_max = zeros(lambda_count, 1);
UR_Hanggi_mean = zeros(lambda_count, 1);
UR_Hanggi_max = zeros(lambda_count, 1);
UR_marginalized_mean = zeros(lambda_count, 1);
UR_marginalized_max = zeros(lambda_count, 1);
UR_oracle_mean = zeros(lambda_count, 1);
UR_oracle_max = zeros(lambda_count, 1);
for l_ind = 1:lambda_count
    %% One-trial fail rate
    % Filter indices from one best period
    indices = data_struct.x_bins_centers_saved{l_ind} >= x_left & data_struct.x_bins_centers_saved{l_ind} <= x_right;
    bin_widths = data_struct.x_bins_centers_saved{l_ind}(indices)';
    bin_widths = bin_widths / sum(bin_widths);
    % Ito
    UR_Ito_mean(l_ind) = sum(data_struct.UR_Ito{l_ind}(indices) .* bin_widths);
    UR_Ito_max(l_ind) = max(data_struct.UR_Ito{l_ind}(indices));
    % Stratonovich
    UR_Stratonovich_mean(l_ind) = sum(data_struct.UR_fwd_Stratonovich{l_ind}(indices) .* bin_widths);
    UR_Stratonovich_max(l_ind) = max(data_struct.UR_fwd_Stratonovich{l_ind}(indices));
    % Hanggi
    UR_Hanggi_mean(l_ind) = sum(data_struct.UR_Hanggi{l_ind}(indices) .* bin_widths);
    UR_Hanggi_max(l_ind) = max(data_struct.UR_Hanggi{l_ind}(indices));
    % Marginalized
    UR_marginalized_mean(l_ind) = sum(data_struct.UR_fwd_marginalized{l_ind}(indices) .* bin_widths);
    UR_marginalized_max(l_ind) = max(data_struct.UR_fwd_marginalized{l_ind}(indices));
    % Oracle
    UR_oracle_mean(l_ind) = sum(data_struct.UR_fwd_divine{l_ind}(indices) .* bin_widths);
    UR_oracle_max(l_ind) = max(data_struct.UR_fwd_divine{l_ind}(indices));
    %% Asymptotic fail rate
    % Ito
    difference = (data_struct.MAP_fD_Ito_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)');
    half_error = (data_struct.MAP_fD_Ito_mean{l_ind}(2, indices) + data_struct.MAP_fD_Ito_mean{l_ind}(3, indices))/2; 
    difference_rel = abs(difference ./ half_error);
    fail = double(difference_rel > 1);
    AFR_Ito_mean(l_ind) = sum(fail .* bin_widths);
%     AFR_Ito_max(l_ind) = max(fail);
    % Stratonovich
    difference = (data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)');
    half_error = (data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(2, indices) + data_struct.MAP_fwd_fD_Stratonovich_mean{l_ind}(3, indices))/2; 
    difference_rel = abs(difference ./ half_error);
    fail = double(difference_rel > 1);
    AFR_Stratonovich_mean(l_ind) = sum(fail .* bin_widths);
    % Hanggi
    difference = (data_struct.MAP_fD_Hanggi_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)');
    half_error = (data_struct.MAP_fD_Hanggi_mean{l_ind}(2, indices) + data_struct.MAP_fD_Hanggi_mean{l_ind}(3, indices))/2; 
    difference_rel = abs(difference ./ half_error);
    fail = double(difference_rel > 1);
    AFR_Hanggi_mean(l_ind) = sum(fail .* bin_widths);
%     AFR_Stratonovich_max(l_ind) = max(fail);
    % Marginalized
    difference = (data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)');
    half_error = (data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(2, indices) + data_struct.MAP_fwd_fD_marginalized_mean{l_ind}(3, indices))/2; 
    difference_rel = abs(difference ./ half_error);
    fail = double(difference_rel > 1);
    AFR_marginalized_mean(l_ind) = sum(fail .* bin_widths);
%     AFR_marginalized_max(l_ind) = max(fail);
    % Divine
    difference = (data_struct.MAP_fwd_fD_divine_mean{l_ind}(1, indices) - data_struct.fD_theor_data{l_ind}(indices)');
    half_error = (data_struct.MAP_fwd_fD_divine_mean{l_ind}(2, indices) + data_struct.MAP_fwd_fD_divine_mean{l_ind}(3, indices))/2; 
    difference_rel = abs(difference ./ half_error);
    fail = double(difference_rel > 1);
    AFR_divine_mean(l_ind) = sum(fail .* bin_widths);
%     AFR_divine_max(l_ind) = max(fail);
end;


%% % Print out
fprintf('\n');
% Custom
% Ito for Ito
fprintf('Average Ito fail rate within [%.2f, %.2f] for Ito simulation only: %.1f%%\n', x_left, x_right, UR_Ito_mean(1) * 100);
% Marginalized for Stratonovich
fprintf('Average marginalized fail rate within [%.2f, %.2f] for Stratonovich simulation only: %.1f%%\n', x_left, x_right, UR_marginalized_mean(2) * 100);
% Oracle for Stratonovich
fprintf('Average oracle fail rate within [%.2f, %.2f] for Stratonovich simulation only: %.1f%%\n', x_left, x_right, UR_oracle_mean(2) * 100);

%%% Fail rate

%% Max
fprintf('\n');
% Ito
fprintf('Max Ito fail rate within [%.2f, %.2f]: %.1f%%\n', x_left, x_right, max(UR_Ito_max) * 100);
% Stratonovich
fprintf('Max Stratonovich fail rate: %.1f%%\n', max(UR_Stratonovich_max) * 100);
% Stratonovich
fprintf('Max Hanggi fail rate: %.1f%%\n', max(UR_Hanggi_max) * 100);
% Marginalized
fprintf('Max marginalized fail rate: %.1f%%\n', max(UR_marginalized_max) * 100);
% Oracle
fprintf('Max oracle fail rate: %.1f%%\n', max(UR_oracle_max) * 100);
% Average
fprintf('\n');
% Ito
fprintf('Average Ito fail rate: %.1f%%\n', mean(UR_Ito_mean) * 100);
% Stratonovich
fprintf('Average Stratonovich fail rate: %.1f%%\n', mean(UR_Stratonovich_mean) * 100);
% Stratonovich
fprintf('Average Hanggi fail rate: %.1f%%\n', mean(UR_Hanggi_mean) * 100);
% Marginalized
fprintf('Average marginalized fail rate: %.1f%%\n', mean(UR_marginalized_mean) * 100);
% Oracle
fprintf('Average oracle fail rate: %.1f%%\n', mean(UR_oracle_mean) * 100);


%%% Asymptotic fail rate
% % % fprintf('\n');
% % % %% Max
% % % % Ito
% % % fprintf('Asymptotic max Ito fail rate within [%.2f, %.2f]: %.1f%%\n', x_left, x_right, max(AFR_Ito_max) * 100);
% % % % Stratonovich
% % % fprintf('Asymptotic max Stratonovich fail rate within [%.2f, %.2f]: %.1f%%\n', x_left, x_right, max(AFR_Stratonovich_max) * 100);
%% Average
fprintf('\n');
% Ito
fprintf('Asymptotic average Ito fail rate: %.1f%%\n', mean(AFR_Ito_mean) * 100);
% Hanggi
fprintf('Asymptotic average Hanggi fail rate: %.1f%%\n', mean(AFR_Hanggi_mean) * 100);
% Stratonovich
fprintf('Asymptotic average Stratonovich fail rate: %.1f%%\n', mean(AFR_Stratonovich_mean) * 100);
% Marginalized
fprintf('Asymptotic average marginalized fail rate: %.1f%%\n', mean(AFR_marginalized_mean) * 100);
% Divine
fprintf('Asymptotic average divine fail rate: %.1f%%\n', mean(AFR_divine_mean) * 100);



1;









