


function print_fail_rates(data_struct)
%% == Print out the mean and maximal fail rates for force inference ==
% The mean was calculated over one best explored period and are weighted
% with the bin size

%% Constants
load_constants;
round_precision = 3;


%% Initialize
output_full_path = strcat(output_data_folder, fail_rates_filename);


%% Print out
fprintf('\nPrinting out average and max. fail rates for force inference:\n');
% Average fail rates
out_cell = num2cell(round(data_struct.UR_fD_bin_mean, round_precision));
out_cell = [lambda_types_names', out_cell];
out_cell_avg = [['> Average fD FR <', conventions_names]; out_cell];
out_cell_avg


% Max fail rates
out_cell = num2cell(round(data_struct.UR_fD_bin_max, round_precision));
out_cell = [lambda_types_names', out_cell];
out_cell_max = [['> Max fD FR <', conventions_names]; out_cell];
out_cell_max


% Write out to the file
writetable(cell2table([out_cell_avg; cell(2, conventions_count + 1); out_cell_max]), output_full_path,...
    'WriteVariableNames', false, 'Delimiter', '\t');



1;









