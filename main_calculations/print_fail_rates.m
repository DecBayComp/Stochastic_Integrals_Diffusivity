


function print_fail_rates(data_struct)
%% == Print out the mean and maximal fail rates for force inference ==
% The mean was calculated over one best explored period and are weighted
% with the bin size

%% Constants
load_constants;


%% Print out
fprintf('\nPrinting out average and max. fail rates for force inference:\n');
% Average fail rates
out_cell = num2cell(data_struct.UR_fD_bin_mean);
out_cell = [lambda_types_names', out_cell];
out_cell = [['> Average fD FR <', conventions_names]; out_cell];
out_cell

% Max fail rates
out_cell = num2cell(data_struct.UR_fD_bin_max);
out_cell = [lambda_types_names', out_cell];
out_cell = [['> Max fD FR <', conventions_names]; out_cell];
out_cell






1;









