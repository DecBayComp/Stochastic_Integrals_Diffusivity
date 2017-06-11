


function print_fail_rates(data_struct)
%% == Print out the mean and maximal fail rates for force inference ==
% The mean was calculated over one best explored period and are weighted
% with the bin size

%% Constants
load_constants;
round_precision = 3;
% round_precision_tex = 0;
tab = char(9);


%% Initialize
output_full_path = strcat(output_figures_folder, fail_rates_filename);


%% Print out
% fprintf('\nPrinting out average and max. fail rates for force inference:\n');
% Average fail rates
out_cell = num2cell(round(data_struct.UR_fD_bin_mean, round_precision));
out_cell = [lambda_types_names', out_cell];
out_cell_avg = [['> Average fD FR <', conventions_names]; out_cell];
out_cell_avg;


% Max fail rates
out_cell = num2cell(round(data_struct.UR_fD_bin_max, round_precision));
out_cell = [lambda_types_names', out_cell];
out_cell_max = [['> Max fD FR <', conventions_names]; out_cell];
out_cell_max;


% Write out to the file
writetable(cell2table([out_cell_avg; cell(2, conventions_count + 1); out_cell_max]), output_full_path,...
    'WriteVariableNames', false, 'Delimiter', '\t');


%% Format the D_UR output
str_output_D_UR = cell(2, 0);
% str_output_D_UR{1} = lambda_types_names;
tmp_str_names = [];
tmp_str_vals = [];
for lambda_type = 1:lambda_types_count
    tmp_str_names = [tmp_str_names, lambda_types_names{lambda_type}, tab];
    tmp_str_vals = [tmp_str_vals,num2str(data_struct.UR_D_bin_mean (lambda_type) * 100, '%.0f'),...
        '(', num2str(data_struct.UR_D_bin_max(lambda_type) * 100, '%.0f'), ')', tab];
end;
% Save dropping the last tab
str_output_D_UR{1} = tmp_str_names(1:end-1);
str_output_D_UR{2} = tmp_str_vals(1:end-1);
str_output_D_UR;


%% Format the fD_UR output to include in the latex file
str_output_fD_UR = cell(conventions_count, 1);
counter = 1;
for convention = [enum_conv_divine, enum_conv_marginalized, enum_conv_Stratonovich, enum_conv_Hanggi, enum_conv_Ito]
    name = conventions_tex_names{convention};
    tmp_str = name;
    for lambda_type = 1:lambda_types_count
        tmp_str = [tmp_str, tab, '&', tab, num2str(data_struct.UR_fD_bin_mean (lambda_type, convention) * 100, '%.0f'),...
            '(', num2str(data_struct.UR_fD_bin_max(lambda_type, convention) * 100, '%.0f'), ')'];
    end;
    str_output_fD_UR{counter} = tmp_str;
    counter = counter + 1;
end;
str_output_fD_UR;


%% Append fD_UR output to the file
% Open file and prepare for output
file = fopen(output_full_path, 'a');
fprintf(file, '\n\nOutput fD fail rate in TEX table format:\n\n\n');
% Output
for convention = 1:conventions_count
    fprintf(file, '%s', str_output_fD_UR{convention});
    fprintf(file, '\n\\\\\n');
end;

% Append D_UR
fprintf(file, '\n\nOutput D fail rate:\n\n\n');
% Output
for i=1:length(str_output_D_UR)
    fprintf(file, '%s', str_output_D_UR{i});
    fprintf(file, '\n\n');
end;


% Close file
fclose(file); 

1;









