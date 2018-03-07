


function [num_files_missing, missing_indices_list] = check_input_sequential(input_data_folder)

% Initialize
num_files_missing = 0;
missing_indices_list = [];

fprintf("Checking consistency of filenames in folder '%s'...\n", input_data_folder);


% Find the highest simulated index
cur_dir = dir([input_data_folder, '*.csv']);
last_filename = cur_dir(end).name;
last_index = extractBetween(last_filename, "sim_data_", ".csv");
last_index = str2num(last_index{1});


% Scan filenames
for file_num = 1:last_index
    
    % Construct full path
    filename = sprintf('sim_data_%09i.csv', file_num);
    output_full_path = strcat(input_data_folder, filename);
    
    if ~exist(output_full_path, 'file')
        num_files_missing = num_files_missing + 1;
        missing_indices_list = [missing_indices_list; file_num];
    end
    
end


str_result = "Checking complete. ";
if num_files_missing > 0
    str_result = strcat(str_result, sprintf("%i files missing.", num_files_missing));
else
    str_result = strcat(str_result, "No problems found.");
end

disp(str_result);









end


