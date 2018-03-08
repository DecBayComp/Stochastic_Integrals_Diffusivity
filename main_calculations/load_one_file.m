% Read one file with a specified number from the input data folder


function [D_case, ksi, x, dx, bl_file_corrupt] = load_one_file(input_data_folder, file_num)

% Constants
load_constants;
KSI_PRECISION = 1e-2;

filename = sprintf('sim_data_%09i.csv', file_num);
% fprintf("Loading trajectory from '%s'. Progress: %i/%i\n", filename, file_num, input_files_count);
full_path = strcat(input_data_folder, filename);
bl_file_corrupt = 0;

try
    input_data = dlmread(full_path, CSV_DELIMITER);

    D_case = input_data(1,1);
    ksi = round(input_data(1,2), -log10(KSI_PRECISION));
    x = input_data(2:N+1, 1);
    dx = input_data(2:N+1,2);
    
catch ME
    fprintf("Error reading from file '%s'. The file may be corrupt. \n", filename);
    disp(strcat("Error: ", getReport(ME)));
    bl_file_corrupt = 1;
%     corrupted_files_count = corrupted_files_count +1;
%     corrupted_files_list = [corrupted_files_list; file_num];
end
        
end