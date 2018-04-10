% Test function


%% Main function of the test suite
function tests = test_select_bins_adaptive_mesh
    tests = functiontests(localfunctions);
end


%% Create a simple data set to test different aspects of the same return results
function setupOnce(test_case)
    % Create a small dataset with known properties (equal spacing)
    L = 1;  % half system size
    N = 1000;   % jumps number
    dx_size = L/(N-1)/2;
    
    x_data = (-1:2/(N-1):1) * L;
    dx_data = ones(N, 1) * dx_size;
    
    % Perform full binning
    max_points_per_bin = -1;    % no limit
    full_binning_data = struct;
    
    [full_binning_data.x_bins_borders, ...
        full_binning_data.x_bins_centers, ...
        full_binning_data.bins_number, ...
        full_binning_data.x_bins_widths,...
        full_binning_data.point_count_in_bins, ...
        full_binning_data.variance_in_bins,...
        full_binning_data.points_binned] = select_bins_adaptive_mesh(x_data, dx_data, max_points_per_bin);
    
    % Impose hard limit of points
    max_points_per_bin = 10;
    limited_binning_data = struct;
    [limited_binning_data.x_bins_borders, ...
        limited_binning_data.x_bins_centers, ...
        limited_binning_data.bins_number, ...
        limited_binning_data.x_bins_widths,...
        limited_binning_data.point_count_in_bins, ...
        limited_binning_data.variance_in_bins,...
        limited_binning_data.points_binned] = select_bins_adaptive_mesh(x_data, dx_data, max_points_per_bin);
    
    limited_binning_data.max_points_per_bin = max_points_per_bin;
    
    % Transfer to test cases
    test_case.TestData.x_data = x_data;
    test_case.TestData.dx_data = dx_data;
    test_case.TestData.full_binning_data = full_binning_data;
    test_case.TestData.limited_binning_data = limited_binning_data;
    
    
end



%% Tests
function test_bin_number(test_case)
    
    exp_bins_number = 100;  % since it is the initial bin number. It should not be increased for this data set
    
    % Full binning
    act_bins_number = test_case.TestData.full_binning_data.bins_number;
    verifyEqual(test_case, act_bins_number, exp_bins_number, "Test bins number for full binning");
    
    % Limited binning
    act_bins_number = test_case.TestData.limited_binning_data.bins_number;
    verifyEqual(test_case, act_bins_number, exp_bins_number, "Test bins number for binning with a hard limit");
end



function test_number_of_points_binned(test_case)
    
    N = length(test_case.TestData.x_data);

    % Full binning
    point_count_in_bins = test_case.TestData.full_binning_data.point_count_in_bins;
    verifyEqual(test_case, sum(point_count_in_bins), N, "Test number of points per bin for full binning");
    
    % Limited binning
    point_count_in_bins = test_case.TestData.limited_binning_data.point_count_in_bins;
    max_points_per_bin = test_case.TestData.limited_binning_data.max_points_per_bin;
    
    tests = point_count_in_bins <= max_points_per_bin;
    verifyTrue(test_case, all(tests));
    
end



function test_small_mean_jump(test_case)
    load_constants;
    
    full_binning_data = test_case.TestData.full_binning_data;
    limited_binning_data = test_case.TestData.limited_binning_data;

    struct_array = {full_binning_data, limited_binning_data};
    for i = 1:length(struct_array)
        cur_struct = struct_array{i};
        
        % Get binned points
        points_binned = cur_struct.points_binned;
        bins_number = cur_struct.bins_number;
        x_bins_widths = cur_struct.x_bins_widths;

        % Calculate mean jump in bins
        dx_mean = zeros(bins_number, 1);
        for bin = 1:bins_number
            dx_mean(bin) = mean(points_binned{bin}(2, :));
        end

        % Compare mean jumps to bin widths
        jump_to_bin_ratios = dx_mean ./ x_bins_widths;

        % Compare this ratio to the minimum ratio defined in constants
        tests = jump_to_bin_ratios < 1/min_bin_to_jump_ratio;
        verifyTrue(test_case, all(tests), "Test that the mean absolute jump is smaller than a certain fraction of bin width");
    end
end


function test_same_points_number(test_case)
    % For symmetric input and full binning, all bins (except maybe the first one) should contain the same number of points

    % Get binned points
    point_count_in_bins = test_case.TestData.full_binning_data.point_count_in_bins;
    
    min_count = min(point_count_in_bins(2:end));
    max_count = max(point_count_in_bins(2:end));
    
    verifyTrue(test_case, min_count == max_count, "Test that for symmetric input all bins contain the same number of points");
end


function test_multiple_limits(test_case)
    x_data = test_case.TestData.x_data;
    dx_data = test_case.TestData.dx_data;
    
    max_points_per_bin = [-1, 10, 100];
    limits_count = length(max_points_per_bin);
    
    % Calculate binning
    [x_bins_borders, x_bins_centers, bins_number, x_bins_widths,...
    point_count_in_bins, variance_in_bins, points_binned] = select_bins_adaptive_mesh(x_data, dx_data, max_points_per_bin);
    
    % Check return values sizes
    assertEqual(test_case, size(bins_number), [1, 1], "Test output size");
    verifyEqual(test_case, size(x_bins_borders), [bins_number, 2], "Test output size");
    verifyEqual(test_case, size(x_bins_centers), [bins_number, 1], "Test output size");
    verifyEqual(test_case, size(x_bins_widths), [bins_number, 1], "Test output size");
    verifyEqual(test_case, size(point_count_in_bins), [limits_count, bins_number], "Test output size");
    verifyEqual(test_case, size(variance_in_bins), [limits_count, bins_number], "Test output size");
    verifyEqual(test_case, size(points_binned), [limits_count, bins_number], "Test output size");
    
    % Check that number of points per bin is correct
    for lim_ind = 2:limits_count
        points_count = point_count_in_bins(lim_ind, :);
        tests = points_count <= max_points_per_bin(lim_ind);
        verifyTrue(test_case, all(tests), "Test that points in bins are hard limited by the input number of points");
    end


end













