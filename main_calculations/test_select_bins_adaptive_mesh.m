% Test function


%% Main function of the test suite
function tests = test_select_bins_adaptive_mesh
    tests = functiontests(localfunctions);
end


%% Create a simple data set to test different aspects of the same return results
function setup(test_case)
    % Create a small dataset with known properties (equal spacing)
    L = 1;  % half system size
    N = 1000;   % jumps number
    dx_size = L/(N-1)/2;
    
    x_data = (-1:2/(N-1):1) * L;
    dx_data = ones(N, 1) * dx_size;
    
    % Perform full binning
    max_points_per_bin = N;
    full_binning_data = struct;
    
    [full_binning_data.x_bins_borders, ...
        full_binning_data.x_bins_centers, ...
        full_binning_data.bins_number, ...
        full_binning_data.x_bins_widths,...
        full_binning_data.point_count_in_bins, ...
        full_binning_data.variance_in_bins,...
        full_binning_data.points_binned] = select_bins_adaptive_mesh(x_data, dx_data, max_points_per_bin);
    
    % Transfer to test cases
    test_case.TestData.x_data = x_data;
    test_case.TestData.dx_data = dx_data;
    test_case.TestData.full_binning_data = full_binning_data;
    
end



%% Tests
function test_full_binning_bin_number(test_case)
    
    exp_bins_number = 100;  % since it is the initial bin number. It should not be increased for this data set
    act_bins_number = test_case.TestData.full_binning_data.bins_number;
    verifyEqual(test_case, act_bins_number, exp_bins_number);
end



function test_full_binning_all_points_binned(test_case)
    
    point_count_in_bins = test_case.TestData.full_binning_data.point_count_in_bins;
    N = length(test_case.TestData.x_data);
    
    verifyEqual(test_case, sum(point_count_in_bins), N);
end



function test_full_binning_small_mean_jump(test_case)
    % Get binned points
    points_binned = test_case.TestData.full_binning_data.points_binned;
    bins_number = test_case.TestData.full_binning_data.bins_number;
    x_bins_widths = test_case.TestData.full_binning_data.x_bins_widths;
    
    % Calculate mean jump in bins
    dx_mean = zeros(bins_number, 1);
    for bin = 1:bins_number
        dx_mean(bin) = mean(points_binned{bin}(2, :));
    end
    
    % Compare mean jumps to bin widths
    jump_to_bin_ratios = dx_mean ./ x_bins_widths;
    
    % Compare this ratio to the minimum ratio defined in constants
    load_constants;
    tests = jump_to_bin_ratios < 1/min_bin_to_jump_ratio;
    verifyTrue(test_case, all(tests));
end


function test_full_binning_same_points_number(test_case)
    % For symmetric input all bins (except maybe the first one) should
    % contain the same number of points

    % Get binned points
    point_count_in_bins = test_case.TestData.full_binning_data.point_count_in_bins;
    
    min_count = min(point_count_in_bins(2:end));
    max_count = max(point_count_in_bins(2:end));
    
    verifyTrue(test_case, min_count == max_count);
end
















