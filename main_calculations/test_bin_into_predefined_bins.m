


%% Main function of the test suite
function tests = test_bin_into_predefined_bins
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
    x_bins_borders = -1:0.1:1;
    
    % Slightly extend to fit all points
    x_bins_borders(1) = x_bins_borders(1) * 1.01;
    x_bins_borders(end) = x_bins_borders(end) * 1.01;
    
    
    x_bins_borders = [x_bins_borders(1:end-1); x_bins_borders(2:end)]';
    
    % Perform full binning
    n_limit = -1;    % no limit
    full_binning_data = struct;
    
    [full_binning_data.elements_in_bins_count,...
        full_binning_data.points_in_bins,...
        full_binning_data.bl_empty_bins] = bin_into_predefined_bins(x_data, dx_data, x_bins_borders, n_limit);
    
    % Impose hard limit of points
    n_limit = 10;
    limited_binning_data = struct;
    [limited_binning_data.elements_in_bins_count,...
        limited_binning_data.points_in_bins,...
        limited_binning_data.bl_empty_bins] = bin_into_predefined_bins(x_data, dx_data, x_bins_borders, n_limit);
    
    limited_binning_data.n_limit = n_limit;
    
    % Transfer to test cases
    test_case.TestData.x_data = x_data;
    test_case.TestData.dx_data = dx_data;
    test_case.TestData.full_binning_data = full_binning_data;
    test_case.TestData.limited_binning_data = limited_binning_data;
    
    
end


%% Tests
function test_number_of_points_binned(test_case)
    data = test_case.TestData.full_binning_data;
    
    N = length(test_case.TestData.x_data);

    % Full binning
    elements_in_bins_count = data.elements_in_bins_count;
    verifyEqual(test_case, sum(elements_in_bins_count), N, "Test number of points per bin for full binning");
    
    % Limited binning
    elements_in_bins_count = test_case.TestData.limited_binning_data.elements_in_bins_count;
    n_limit = test_case.TestData.limited_binning_data.n_limit;
    
    tests = elements_in_bins_count == n_limit;
    verifyTrue(test_case, all(tests), "Test number of points per bin for hard-limit binning");
    
end


function test_same_points_number(test_case)
    % For symmetric input and full binning, all bins (except maybe the last one) should contain the same number of points

    % Get binned points
    elements_in_bins_count = test_case.TestData.full_binning_data.elements_in_bins_count;
    
    min_count = min(elements_in_bins_count(1:end-1));
    max_count = max(elements_in_bins_count(1:end-1));
    
    verifyTrue(test_case, min_count == max_count, "Test that for symmetric input all bins contain the same number of points");
end


function test_lie_within_bins(test_case)
    data = test_case.TestData.full_binning_data;
    N = length(test_case.TestData.x_data);

%     % Full binning
%     elements_in_bins_count = data.elements_in_bins_count;
%     verifyEqual(test_case, sum(elements_in_bins_count), N, "Test number of points per bin for full binning");
%     
%     % Limited binning
%     elements_in_bins_count = test_case.TestData.limited_binning_data.elements_in_bins_count;
%     n_limit = test_case.TestData.limited_binning_data.n_limit;
%     
%     tests = elements_in_bins_count == n_limit;
%     verifyTrue(test_case, all(tests), "Test number of points per bin for hard-limit binning");
    
end






