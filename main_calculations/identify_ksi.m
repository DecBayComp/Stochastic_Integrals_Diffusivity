%% Identify and enumerate simulated ksi values


function [ksi_array, ksi_count, trials_ksi_type, trial_first_ksi_type_index] = identify_ksi(trials_ksi)

trials = length(trials_ksi);

ksi_array = unique(trials_ksi);
ksi_count = length(ksi_array);
trials_ksi_type = zeros(trials, 1);
trial_first_ksi_type_index = zeros(ksi_count, 1);
for ksi_ind = 1:ksi_count
    trials_ksi_type(trials_ksi == ksi_array(ksi_ind)) = ksi_ind;

    % Find first occurrence
    trial_first_ksi_type_index(ksi_ind) = find(trials_ksi_type == ksi_ind, 1);
end
    
end