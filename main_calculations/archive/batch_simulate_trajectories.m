

%% Constants
load_constants;
% bl_batch_start = 1;


for D_case_number = 1:max_D_case_number
    for f_case_number = 1:max_f_case_number
%         bl_simulation_succeeded = false;
%         for try_num = 1:SIMULATION_TRIES_PER_CASE
            simulate_trajectories_only;
%             if bl_simulation_succeeded
%                 break;
%             end;
%         end;
    end;
end;

