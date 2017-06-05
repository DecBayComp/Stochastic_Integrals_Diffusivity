



l_ind = 2;
bin = 6;


fD_abs_max = 200;
fD_steps_number = 2^12;
step = 2*fD_abs_max/(fD_steps_number - 1);
fD_mesh = -fD_abs_max:step:fD_abs_max;
% fD_steps_number = length(fD_mesh);

data = zeros(1, fD_steps_number);
for ind = 1:fD_steps_number    
    data(ind) = bin_fD_simple_Hanggi_posterior_func (l_ind, bin, fD_mesh(ind));
    if ~mod(ind, 100)
        fprintf('Progress: %i/%i\n', ind, fD_steps_number);
    end
end;


trapz(fD_mesh, data)

% Plot
figure(3);
plot(fD_mesh, data);





