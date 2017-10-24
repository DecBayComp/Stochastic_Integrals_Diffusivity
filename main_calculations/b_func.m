%% Return b(x) based through querying D(x)


function b = b_func(selected_D_case, x, L)

[d, ~, ~] = D_func(selected_D_case, x, L);
b = sqrt(2 * d);