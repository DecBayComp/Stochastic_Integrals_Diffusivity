
%% This function bins points to plot density distributions

x_bins_number = 50;
bin_width = L / x_bins_number;
bins_borders = x_min:bin_width:x_max;
bins_centers = (bins_borders(1:end-1) + bins_borders(2:end))/2;

%% Binning
binned_density = zeros(lambda_count, x_bins_number);
str_legend = cell(1, lambda_count);
for l_ind = 1:lambda_count
    % Binning
    binned_density(l_ind, :) = histcounts(x_lambda(l_ind,:), bins_borders);
    % Normalizing
    binned_density(l_ind, :) = binned_density(l_ind, :) / mean(binned_density(l_ind, :));
    % Preparing legend
    str_legend{l_ind} = sprintf('\\lambda = %.2f (%s)', lambda_array(l_ind), lambda_names_array{l_ind});
end;

%% Preparing to plot the density histogram
figure(2);
clf;
plot(bins_centers/L, binned_density, 'LineWidth', 2);
xlabel('x/L');
ylabel('Normalized density of points');

%% Adding legend
legend(str_legend, 'Location', 'best');






