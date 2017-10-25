%% Fill bins background to represent the bin size


function color_bins(bin_borders, y_lim_vec, bin_color)

%% Initialize
bins_number = size(bin_borders, 2);


%% Plot
bl_color = true;
for bin = 1:bins_number
	% Calculate rectangle position
	rect = [bin_borders(1, bin), y_lim_vec(1),  bin_borders(2, bin) - bin_borders(1, bin), y_lim_vec(2) - y_lim_vec(1)];
	
	% Determine the color
	if bl_color
		cur_color = bin_color;
	else
		cur_color = 'none';
	end;
	bl_color = ~bl_color;
	
	% Color the bin
 	h_rec = rectangle('Position', rect, 'LineStyle', 'none', 'FaceColor',  cur_color);
	
	% Send the rectangle behind all
	uistack(h_rec, 'bottom');
	
% 	% Set transparency
%  	set(h_rec, 'FaceAlpha', 0.5);
end;

% Put axes on top
set(gca, 'Layer', 'top');


