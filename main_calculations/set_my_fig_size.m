
%% Set my figure size. If the second argument is provided, it sets the [width, height]
% of the image in px

function set_my_fig_size(fig_hand, dimensions)

% Constant
a = 800;


%% Selecting the behavior based on whether the second argument if provided
switch nargin
    case 1
        fig_size = [1.4*a, a];
    case 2
        fig_size = dimensions;
end;

set(fig_hand, 'Units', 'points');
current_position = get(fig_hand, 'Position');

set(fig_hand, 'Position', ([current_position(1:2), fig_size]));













