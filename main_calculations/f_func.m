

function [f_func_value, U_func_value] = f_func (f_case_number, x, L)

%% Constants
f_shift_7 = 10.0;
f_weak = 40.0;	% fN
f_strong = -120.0;	% fN


%% Select function
switch f_case_number
    case 1 % Flat f0 = 0
        f0 = 0.0;
        f_func_local = @(x) 0 .* x + f0;
        U_func = @(x) -f0 * x;
    case 2 % Flat f0 = 0
        f0 = 10.0;
        f_func_local = @(x) 0 .* x + f0;
        U_func = @(x) -f0 * x;
    case 3 % Flat f0
        f0 = -10.0;
        f_func_local = @(x) 0 .* x + f0;
        U_func = @(x) -f0 * x;
    case 4 % Linearly growing force 
        f_slope = 10.0;
        f_func_local = @(x) f_slope .* (2*x./L);
        U_func = @(x) -f_slope/L .* x.^2;
    case 5 % Linearly decreasing force
        f_slope = -5;
        f_shift = 5;
        f_func_local = @(x) f_shift +  f_slope .* (2*x./L);
        U_func = @(x) - f_shift .* x -f_slope/L .* x.^2;
    case 6 % Linear
        f_shift = 1.5;
        f_slope = 1.0;
        f_func_local = @(x) f_shift +  f_slope .* (x./L);
        U_func = @(x) - f_shift .* x - f_slope/L .* x.^2 / 2;
    case 7 % Flat
        f_shift = f_shift_7;
        f_func_local = @(x) f_shift +  0.0 .* x	;
        U_func = @(x) - f_shift .* x;
	case 8 % Two forces
        f_func_local = @(x) f_weak * (x<=0) +  f_strong * (x>0);
        U_func = @(x) -f_weak * x/L .* (x<=0) -  f_strong * x/L .* (x>0);
end;


%% Return the result
f_func_value = f_func_local(x);
U_func_value = U_func(x);




