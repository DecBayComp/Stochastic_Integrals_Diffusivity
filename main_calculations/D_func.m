 

function [D_func_value, D_prime_func_value, D_scnd_der_func_value, D_antider_value] = D_func (D_case_number, x, L)


%% Iniitialize
D_scnd_der_func_value = zeros('like', x);
D_antider_func = @(x) NaN .* x;
D_0 = 1/25/pi;  % 1e-2;	% um^2/s


%% Select function
switch D_case_number
    case 1  % Parabolic diffusivity
        D_min = 1;
        D_max = 2;
        D_func_local = @(x) D_min + (D_max - D_min) * (2*x/L).^2;
        D_prime_func = @(x) (D_max - D_min) * 8 * x / L^2;
            %% Test replacement
%         D_func_local = @(x) D_min + 0.*x;
%         D_prime_func = @(x) 0.*x;
        
    case 2 % Linear D
        D_min = 1;
        D_max = 2;
        D_func_local = @(x) D_min + (x/L + 1/2) * (D_max - D_min);
        D_prime_func = @(x) (D_max - D_min)/L;    
    
    case 3 % Discontinuity in D
        D_min = 1;
        D_max = 2;
        D_func_local = @(x) D_min .* (x<0) + D_max .* (x>=0);
        D_prime_func = @(x) 0 .* x;
    
    case 4 % Exponentially growing D
        D_min = 1;
        D_max = 2;
        alpha = 7*log(2);
        
        A = (D_max - D_min) / (exp(alpha) - 1);
        B = D_min - A;
        D_func_local = @(x) A * exp(alpha * (x + 1/2)) + B;
        D_prime_func = @(x) alpha * A * exp(alpha * (x + 1/2));
    case 5 % Linear
        D_shift = 1.5;
        D_slope = 1.0;
        D_func_local = @(x) D_shift + D_slope * x/L;
        D_prime_func = @(x) D_slope * 1.0/L;
	case 6 % Oscillating
        w = 10.0;	% um^(-1)
        D_func_local = @(x) D_0/2.0 * (2 + sin(pi * w * x));
        D_prime_func = @(x) D_0/2.0 * pi * w * cos(pi * w * x);
        D_scnd_der_func = @(x) - D_0/2.0 * (pi * w).^2 * sin(pi * w * x);
        D_antider_func = @(x) D_0 * (x - cos(pi * w * x) / (2 * pi * w));
		
		D_scnd_der_func_value = D_scnd_der_func(x);
		
		
% % %         % Overriding boundary conditions
% % %         str_mode = 'left_wall_only'; 

end;


%% Return the result
D_func_value = D_func_local(x);
D_prime_func_value = D_prime_func(x);
D_antider_value = D_antider_func(x);



