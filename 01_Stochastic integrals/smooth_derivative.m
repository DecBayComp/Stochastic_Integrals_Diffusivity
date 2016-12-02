

function polynomial_derivative = smooth_derivative(polynomial)


% Extracting polynomial coefficients
[breaks, coefs, l, k, ~] = unmkpp(polynomial);
% Differentiating
dif_matrix = ones(l, 1) * ((k-1):-1:0);
dif_coefs = coefs .* dif_matrix;
dif_coefs = dif_coefs(:, 1:k-1);
% Reconstructing the derivative polynomial
polynomial_derivative = mkpp(breaks, dif_coefs);




