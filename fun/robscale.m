function [robS] = robscale(X, c2, delta, epsilon)
%--------------------------------------------------------------------------
% robscale - Robust scale estimator
%--------------------------------------------------------------------------
%
% Inputs:
%   X: Data matrix (rows represent observations, columns represent variables)
%   c2: Tuning constant for robustness (default: 2.5)
%   delta: Scaling constant for robustness (default: 0.844472)
%   epsilon: Small positive value to avoid division by zero (default: 1e-12)
%
% Output:
%   robS: Robust scale estimator for each variable in X
%
% Example:
%   X = [1 2 3; 4 5 NaN; 6 7 8]; % Data matrix with missing values
%   robS = robscale(X); % Calculate robust scale estimator
%
% Note: The output robS is a row vector containing the robust scale estimator for each variable in X.

arguments
    X double
    c2 (1,1) double = 2.5
    delta (1,1) double = 0.844472
    epsilon (1,1) double = 1e-12
end

X = X - robloc(X); % Center X using the robust local location estimator
d = size(X, 2); % Number of variables
robS = nan(1, d); % Initialize output vector

for j = 1:d
    x = X(:, j); % Extract column j from X
    x = x(~isnan(x)); % Remove NaN values
    n = length(x); % Number of non-NaN values
    
    s0 = median(abs(x)); % Robust scale estimate
    
    if c2 * s0 < epsilon
        robS(j) = s0; % Use median absolute deviation as the scale estimate
    else
        x = x / s0; % Rescale x
        rho = x .^ 2;
        rho(rho > (c2 ^ 2)) = c2 ^ 2; % Truncate squared values
        robS(j) = s0 * sqrt(sum(rho) / (n * delta)); % Robust scaling factor
    end
end

end
