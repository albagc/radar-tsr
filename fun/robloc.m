function [robL] = robloc(X, c1, epsilon)
%--------------------------------------------------------------------------
% robloc - Robust local location estimator
%--------------------------------------------------------------------------
%
% Inputs:
%   X: Data matrix (rows represent observations, columns represent variables)
%   c1: Tuning constant for robustness (default: 3)
%   epsilon: Small positive value to avoid division by zero (default: 1e-12)
%
% Output:
%   robL: Robust local location estimator for each variable in X
%
% Example:
%   X = [1 2 3; 4 5 NaN; 6 7 8]; % Data matrix with missing values
%   robL = robloc(X); % Calculate robust local location estimator
%
% Note: The output robL is a row vector containing the robust local location estimator for each variable in X.

arguments
    X double
    c1 (1,1) double = 3
    epsilon (1,1) double = 1e-12
end

d = size(X, 2); % Number of variables
robL = nan(1, d); % Initialize output vector

for j = 1:d
    x = X(:, j); % Extract column j from X
    x = x(~isnan(x)); % Remove NaN values
    
    m1 = median(x); % Median of the non-NaN values
    s1 = c1 * median(abs(x - m1)); % Robust scale estimate
    
    if s1 > epsilon
        x1 = abs(x - m1) / s1;
        w = 1 - x1 .* x1;
        w = ((abs(w) + w) / 2) .^ 2;
        robL(j) = sum(x .* w) / sum(w); % Robust weighted average
    else
        robL(j) = m1; % Use median as the location estimate
    end
end

end
