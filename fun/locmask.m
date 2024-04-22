function [mvec] = locmask(X, M, est)
%--------------------------------------------------------------------------
% locmask - Calculate mean not using mask elements set to 1
%--------------------------------------------------------------------------
%   [mvec] = locmask(X, M, est) calculates the mean of the available data
%   in X, excluding the elements where M is set to 1 (or true).
%
%   Input:
%   - X: Data matrix.
%   - M: Mask matrix specifying which elements to exclude from the mean
%        calculation. It can be a logical matrix or a matrix of the same
%        size as X with 1s (or true) indicating the elements to exclude.
%        Default is isnan(X), which excludes NaN values.
%   - est: Estimation method. 'ls' for least squares estimation (default),
%          'rob' for robust estimation using the robloc function.
%
%   Output:
%   - mvec: Mean vector calculated from the available data.
%
%   Example:
%   X = [1, NaN, 3; 4, 5, 6; NaN, 8, 9];
%   M = isnan(X);
%   est = 'ls';
%   mvec = locmask(X, M, est);
%
%--------------------------------------------------------------------------

arguments
   X double
   M {double, logical} = isnan(X)
   est {string, char} = 'ls'
end

Z = X;
[n, p] = size(X);

if strcmp(est, 'ls')
    Z(M) = 0;
    mvec = sum(Z) / (n - sum(M));
    
elseif strcmp(est, 'rob')
    mvec = zeros(1, size(Z, 2));
    
    for j = 1:p
        z1 = Z(M(:, j) == 0, j);
        mvec(j) = robloc(z1, 3, 1e-12);
    end
end

end
