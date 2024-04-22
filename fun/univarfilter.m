function [univout, uvfilt] = univarfilter(X, alpha, est, M, mode, xcent, xscale)
%--------------------------------------------------------------------------
% univarfilter - Univariate outlier filter assuming normal distribution of variables
%--------------------------------------------------------------------------
%
% Inputs:
%   X: Data matrix
%   alpha: Significance level (default: 0.05)
%   est: Estimation method ('rob' for robust estimation, 'ls' for least squares estimation, default: 'rob')
%   M: Mask matrix specifying which elements to exclude from the filter (default: true(size(X)))
%   mode: Filter mode ('auto' or 'manual', default: 'auto')
%   xcent: Flag indicating whether to center the variables (default: true)
%   xscale: Flag indicating whether to scale the variables (default: true)
%
% Outputs:
%   univout: Matrix indicating the univariate outliers
%   uvfilt: Structure containing the filter parameters (m: mean, s: standard deviation, c: cutoff value)
%
% Example:
%   X = [1, 2, 3; 4, 5, 6; 7, 8, 9]; % Data matrix
%   alpha = 0.05; % Significance level
%   est = 'rob'; % Estimation method
%   M = true(size(X)); % Mask matrix
%   mode = 'auto'; % Filter mode
%   xcent = true; % Flag for centering
%   xscale = true; % Flag for scaling
%   [univout, uvfilt] = univarfilter(X, alpha, est, M, mode, xcent, xscale); % Apply univariate filter

arguments
    X double
    alpha (1,1) double = 0.05
    est {char, string} = "rob"
    M logical = true(size(X))
    mode {char, string} = "auto"
    xcent logical = true
    xscale logical = true
end

% Calculate cutoff value for bilateral test
c_uniout = norminv(1 - alpha/2);

% Center and scale the variables if necessary
if xcent
    mX = locmask(X, M, est);
else
    mX = zeros(1, size(X, 2));
end

if xscale
    sX = scalemask(X, M, est);
else
    sX = ones(1, size(X, 2));
end

Xaux = (X - mX) ./ (ones(size(X, 1), 1) * sX);
Xaux(M) = 0;
univout = abs(Xaux) > c_uniout;

% Handle cases with too few non-missing values
if any(sum(~logical(M + univout)) <= 3)
    if strcmp(mode, "manual")
        figure, subplot(131), imagesc(M), title('Missing entries'),
        xlabel('Variables'), ylabel('Observations')
        
        subplot(132), imagesc(univout), title('Big univariate outliers'),
        xlabel('Variables'), ylabel('Observations')
        
        deleting_mask = false(size(M));
        deleting_mask(:, sum(~logical(M + univout)) <= 3) = true;
        
        subplot(133), imagesc(deleting_mask), title('Deleted variables'),
        xlabel('Variables'), ylabel('Observations'), colormap(flipud(flag))
        
        disp("The following variables will be deleted if univariate outliers are set to NaN: " + ...
            strjoin(string(find(sum(~logical(M + univout)) <= 3)), ","))
        
        cb2na = input("Set UV outliers to NaN anyway [1(y)/0(n)]? ");
        if cb2na ~= 1
            univout(:, sum(~logical(M + univout)) <= 3) = false;
        end
    else
        univout(:, sum(~logical(M + univout)) <= 3) = false;
    end
end

% Check if the proportion of outliers is less than or equal to alpha
if round(sum(univout(:)) / (size(X, 1) * size(X, 2)), 2) <= alpha
    univout = false(size(X));
end

% Store filter parameters in the output structure
uvfilt.m = mX;
uvfilt.s = sX;
uvfilt.c = c_uniout;
end
