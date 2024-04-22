function [X, pro] = xpostpro(Z, m, s, xscale)
%--------------------------------------------------------------------------
% xpostpro - Postprocessing function to transform standardized data back to the original scale
%--------------------------------------------------------------------------
%
% Inputs:
%   Z: Standardized data matrix
%   m: Mean vector of the original data
%   s: Standard deviation vector of the original data (default: [])
%   xscale: Flag indicating whether scaling was applied (default: ~isempty(s))
%
% Outputs:
%   X: Transformed data matrix in the original scale
%   pro: Postprocessing method used ('cent' for centering, 'autosc' for centering and scaling)
%
% Example:
%   Z = [0.5, 1.0, -0.5; -1.0, 0.0, 1.0]; % Standardized data matrix
%   m = [10, 20, 30]; % Mean vector of the original data
%   s = [2, 2, 2]; % Standard deviation vector of the original data
%   xscale = true; % Flag indicating scaling was applied
%   [X, pro] = xpostpro(Z, m, s, xscale); % Transform data back to the original scale

arguments
    Z double
    m double
    s double = []
    xscale logical = ~isempty(s)
end

if and(~isempty(s), xscale)
    X = (Z .* (ones(size(Z, 1), 1) * s)) + m;
    pro = "autosc";
else
    X = Z + m;
    pro = "cent";
end

end
