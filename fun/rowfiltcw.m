function [m0, m0filt] = rowfiltcw(Mcw, M, alpha, est)
%--------------------------------------------------------------------------
% rowfiltcw - Row filter based on cellwise outliers
%--------------------------------------------------------------------------
%
% Inputs:
%   Mcw: Matrix indicating cellwise outliers (rows represent observations, columns represent variables)
%   M: Matrix indicating missing values (rows represent observations, columns represent variables)
%   alpha: False positive rate (default: 0.01)
%   est: Estimation method for cutoff calculation ('rob' or 'ls', default: 'rob')
%
% Outputs:
%   m0: Logical vector indicating rows filtered as outliers
%   m0filt: Filter information including the test used and the cutoff value
%
% Example:
%   Mcw = [false false; true false; false true]; % Matrix indicating cellwise outliers
%   M = [false true; false false; true false]; % Matrix indicating missing values
%   alpha = 0.01; % False positive rate
%   est = 'rob'; % Estimation method
%   [m0, m0filt] = rowfiltcw(Mcw, M, alpha, est); % Apply row filter

arguments
    Mcw {double,logical}
    M {double,logical}
    alpha double = 0.01
    est {char, string} = 'rob'
end

% If the total amount of cellwise outliers is above the false positive rate assumed (alpha)
if (sum(Mcw(:)) / length(Mcw(:))) > alpha
    k = size(Mcw, 2);
    [m0, m0filt.test, m0filt.ccwpre] = pcwcutoff(Mcw, k, alpha, est);
else
    [~, m0filt.test, m0filt.ccwpre] = pcwcutoff(Mcw, size(Mcw, 2), alpha, est);
    m0 = false(size(Mcw, 1), 1);
end

end
