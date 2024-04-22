function [rwfilt, test, PCI] = pcwcutoff(cwest, dobs, alpha, est)
%--------------------------------------------------------------------------
% pcwcutoff - Perform the cutoff test for the estimated percentage of
% cellwise outliers
%--------------------------------------------------------------------------
%
% Inputs:
%   cwest: estimated cellwise matrix.
%   dobs: Number of observations.
%   alpha: Significance level for the cutoff test (optional, default is 0.05).
%   est: Estimation method, either 'ls' (least squares) or 'rob' (robust)
%        (optional, default is 'ls').
%
% Outputs:
%   rwfilt: Logical vector indicating which observations pass the cutoff test.
%   test: Test used for the cutoff ('norm' for normal approximation, 'bin' for binomial approximation).
%   PCI: Percentage Confidence Interval used to set the cut-off value
%
% Example:
%   cwest = load('pcw.mat');
%   dobs = size(cwest, 2);
%   alpha = 0.05;
%   est = 'ls';
%   [rwfilt, test, PCI] = pcwcutoff(cwest, dobs, alpha, est);
%
arguments
    cwest {double,logical}
    dobs double = size(cwest, 2)
    alpha (1,1) double = 0.05
    est {char,string} = 'ls'
end

n = size(cwest, 1);
% Percentage of cw per row:
pcw = sum(cwest, 2) ./ dobs;

switch est
    case 'ls'
        phat = mean(pcw);
    case 'rob'
        phat = robloc(pcw);
end

% Initialize logical m_rw^(0)
rwfilt = false(n, 1);

min_cpre_cw = max(alpha, 2 / dobs);

% If all observations meet the criteria for the normal approximation:
if any(or((pcw * dobs) > 10, (dobs * (1 - pcw)) > 10))
    zci = norminv(1 - alpha);
    se = sqrt(phat * (1 - phat) ./ dobs);
    cpre_rw = phat + zci * se;
    test = 'norm';
else
    xcpre_rw = binoinv(1 - alpha, dobs, mean(pcw));
    cpre_rw = xcpre_rw / dobs;
    test = 'bin';
end

if cpre_rw < min_cpre_cw
    cpre_rw = min_cpre_cw;
end

[spcw, ipcw] = sort(pcw, 'ascend');

ipcw(or(spcw < cpre_rw, spcw == 0)) = [];
nfa = round(alpha * n);

if ~isempty(ipcw)
    if length(ipcw) > nfa
        rwfilt(ipcw) = true;    
    end
end

PCI = cpre_rw;

end
