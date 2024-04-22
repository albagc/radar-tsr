function [limht2ph1, limht2ph2] = uclht2(T2, A, alpha, mode)
%--------------------------------------------------------------------------
% uclht2 - Upper Confidence Limit for Hotelling's T^2 Test
%--------------------------------------------------------------------------
%
% Inputs:
%   T2: Hotelling's T^2 statistic
%   A: Number of variables (columns) in the data matrix
%   alpha: Significance level (default: 0.05)
%   mode: Calculation mode for the confidence limits ('theor' or 'pctiles', default: 'theor')
%
% Outputs:
%   limht2ph1: Upper confidence limit for Hotelling's T^2 test (theoretical mode) or
%              estimated percentile limit (pctiles mode)
%   limht2ph2: Upper confidence limit for Hotelling's T^2 test (theoretical mode) or
%              estimated percentile limit (pctiles mode)
%
% Example:
%   T2 = 10; % Hotelling's T^2 statistic
%   A = 5; % Number of variables
%   alpha = 0.05; % Significance level
%   mode = 'theor'; % Calculation mode
%   [limht2ph1, limht2ph2] = uclht2(T2, A, alpha, mode); % Calculate upper confidence limits

arguments
    T2 double
    A (1,1) double
    alpha (1,1) double = 0.05
    mode {char, string} = 'theor'
end

switch mode
    case 'theor'
        n = numel(T2);
        
        B_alpha = betainv(1 - alpha, A / 2, (n - A - 1) / 2);
        limht2ph1 = (n - 1) ^ 2 / n * B_alpha;
        
        F_alpha = finv(1 - alpha, A, n - A);
        limht2ph2 = (n ^ 2 - 1) * A / (n * (n - A)) * F_alpha;
        
    case 'pctiles'
        n = length(T2);
        perclim = 1 - alpha;
        numb = max(1, floor(perclim * n));
        rest = n - numb;
        LIM = nan(1000, 1);
        for r = 1:1000
            rand_n = randperm(n);
            sortht2 = sort(T2(rand_n(1:numb)), 'descend');
            LIM(r) = sortht2(rest + 1) + 0.001;
        end
        limht2ph1 = median(LIM);
        limht2ph2 = limht2ph1;
end
end
