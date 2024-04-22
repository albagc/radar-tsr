function [limspe] = uclspe(SPE, alpha, mode, E)
%--------------------------------------------------------------------------
% uclspe - Upper Confidence Limit for SPE (Squared Prediction Error)
%--------------------------------------------------------------------------
%
% Inputs:
%   SPE: Squared Prediction Error values
%   alpha: Significance level (default: 0.05)
%   mode: Calculation mode for the confidence limit ('theor', 'pctiles', or 'nomikmac', default: 'theor')
%   E: Eigenvector matrix (required for 'theor' mode)
%
% Outputs:
%   limspe: Upper confidence limit for SPE
%
% Example:
%   SPE = [0.5, 0.7, 0.9]; % Squared Prediction Error values
%   alpha = 0.05; % Significance level
%   mode = 'theor'; % Calculation mode
%   E = [1, 2, 3; 4, 5, 6; 7, 8, 9]; % Eigenvector matrix
%   limspe = uclspe(SPE, alpha, mode, E); % Calculate upper confidence limit

arguments
    SPE double
    alpha (1,1) double = 0.05
    mode {char, string} = 'theor'
    E double = []
end

switch mode
    case 'theor'
        LambdaE = eig(cov(E));
        
        % Calculate SPE upper control limit
        theta1 = sum(LambdaE);
        theta2 = sum(LambdaE.^2);
        theta3 = sum(LambdaE.^3);
        h0 = 1 - 2 * theta1 * theta3 / (3 * theta2^2);
        z_alpha = norminv(1 - alpha);
        spe1 = z_alpha * sqrt(2 * theta2 * (h0^2)) / theta1;
        spe2 = theta2 * h0 * (h0 - 1) / (theta1^2);
        limspe = theta1 * ((spe1 + 1 + spe2)^(1 / h0));
        
    case 'pctiles'
        n = length(SPE);
        perclim = 1 - alpha;
        numb = max(1, floor(perclim * n));
        rest = n - numb;
        LIM = nan(1000, 1);
        for r = 1:1000
            rand_n = randperm(n);
            sortspe = sort(SPE(rand_n(1:numb)), 'descend');
            LIM(r) = sortspe(rest + 1) + 0.001;
        end
        limspe = median(LIM);
        
    case 'nomikmac'
        b = mean(SPE);
        v = var(SPE);
        dof = 2 * b^2 / v;
        xchi = chi2inv(1 - alpha, dof);
        limspe = v / (2 * b) * xchi;
end
end
