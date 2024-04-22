function [Ahat, Asv, Acumsv] = predncomp(X, cumss, svthres, mode, Amax, printtxt)
%--------------------------------------------------------------------------
% predncomp - Predict the number of principal components based on cumulative explained variance
%--------------------------------------------------------------------------
%
% Inputs:
%   X: Data matrix with observations in rows and variables in columns.
%   cumss: Desired cumulative explained variance (optional, default is 0.9).
%   svthres: Singular value threshold (optional, default is 1).
%   mode: Mode of prediction, either 'auto' or 'manual' (optional, default is 'auto').
%   Amax: Maximum number of principal components to consider (optional, default is the number of variables in X).
%   printtxt: Logical value indicating whether to print additional information (optional, default is false).
%
% Outputs:
%   Ahat: Predicted number of principal components.
%   Asv: Number of principal components based on singular value threshold.
%   Acumsv: Number of principal components to reach the desired cumulative explained variance.
%
% Example:
%   X = load('data.mat');
%   cumss = 0.9; % Desired cumulative explained variance
%   svthres = 1; % Singular value threshold
%   mode = 'auto'; % Mode of prediction
%   Amax = size(X, 2); % Maximum number of principal components
%   printtxt = true; % Print additional information
%   [Ahat, Asv, Acumsv] = predncomp(X, cumss, svthres, mode, Amax, printtxt);
%
% Note: The X input should be a data matrix with observations in rows and variables in columns.
arguments
    X double
    cumss (1,1) double = 0.9
    svthres (1,1) double = 1
    mode char = 'auto'
    Amax (1,1) double = size(X,2)
    printtxt logical = false
end

[n, k] = size(X);
vX = sum(var(X));

if n > k
    ss1 = (svd(X, 0).^2) / (n - 1);
else
    ss1 = (svd(X', 0).^2) / (n - 1);
end

Amax = min([Amax, size(X, 2), length(ss1)]);
Asv = nan;
Acumsv = sum(round((cumsum((ss1)) / vX), 2) <= cumss) + 1;

if printtxt
    disp("Suggested number of PCs:")
    disp(strcat(" - Minimum PCs to reach cumulative variance >", string(cumss*100), " % = ", string(Acumsv)))
end

switch mode
    case 'auto'
        Avalues = [Amax, Acumsv];
        Avalues(Avalues <= 0) = [];
        Ahat = min(Avalues);
    case 'manual'
        disp("Suggested number of PCs:")
        disp(strcat(" - Exp. Variance (%) ", strjoin(string(ss1(1:min([Amax, length(ss1)])) / vX * 100), ", ")))
        disp(strcat(" - Minimum PCs to reach cumulative variance >", string(cumss*100), " % = ", string(Acumsv)))

        afig = figure;
        
        subplot(211)
        plot(ss1(1:Amax), 'o-')
        title('Scree plot')
        grid on
        xlabel('PCs')
        
        subplot(212)
        bar(ss1(1:Amax) / vX * 100)
        grid on
        xlabel('PCs')
        hold on
        plot(cumsum(ss1(1:Amax) / vX * 100), 'ro-')
        title('Variance Explained')
        ylabel('Variance (%)')
        ylim([0, 105])
        hold on
        yticks(0:20:100)
        yline(100, 'k--')

        Ahat = input("Select the number of PCs: ");
        close(afig)
end

end
