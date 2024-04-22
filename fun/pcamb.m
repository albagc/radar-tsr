function [T, E, SPE, T2, pcamodel, Xrec, T2cont] = pcamb(X, A, alpha, prepro, est, limitscrit_spe, limitscrit_t2, threshold, acompmode, thresholdmode)
%--------------------------------------------------------------------------
% pcamb - Principal Component Analysis Monitoring (PCAM) with different estimation methods
%--------------------------------------------------------------------------
%   [T, E, SPE, T2, pcamodel, Xrec, T2cont] = pcamb(X, A, alpha, prepro, est, limitscrit_spe, limitscrit_t2, threshold, acompmode, thresholdmode)
%   performs Principal Component Analysis Monitoring (PCAM) on the input data matrix X
%   using different estimation methods. It returns the scores, residuals,
%   Squared Prediction Error (SPE), Hotelling's T^2 statistic, PCA model,
%   reconstructed data, and T^2 contributions for each observation.
%
%   Inputs:
%   - X: Data matrix with observations in rows and variables in columns.
%   - A: Number of principal components to retain or the desired cumulative
%        explained variance (A must be between 0 and the number of variables).
%   - alpha: Significance level for control limits (default is 0.05).
%   - prepro: Preprocessing method, either 'cent' for centering, 'autosc' for
%             autoscaling, or 'none' for no preprocessing (default is 'none').
%   - est: Estimation method, either 'ls' (least squares) or 'rob' (robust)
%          (default is 'ls').
%   - limitscrit_spe: Control limit calculation method for SPE, either 'nomikmac'
%                     (Nomikos and MacGregor's approximation) or 'threshperc'
%                     (threshold percentage of control samples) (default is 'nomikmac').
%   - limitscrit_t2: Control limit calculation method for T^2, either 'theor'
%                    (theoretical limits) or 'spec' (using specified control limits)
%                    (default is 'theor').
%   - threshold: Threshold value used for determining the number of principal components
%                when A is specified as a cumulative explained variance (default is 0.8).
%   - acompmode: Automatic mode for determining the number of principal components,
%                either 'auto' (based on threshold) or 'manual' (A is manually specified)
%                (default is 'auto').
%   - thresholdmode: Threshold mode used for determining the number of principal components
%                    when A is specified as a threshold value, either 'cumsum' (cumulative sum
%                    of explained variance) or 'eachpc' (individual explained variance)
%                    (default is 'cumsum').
%
%   Outputs:
%   - T: Scores matrix.
%   - E: Residuals matrix.
%   - SPE: Squared Prediction Error for each observation.
%   - T2: Hotelling's T^2 statistic for each observation.
%   - pcamodel: PCA model struct containing various model elements.
%   - Xrec: Reconstructed data matrix using the retained principal components.
%   - T2cont: T^2 contributions for each observation.
%
%   Example:
%   X = load('data.mat');
%   A = 3; % Retain 3 principal components
%   alpha = 0.01; % Significance level of control limits
%   prepro = 'cent'; % Centering preprocessing
%   est = 'rob'; % Robust estimation method
%   limitscrit_spe = 'nomikmac'; % Nomikos and MacGregor's approximation for SPE control limits
%   limitscrit_t2 = 'theor'; % Theoretical control limits for T^2
%   threshold = 0.85; % Threshold for determining the number of principal components
%   acompmode = 'auto'; % Automatic mode for determining the number of principal components
%   thresholdmode = 'cumsum'; % Threshold mode for determining the number of principal components
%   [T, E, SPE, T2, pcamodel, Xrec, T2cont] = pcamb(X, A, alpha, prepro, est, limitscrit_spe, limitscrit_t2, threshold, acompmode, thresholdmode);
%
%--------------------------------------------------------------------------
arguments
    X double
    A double {mustBeInRangeSize(A,X,2)}
    alpha double {mustBeInRange(alpha, [0,1])} = 0.05
    prepro string {mustBeMember(prepro, ["cent", "autosc", "none"])} = "none"
    est {char,string} = "ls"
    limitscrit_spe {char,string} = "nomikmac"
    limitscrit_t2 {char,string} = "theor"
    threshold double = 0.8
    acompmode {char,string} = "auto"
    thresholdmode {char,string} = "cumsum"
end
% PCA Model Building
switch est
    case "ls"
        m = mean(X);
        s = std(X);
    case "rob"
        m = robloc(X);
        s = robscale(X);
end

[n, k] = size(X);

% Preprocessing
if strcmp(prepro, 'cent')
    Xaux = X - m;
elseif strcmp(prepro, 'autosc')
    Xaux = (X - m) ./ (ones(n, 1) * s);
elseif strcmp(prepro, 'none')
    Xaux = X;
end

% Singular Value Decomposition
if n < k
    [V, D, ~] = svd(Xaux', 0);
else
    [~, D, V] = svd(Xaux, 0);
end

vexp = diag(D).^2 / (n - 1) / sum(var(Xaux));

% Determine the number of principal components to retain
if and(A < 1, strcmp(acompmode, "auto"))
    if thresholdmode == "cumsum"
        A = min([find(cumsum(vexp) >= threshold, 1, 'first'), 10]);
    elseif thresholdmode == "eachpc"
        A = min([sum(vexp >= threshold), 10]);
    end
elseif strcmp(acompmode, "manual")
    A = predncomp(Xaux, threshold, 1, 'manual', 10);
end

% Compute scores, residuals, and SPE
P = V(:, 1:A);
d = (diag(D).^2) ./ (n - 1);
Tfull = Xaux * P;
T = Xaux * P;
E = Xaux - (T * P');
SPE = sum(E.^2, 2);

% Compute T^2 statistic and T^2 contributions
switch est
    case "ls"
        Lambda = var(T);
    case "rob"
        Lambda = robscale(T).^2;
end

lambda = Lambda(1:A);
T2cont = (T.^2) ./ repmat(lambda, n, 1);
T2 = sum(T.^2 ./ repmat(lambda, n, 1), 2);

% Reconstruct data
if strcmp(prepro, 'cent')
    Xrec = (T * P') + m;
elseif strcmp(prepro, 'autosc')
    Xrec = ((T * P') .* repmat(s, n, 1)) + m;
elseif strcmp(prepro, 'none')
    Xrec = T * P';
end

% Calculate control limits
cl_spe = uclspe(SPE, alpha, limitscrit_spe, E);
[cl_t2, cl_t2_ph2] = uclht2(T2, A, alpha, limitscrit_t2);

% Limits for scores
z = ((n - 1) * (n - 1) / n) * betainv(1 - alpha, 1, (n - 3) / 2);
limits_t = struct();

switch est
    case "ls"
        for i = 1:A
            limits_t.(strcat('pc', string(i))) = [-sqrt(var(Tfull(:, i)) * z), sqrt(var(Tfull(:, i)) * z)];
        end
    case "rob"
        for i = 1:A
            limits_t.(strcat('pc', string(i))) = [-sqrt((robscale(Tfull(:, i)).^2) * z), sqrt((robscale(Tfull(:, i)).^2) * z)];
        end
end

% Output struct with PCA model elements
pcamodel.m = m;
pcamodel.s = s;
pcamodel.P = P;
pcamodel.eigv = d(1:A);
pcamodel.vexp = d ./ sum(var(Xaux));
pcamodel.sstot = sum(var(Xaux));
pcamodel.eigvall = d;
pcamodel.Pfull = V;
pcamodel.lambda = lambda;
pcamodel.limspe = cl_spe;
pcamodel.limt2 = cl_t2;
pcamodel.limt2_ph2 = cl_t2_ph2;
pcamodel.prepro = prepro;
pcamodel.ncomp = A;
pcamodel.alpha = alpha;
pcamodel.n = n;
pcamodel.S = cov(X);
pcamodel.limits_t = limits_t;
end

function mustBeInRange(arg, b)
    if (arg < b(1)) || (arg > b(2))
        error(['Value assigned to Data is not in range ', num2str(b(1)), '...', num2str(b(2))])
    end
end

function mustBeInRangeSize(arg, B, dimB)
    if (arg < 0) || (arg > size(B, dimB))
        error(['Value assigned to Data is not in range ', num2str(2), '...', num2str(size(B, dimB))])
    end
end
