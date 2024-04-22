function [X_imp, X_rec] = tsrme(X, m, S, P, ncomp, prepro, s)
%--------------------------------------------------------------------------
% tsrme - Time Series Regression Missing Values Estimation (TSR-ME)
%--------------------------------------------------------------------------
%
% Inputs:
%   X: Data matrix with missing values
%   m: Mean vector of the training data
%   S: Covariance matrix of the training data (default: [])
%   P: Loadings matrix of the training data (default: [])
%   ncomp: Number of principal components to consider (default: size(P,1))
%   prepro: Data preprocessing method ('autosc', 'cent', 'none', default: 'autosc')
%   s: Standard deviation vector of the training data (default: sqrt(diag(S)))
%
% Outputs:
%   X_imp: Imputed data matrix with missing values replaced
%   X_rec: Reconstructed data matrix with missing values retained
%
% Example:
%   X = [1, NaN, 3; 4, 5, 6; NaN, 8, 9]; % Data matrix with missing values
%   m = [2, 4, 7]; % Mean vector of the training data
%   S = [1, 0.5, 0.3; 0.5, 1, 0.2; 0.3, 0.2, 1]; % Covariance matrix of the training data
%   P = [0.4, 0.5, 0.6; 0.3, 0.2, 0.1]; % Loadings matrix of the training data
%   ncomp = 2; % Number of principal components
%   prepro = 'autosc'; % Data preprocessing method
%   s = [1, 0.5, 0.8]; % Standard deviation vector of the training data
%   [X_imp, X_rec] = tsrme(X, m, S, P, ncomp, prepro, s); % Perform TSR-ME

arguments
   X double
   m 
   S double = []
   P double = []
   ncomp (1,1) double = size(P,1)
   prepro {string, char} = 'autosc'
   s double = sqrt(diag(S))
end

% Check if input is a PCA model struct
if isstruct(m)
    pcamodel = m;
    m = pcamodel.m;
    s = pcamodel.s;
    S = pcamodel.S;
    P = pcamodel.P;
    ncomp = pcamodel.ncomp;
    prepro = pcamodel.prepro;
end

% PCA model elements
n = size(X, 1);
if isfield(pcamodel, 'spre')
    Xini = X;
    X = X ./ (ones(n, 1) * pcamodel.spre);
end

if strcmp(prepro, 'cent')
    X_aux = X - m;
elseif strcmp(prepro, 'autosc')
    X_aux = (X - m) ./ repmat(s, n, 1);
elseif strcmp(prepro, 'none')
    X_aux = X;
end

% TSR-ME
M = isnan(X);
Znew = X_aux;
P = P(:, 1:ncomp);
for i = 1:size(X_aux, 1)
    z = X_aux(i, :);
    zstar = z(M(i, :) == 0);
    Saststar = S(M(i, :) == 1, M(i, :) == 0);
    Pstar = P(M(i, :) == 0, :);
    Star2 = S(M(i, :) == 0, M(i, :) == 0);
    zast = Saststar * Pstar * pinv(Pstar' * Star2 * Pstar) * Pstar' * zstar';
    Znew(i, M(i, :) == 1) = zast;
end

% Come back to the original space
if strcmp(prepro, 'cent')
    X_rec = xpostpro((Znew * P) * P', m);
    X_imp = xpostpro(Znew, m);
elseif strcmp(prepro, 'autosc')
    X_rec = xpostpro((Znew * P) * P', m, s);
    X_imp = xpostpro(Znew, m, s);
elseif strcmp(prepro, 'none')
    X_rec = (Znew * P) * P';
    X_imp = Znew;
end

if isfield(pcamodel, 'spre')
    X_rec = X_rec .* (ones(n, 1) * pcamodel.spre);
    X_imp = X_imp .* (ones(n, 1) * pcamodel.spre);
end

end
