function [Xrec, P, mX, lambdaT, vexp] = pcals(X, A, est)
%--------------------------------------------------------------------------
% pcals - Principal Component Analysis (PCA) with different estimation methods
%--------------------------------------------------------------------------
%   [Xrec, P, mX, lambdaT, vexp] = pcals(X, A, est) performs Principal Component
%   Analysis (PCA) on the input data matrix X using different estimation methods.
%   It returns the reconstructed data, principal component loadings, sample
%   means, eigenvalues of the principal components, and the explained variance.
%
%   Inputs:
%   - X: Data matrix with observations in rows and variables in columns.
%   - A: Number of principal components to retain, or the desired cumulative
%        explained variance (A must be between 0 and 1).
%   - est: Estimation method, either 'ls' (least squares) or 'rob' (robust)
%          (optional, default is 'ls').
%
%   Outputs:
%   - Xrec: Reconstructed data matrix using the retained principal components.
%   - P: Principal component loadings.
%   - mX: Sample means of the variables.
%   - lambdaT: Eigenvalues of the retained principal components.
%   - vexp: Explained variance by each principal component.
%
%   Example:
%   X = load('data.mat');
%   A = 0.9; % Retain principal components that explain 90% variance
%   est = 'rob'; % Use robust estimation method
%   [Xrec, P, mX, lambdaT, vexp] = pcals(X, A, est);
%
%--------------------------------------------------------------------------
arguments
    X
    A
    est = "ls"
end
switch est
    case "ls"
        mX = mean(X);
    case "rob"
        mX = robloc(X);
end

[n, p] = size(X);
Xc = X - ones(n, 1) * mX;

if n < p
    [P, D, U] = svd(Xc', 0);
else
    [U, D, P] = svd(Xc, 0);
end

vexp = diag(D).^2 / (n - 1) / sum(var(Xc));

if A < 1
    A = find(cumsum(vexp) >= A, 1, 'first');
end

T = U * D;
T = T(:, 1:A);
P = P(:, 1:A);
Xrec = T * P';
Xrec = Xrec + ones(n, 1) * mX;

switch est
    case "ls"
        lambdaT = var(T);
    case "rob"
        lambdaT = robscale(T).^2;
end

end
