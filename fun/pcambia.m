function [X, m, S, It, diff, Xrec] = pcambia(X, A, maxiter, tol)
%--------------------------------------------------------------------------
% pcambia - Principal Component Analysis (PCA) Missing Data Imputation
%--------------------------------------------------------------------------
%   [X, m, S, It, diff, Xrec] = pcambia(X, A, maxiter, tol) performs PCA-based
%   missing data imputation on the input data matrix X. It replaces the missing
%   values with estimated values based on PCA reconstruction using A components.
%
%   Inputs:
%   - X: Data matrix with NaNs for missing data.
%   - A: Number of principal components to use for reconstruction.
%   - maxiter: Maximum number of iterations for imputation convergence.
%   - tol: Tolerance value for convergence (default is 1e-6).
%
%   Outputs:
%   - X: Data matrix with imputed values.
%   - m: Estimated mean vector of X.
%   - S: Estimated covariance matrix of X.
%   - It: Number of iterations until convergence.
%   - diff: Convergence criterion (mean squared difference between imputed values).
%   - Xrec: PCA reconstruction of X using A components.
%
%   Example:
%   X = load('data.mat');
%   A = 3; % Number of principal components for reconstruction
%   maxiter = 100; % Maximum number of iterations
%   tol = 1e-6; % Tolerance for convergence
%   [X, m, S, It, diff, Xrec] = pcambia(X, A, maxiter, tol);
%
%--------------------------------------------------------------------------

[n, p] = size(X);
mis = isnan(X);
[r, c] = find(mis);
X(mis) = 0;
meanc = sum(X) ./ (n - sum(mis));

for k = 1:length(r)
    X(r(k), c(k)) = meanc(c(k));
end

diff = 100;
It = 0;

while and(It < maxiter, diff > tol)
    It = It + 1;
    Xmis = X(mis);
    mX = mean(X);
    Xc = X - ones(n, 1) * mX;

    if n < p
        [P, D, U] = svd(Xc', 0);
    else
        [U, D, P] = svd(Xc, 0);
    end

    T = U * D;
    T = T(:, 1:A);
    P = P(:, 1:A);
    Xrec = T * P';
    Xrec = Xrec + ones(n, 1) * mX;
    X(mis) = Xrec(mis);
    d = (X(mis) - Xmis).^2;
    diff = mean(d);
end

S = cov(X);
m = mean(X);

