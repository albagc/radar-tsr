function [Xrec,T,E,SPE,T2,T2matrix] = pcame(X, pcamodel)
%--------------------------------------------------------------------------
% pcame - Project data onto a PCA model
%--------------------------------------------------------------------------
%
% Inputs:
%   X: Data matrix with observations in rows and variables in columns.
%   pcamodel: PCA model structure obtained from pcamb or pcals function.
%
% Outputs:
%   Xrec: Reconstructed data matrix using the PCA model.
%   T: Scores (latent variables) of the projected data.
%   E: Residuals (prediction errors) of the projected data.
%   SPE: Squared prediction error for each observation.
%   T2: Hotelling's T-squared statistic for each observation.
%   T2matrix: Matrix of T2 statistics for each observation and principal component.
%
% Example:
%   X = load('data.mat');
%   pcamodel = load('pcamodel.mat');
%   [Xrec, T, E, SPE, T2, T2matrix] = pcame(X, pcamodel);
%
% Note: The pcamodel struct should contain the necessary fields obtained from the pcamb or pcals function.

arguments
    X double
    pcamodel struct
end
% Project the data onto the PCA model in pcamodel struct
n = size(X,1);
if isfield(pcamodel,'spre')
    Xini = X;
    X = X./(ones(n,1)*pcamodel.spre);
end
switch pcamodel.prepro
    case 'cent'
        Xaux = X - pcamodel.m;
    case 'autosc'
        Xaux = (X - pcamodel.m) ./ (ones(n,1)*pcamodel.s);
    case 'none'
        Xaux = X;
end
P = pcamodel.P(:,1:pcamodel.ncomp);
T = Xaux * P;

switch pcamodel.prepro
    case 'cent'
        Xrec = T * P' + pcamodel.m;
    case 'autosc'
        Xrec = ((T * P').* (ones(n,1)*pcamodel.s)) + pcamodel.m ;
    case 'none'
        Xrec = T * P';
end
if isfield(pcamodel,'spre')
    Xrec = Xrec.*(ones(n,1)*pcamodel.spre);
end
E = Xaux - T * P';
SPE = sum(E .^ 2, 2);
T2matrix = (T .^ 2) ./ repmat(pcamodel.lambda(1:pcamodel.ncomp), size(T,1),1);
T2 = sum(T2matrix, 2);
end
