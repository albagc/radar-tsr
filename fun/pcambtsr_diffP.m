function [X,m,S,It,difvec,Xrec,pcamodel,SPE,T2,Mrw,Mcw]=pcambtsr_diffP(X,A,xscale,maxiter,tol,alpha,acompmode)
%--------------------------------------------------------------------------
% pcambtsr_diffP - Principal Component Analysis with Missing Data Imputation 
% using TSR
%--------------------------------------------------------------------------
%
% This function performs PCA on a dataset with missing values using the Two-Stage Regression (TSR) approach
% for imputing missing values. It iteratively estimates missing values based on the observed values and updates
% the PCA model until convergence. This version allows specifying the number of principal components to retain (A)
% and the imputation component mode (acompmode).
%
% Inputs:
%   X: Data matrix with missing values. Rows represent observations and columns represent variables.
%   A: Number of principal components to retain in the PCA model.
%   xscale (optional): Logical flag indicating whether to scale the data before performing PCA. Default is false.
%   maxiter (optional): Maximum number of iterations for the TSR algorithm. Default is 5000.
%   tol (optional): Convergence tolerance for the TSR algorithm. Default is 5e-3.
%   alpha (optional): Significance level for computing control limits of SPE and T2 statistics. Default is 0.01.
%   acompmode (optional): Imputation component mode. Default is 'auto'.
%
% Outputs:
%   X: Imputed data matrix with missing values replaced.
%   m: Estimated mean vector of X.
%   S: Estimated covariance matrix of X.
%   It: Number of iterations performed by the TSR algorithm.
%   difvec: Vector containing the differences between successive iterations to track the convergence of the TSR algorithm.
%   Xrec: PCA reconstruction of X using A principal components.
%   pcamodel: Structure containing the PCA model with the following fields:
%       P: Loadings matrix of the PCA model.
%       m: Estimated mean vector of X.
%       s: Estimated standard deviation vector of X.
%       lambda: Eigenvalues of the principal components.
%       ncomp: Number of principal components retained.
%       prepro: Preprocessing method used (either 'cent' or 'autosc').
%       cSPE: Control limit for SPE statistic.
%       cT2: Control limit for T2 statistic.
%   SPE: Squared prediction error statistics for each observation.
%   T2: Hotelling's T2 statistics for each observation.
%   Mrw: Logical vector indicating whether an observation is flagged as an outlier based on SPE or T2 control limits.
%   Mcw: Logical matrix indicating whether a variable is flagged as an outlier based on the contribution to SPE.
%
% Example:
%   X = load('data.mat');
%   A = 3; % Number of principal components to retain
%   xscale = true; % Perform scaling of variables
%   maxiter = 1000; % Maximum number of iterations
%   tol = 1e-4; % Convergence tolerance
%   alpha = 0.05; % Significance level for outlier detection
%   acompmode = 'auto'; % Automatic mode selection for number of principal components
%   [X, m, S, It, difvec, Xrec, pcamodel, SPE, T2, Mrw, Mcw] = pcambtsr_diffP(X, A, xscale, maxiter, tol, alpha, acompmode);
%
% References:
%   - Folch-Fortuny, A., Arteaga, F., Ferrer, A., & Forina, M. (2015). Missing Data
%     Imputation Toolbox for MATLAB.

arguments
    X double
    A (1,1)double
    xscale logical = false
    maxiter (1,1)double = 5000
    tol (1,1)double = 5e-3
    alpha double= 0.01
    acompmode = "auto"
end

Xini = X;
[n,p]=size(X);
mis=isnan(X);

if xscale
    sX = std(X,'omitnan');
    X = X./(repmat(sX,n,1));
else 
    sX = ones(1,p);
end

[r, c]=find(isnan(X));
X(mis)=0;
meanc=sum(X)./(n-sum(mis));
for k=1:length(r)
    X(r(k),c(k))=meanc(c(k));
end
diff=100;
difvec = nan(maxiter,1);
It=0;
% Before the loop
mX=mean(X);
Xc=X-mX;
S=cov(Xc);
if A<1
    A = predncomp(Xc, 0.8, 1, acompmode, 10);
end
if n>p, [~, ~, V]=svd(Xc,0); else, [V, ~, ~]=svd(Xc',0); end
V=V(:,1:A);
while It<maxiter && diff>tol
    
    Xmis=X(mis);
    Pold = V;
    for i=1:n           % for each row
        if sum(mis(i,:))>0     % if there are missing values
            L=V(~mis(i,:),1:min(A,sum(~mis(i,:)))); % L is the key matrix
            S11=S(~mis(i,:),~mis(i,:));
            S21=S(mis(i,:),~mis(i,:));
            z1=Xc(i,~mis(i,:))';
            z2=S21*L*pinv(L'*S11*L)*L'*z1;
            Xc(i,mis(i,:))=z2';
        end
    end
    X = Xc + mX;
    mX = mean(X);
    Xc = X-mX;
    S = cov(Xc);
    if n>p, [~, ~, V]=svd(Xc,0); else, [V, ~, ~]=svd(Xc',0); end
    V = V(:,1:A);
    Pnew = V;
    % As in MacroPCA:
%     diff = acos(sqrt(min(eig(Pnew'*Pold*Pold'*Pnew))));
    % Compute the scalar product between all the PCs
    diff = acos(min(sqrt(diag(Pnew'*Pold).^2)));
    It=It+1;
    difvec(It) = diff;
    
end
X = (Xc + mX).*(repmat(sX,n,1));

if xscale
    pcaprep = "autosc";
else
    pcaprep = "cent";
end
[T,E,SPE,T2,pcamodel,Xrec,T2cont] = pcamb(X, A, alpha, pcaprep, "ls");
m = pcamodel.m;
S = pcamodel.S;
pcamodel.cSPE = pcamodel.limspe*2.5;
pcamodel.cT2 = pcamodel.limt2*2.5;
Mrw = or(SPE > pcamodel.cSPE, T2 > pcamodel.cT2);
Mcw = abs(E./repmat(std(E),size(E,1),1)) > norminv(1-(alpha/2));

end
