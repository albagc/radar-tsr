function [X,m,S,It,difvec,Xrec,pcamodel,SPE,T2,Mrw,Mcw]=pcambtsr(X,A,xscale,maxiter,tol,alpha)
%--------------------------------------------------------------------------
% pcambtsr - Principal Component Analysis (PCA) with Missing Data Imputation using TSR
%--------------------------------------------------------------------------
% 
% This function performs PCA on a dataset with missing values using the
% Two-Stage Regression (TSR) approach for imputing missing values. It
% iteratively estimates missing values based on the observed values and
% updates the PCA model until convergence.
%
% Syntax:
%   [X, m, S, It, difvec, Xrec, pcamodel, SPE, T2, Mrw, Mcw] = pcambtsr(X, A, xscale, maxiter, tol, alpha)
%
% Input arguments:
%   - X: Data matrix with missing values. Rows represent observations and
%        columns represent variables.
%   - A: Number of principal components to retain in the PCA model.
%   - xscale (optional): Logical flag indicating whether to scale the data
%        before performing PCA. Default is false.
%   - maxiter (optional): Maximum number of iterations for the TSR algorithm.
%        Default is 5000.
%   - tol (optional): Convergence tolerance for the TSR algorithm. Default is 10e-10.
%   - alpha (optional): Significance level for computing control limits of SPE and T2
%        statistics. Default is 0.01.
%
% Output arguments:
%   - X: Imputed data matrix with missing values replaced.
%   - m: Estimated mean vector of X.
%   - S: Estimated covariance matrix of X.
%   - It: Number of iterations performed by the TSR algorithm.
%   - difvec: Vector containing the differences between successive iterations
%        to track the convergence of the TSR algorithm.
%   - Xrec: PCA reconstruction of X using A principal components.
%   - pcamodel: Structure containing the PCA model with the following fields:
%       - P: Loadings matrix of the PCA model.
%       - m: Estimated mean vector of X.
%       - s: Estimated standard deviation vector of X.
%       - lambda: Eigenvalues of the principal components.
%       - ncomp: Number of retained principal components.
%       - prepro: Preprocessing method applied to the data ('cent' or 'autosc').
%       - cSPE: Control limit for Squared Prediction Error (SPE) statistic.
%       - cT2: Control limit for T^2 statistic.
%   - SPE: Squared Prediction Error (SPE) statistic for each observation.
%   - T2: T^2 statistic for each observation.
%   - Mrw: Logical vector indicating outliers based on the SPE control limit.
%   - Mcw: Logical matrix indicating outliers based on the T^2 control limit.
%
% References:
%   - Folch-Fortuny, A., Arteaga, F., & Ferrer, A. (2015). Missing Data
%     Imputation Toolbox v1.0. Retrieved from https://www.mathworks.com/matlabcentral/fileexchange/49172-missing-data-imputation-toolbox
%   - Savorani, F., & Tomasi, G. (2011). An introduction to
%     N-way toolbox for MATLAB. Chemometrics and Intelligent Laboratory Systems, 118, 24-32.
%
arguments
    X double
    A (1,1)double
    xscale logical = false
    maxiter (1,1)double = 5000
    tol (1,1)double = 10e-10
    alpha double= 0.01
end

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
while It<maxiter && diff>tol
    It=It+1;
    Xmis=X(mis);
    mX=mean(X);
    Xc=X-mX;
    S=cov(Xc);
    if A==0
        A = predncomp(Xc, 0.8,1,'auto',10);
    end
    if n>p
        [~, ~, V]=svd(Xc,0); 
    else 
        [V, ~, ~]=svd(Xc',0); 
    end
    V=V(:,1:A);
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
    d=(X(mis)-Xmis).^2;
    diff=mean(d);
    difvec(It) = diff;
end
X = Xc.*(repmat(sX,n,1)) + mX;
m = mean(X);
if xscale
    sX = std(X);
end
S = cov(Xc);
[~, ~, v]=svd(S,0);
P=v(:,1:A);
T = ((X - m)./repmat(sX,n,1))*P;
Xrec= T*P';
E = ((X - m)./repmat(sX,n,1)) - Xrec;
Xrec=(Xrec.*repmat(sX,n,1)) + m;
if xscale
    pcaprep = "autosc";
else
    pcaprep = "cent";
end

SPE = sum(E.^2,2);
T2 = sum((T.^2)./repmat(var(T),size(T,1),1),2);
clspe = uclspe(SPE, alpha, "nomikmac", E);
clt2 = uclht2(T2, A, alpha, "theor");
pcamodel = struct('P',P,'m',m,'s',std(X),'lambda',var(T),'ncomp',A, ...
    'prepro',pcaprep,...
    'cSPE',2.5*clspe,'cT2',2.5*clt2);
Mrw = or(SPE > 2.5*clspe, T2 > 2.5*clt2);
Mcw = abs(E./repmat(std(E),size(E,1),1)) > norminv(1-(alpha/2));

end
