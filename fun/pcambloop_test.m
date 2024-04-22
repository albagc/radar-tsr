function [tsrmbloop, SPE, T2, pcamodel, output_test] = pcambloop_test(X, A, Mmd, Mmdcw, obsid, alpha, est, ...
    tolangle, toldiff, maxit, deleterw, limspecrit, limt2crit, xscale,...
    min_ntpctge, acompmode, mtrue, Xnai)
%--------------------------------------------------------------------------
% pcambloop - Principal Component Analysis (PCA) Missing Data Imputation
%              with Robust Row and Column Weighting and Looping
%--------------------------------------------------------------------------
%   [tsrmbloop, SPE, T2, pcamodel] = pcambloop(X, A, Mmd, Mmdcw, obsid, alpha, est,
%   tolangle, toldiff, maxit, deleterw, limspecrit, limt2crit, xscale, min_ntpctge,
%   acompmode, mtrue, Xnai) performs PCA-based missing data imputation with
%   robust row and column weighting and looping. It replaces the missing
%   values in the input data matrix X with estimated values based on PCA
%   reconstruction using A components. The algorithm iteratively updates the
%   imputed values, applies robust row and column weighting, and checks for
%   convergence based on the change in the P matrix.
%
%   Inputs:
%   - X: Data matrix with NaNs for missing data.
%   - A: Number of principal components to use for reconstruction.
%   - Mmd: Mask matrix indicating missing values in X (default: isnan(X)).
%   - Mmdcw: Mask matrix indicating missing values to be treated by robust
%            column weighting (default: isnan(X)).
%   - obsid: Indices of observations in X (default: 1:size(X,1)).
%   - alpha: Significance level for robust weighting (default: 0.05).
%   - est: Estimation method for missing data imputation ('ls' or 'rob',
%          default: 'ls').
%   - tolangle: Tolerance for convergence based on change in P matrix
%               (default: 5e-3).
%   - toldiff: Tolerance for convergence based on change in imputed values
%              (default: 1e-6).
%   - maxit: Maximum number of iterations (default: 5000).
%   - deleterw: Flag indicating whether to delete rows with outliers (default: false).
%   - limspecrit: Limit criterion for SPE (default: 'nomikmac').
%   - limt2crit: Limit criterion for T2 (default: 'theor').
%   - xscale: Flag indicating whether to scale X (default: false).
%   - min_ntpctge: Minimum percentage of non-outlier observations required
%                  for each principal component (default: 0.5).
%   - acompmode: Mode for determining the number of components (default: 'auto').
%   - mtrue: Mask vector indicating which rows contain true missing values
%            (default: false(size(X,1),1)).
%   - Xnai: Data matrix with NaNs for non-available data (default: X).
%
%   Outputs:
%   - tsrmbloop: Structure containing the results of the missing data
%                imputation with robust weighting and looping.
%   - SPE: Squared prediction error.
%   - T2: Hotelling's T2 statistic.
%   - pcamodel: PCA model structure with P matrix, mean vector, standard
%               deviation vector, lambda (variance) values, number of
%               components, and preprocessing information.
%
%   Example:
%   X = load('data.mat');
%   A = 3; % Number of principal components for reconstruction
%   Mmd = isnan(X); % Mask matrix for missing values
%   Mmdcw = Mmd; % Mask matrix for robust column weighting
%   obsid = 1:size(X,1); % Indices of observations
%   alpha = 0.05; % Significance level for robust weighting
%   est = 'ls'; % Estimation method for missing data imputation
%   tolangle = 5e-3; % Tolerance for convergence based on change in P matrix
%   toldiff = 1e-6; % Tolerance for convergence based on change in imputed values
%   maxit = 5000; % Maximum number of iterations
%   deleterw = false; % Flag to delete rows with outliers
%   limspecrit = 'nomikmac'; % Limit criterion for SPE
%   limt2crit = 'theor'; % Limit criterion for T2
%   xscale = false; % Flag to scale X
%   min_ntpctge = 0.5; % Minimum percentage of non-outlier observations per PC
%   acompmode = 'auto'; % Mode for determining the number of components
%   mtrue = false(size(X,1),1); % Mask vector indicating true missing values
%   Xnai = X; % Data matrix with NaNs for non-available data
%   [tsrmbloop, SPE, T2, pcamodel] = pcambloop(X, A, Mmd, Mmdcw, obsid, alpha, est,
%   tolangle, toldiff, maxit, deleterw, limspecrit, limt2crit, xscale, min_ntpctge,
%   acompmode, mtrue, Xnai);
%
%--------------------------------------------------------------------------

arguments
    X double
    A {double,char}
    Mmd {logical, double} = isnan(X)
    Mmdcw {logical, double} = isnan(X)
    obsid double = 1:size(X,1)
    alpha (1,1)double = 0.05
    est {string, char} = 'ls'
    tolangle (1,1)double = 5e-3
    toldiff (1,1)double = 1e-6
    maxit (1,1)double = 5000
    deleterw (1,1)logical = false
    limspecrit {string, char} = 'nomikmac'
    limt2crit {string, char} = 'theor'
    xscale logical = false
    min_ntpctge (1,1)double = 0.5
    acompmode {string,char} = 'auto'
    mtrue = false(size(X,1),1)
    Xnai = X
end
X0 = X;
Xini = X;
[n0,p0] = size(X0);

% Compute scales
if xscale
    sx0 = scalemask(X,Mmd,est);
    prepro = "autosc";
    % Standardization (mandatory and undone at the very END)
    X = X./repmat(sx0, n0, 1);
else
    sx0 = ones(1,p0);
    prepro = "cent";
end
colmis = find(or(or(sum(~isnan(X))<=3, sum(isinf(X)))>0, (sx0==0)));

if ~isempty(colmis)
    remvar = true;
    X(:,colmis) = [];
    X0(:,colmis) = [];
    sx0(colmis) = [];
    Xnai(:,colmis) = [];
    Mmd(:,colmis) = [];
    Mmdcw(:,colmis) = [];
    if length(colmis) == 1
        disp(strcat("Variable ", strjoin(string(colmis), ", "), ...
            " was missing for all observations and was not included in the model"))
    else
        disp(strcat("Variables: ", strjoin(string(colmis), ", "), ...
         " were missing for all observations and were not included in the model"))
    end
else
    remvar = false;
end

[n,p] = size(X);
M0 = Mmd;
miscw0 = Mmdcw;
mX = locmask(X,Mmd,est);

% Preprocessing with observed mean centering (or scaling)
Meanc = repmat(mX,n,1);
X(Mmdcw) = Meanc(Mmdcw);
Xc = X-Meanc;

pcaprep = "cent";

if n > p, [~, ~, Pold] = svd(Xc,0); else, [Pold,~,~] = svd(Xc',0); end
if or(ischar(A), A == 0)
    A = predncomp(Xc, 0.8,1,acompmode,10);
    min_ntpctge = max([A/size(Xc,1), min_ntpctge]);
end

Pold = Pold(:,1:A);
%
output_test.Pfirst = Pold;
%
loopdiff = 100;
loopanglep = 1;
del_rows = true;
it = 0;
mx = mX; 
hP = nan(1,maxit);
hD = nan(1,maxit);
nobs = nan(1,maxit);
ncells = nan(1,maxit);
Mrwloop = false(size(X,1),1);
obsidloop = 1:size(X,1);
Xnai_it = Xnai./repmat(sx0,n,1);
% Check that the minimum amount of observations specified is feasible
% according to the number of PCs specified
if A > 0
    min_ntpctge = max([A/size(Xc,1), min_ntpctge]);
end

tsr_convergence = false;

while and(or(~tsr_convergence, del_rows), it < maxit)
    nloop = size(X,1);
    nobs(it+1) = nloop/n0*100;
    ncells(it+1) = sum(~(Mmdcw(:)))/(n0*p0)*100;
    mis = Mmdcw;
    Xmis = X(mis); % Missing values (original),
    Xc = X-mx;
    if it==1
        if any(T2 > pcamodel.limt2)
            S = cov(Xc(T2loop<prctile(T2,50),:));
        end
    else
        S = cov(Xc);
    end
    V = Pold(:, 1:A);
    for i = 1:nloop           % for each row
        if sum(mis(i,:))>0     % if there are missing values
            L = V(~mis(i,:),1:min(A,sum(~mis(i,:)))); % L is the key matrix
            S11 = S(~mis(i,:),~mis(i,:));
            S21 = S(mis(i,:),~mis(i,:));
            z1 = Xc(i,~mis(i,:))';
            z2 = S21*L*pinv(L'*S11*L)*L'*z1;
            Xc(i,mis(i,:)) = z2';
        end
    end
    Ximp = Xc + mx;
    if it==1
        if any(T2 > pcamodel.limt2)
            [~,~,~,~,pcamodel] = pcamb(Ximp(T2loop<prctile(T2,50),:), A, alpha, pcaprep, est, limspecrit, limt2crit);
        end
    else
        
        [~,~,~,~,pcamodel] = pcamb(Ximp, A, alpha, pcaprep, est, limspecrit, limt2crit);
        %
        output_test.Ximp_first = Ximp;
        output_test.pca_first = pcamodel;
        %
    end
    [Xrec,T,E,SPE,T2] = pcame(Ximp, pcamodel);
    cwit = univarfilter(Xnai_it - Xrec, alpha, 'rob', false(size(Xnai_it)), acompmode);
    [rw_cwloop, rw_cw_filt] = rowfiltcw(cwit, Mmd, alpha, "rob");
    Mmdcw = Mmd;
    if rw_cw_filt.ccwpre > alpha
        Mmdcw(~rw_cwloop,:) = cwit(~rw_cwloop,:);
    end
    mx = pcamodel.m;
    Pnew = pcamodel.P;
    A = pcamodel.ncomp;
    [rwloop, ~, rwt2] = rw_detect(SPE, T2, pcamodel.limspe, pcamodel.limt2, 1, 2.5, alpha);
    rwloop(sum(cwit,2)==1) = false;
    T2loop = T2(~rwloop);
    del_rows = any(rwloop);
    if and(del_rows, deleterw)
        if (size(Ximp(~rwloop,:),1) >= round(min_ntpctge*n0))
            Mrwloop(obsidloop(rwloop)) = true;
            obsid = obsid(~rwloop); % Delete them from the rows vector
            obsidloop = obsidloop(~rwloop);
            X = Ximp(~rwloop,:);
            Xnai_it = Xnai_it(~rwloop,:);
            Mmd = Mmd(~rwloop,:);
            Mmdcw = Mmdcw(~rwloop,:);
        else
            X = Ximp;
        end
    else
        X = Ximp;
    end
    it = it + 1;

    loopanglep = acos(min(sqrt(diag(Pnew'*Pold).^2)));
    d = mean((Ximp(mis) - Xmis).^2);
    tsr_convergence = loopanglep <= tolangle;
    loopdiff = mean(d);
    
    hP(it) = loopanglep;
    hD(it) = loopdiff;
    Pold = Pnew;
end
%
output_test.pca_end = pcamodel;
%
Xnaimp = tsrme(X0./repmat(sx0,n,1), pcamodel);
[Xrec,T,E,SPE,T2] = pcame(Xnaimp, pcamodel);
Rall = Xnai./repmat(sx0,n,1) - Xnaimp;
cwloop = univarfilter(Rall, alpha, "rob", M0, acompmode);
cwloop(obsidloop,:) = and(Mmdcw, ~Mmd);
cwloop(sum(cwloop,2)/size(cwloop,2) > rw_cw_filt.ccwpre,:) = false;
Xcell = X0./repmat(sx0,n,1);
Xcell(cwloop) = nan;
Xcellimp = tsrme(Xcell, pcamodel);
[Xrec,T,E,SPE,T2] = pcame(Xcellimp, pcamodel);
[~,rwloop] =  rw_detect(SPE, T2, pcamodel.limspe, pcamodel.limt2, 2.5, 2.5, alpha);
pcamodel.m = pcamodel.m.*sx0;
pcamodel.s = sx0;
pcamodel.prepro = prepro;
L = Pnew;
if remvar
    Xremvars = Xini;
    Xremvars(obsidloop,setdiff(1:p0, colmis)) = X.*repmat(sx0,size(X,1),1);
else
    Xremvars = X.*repmat(sx0,size(X,1),1);
end
tsrmbloop.varsid = setdiff(1:p0, colmis);
tsrmbloop.X = X.*repmat(sx0,size(X,1),1);
tsrmbloop.Xallvars = Xremvars;
tsrmbloop.L = L;
tsrmbloop.Xrec = Xrec;
tsrmbloop.SPE = SPE;
tsrmbloop.E = E;
tsrmbloop.T2 = T2;
tsrmbloop.pcamodel = pcamodel;
tsrmbloop.obsid = obsid;
tsrmbloop.rwloop = rwloop;
tsrmbloop.cwloop = cwloop;
tsrmbloop.pcwfilt = rw_cw_filt;
tsrmbloop.loopanglep = loopanglep;
tsrmbloop.loopdiff = loopdiff;
tsrmbloop.it = it;
tsrmbloop.cspe = 2.5*pcamodel.limspe;
tsrmbloop.ct2 = 2.5*pcamodel.limt2;
tsrmbloop.ct2_ph2 = 2.5*pcamodel.limt2_ph2;
tsrmbloop.h.maxp = hP;
tsrmbloop.h.diff = hD;
tsrmbloop.h.nobs = nobs;
tsrmbloop.h.ncells = ncells;
end