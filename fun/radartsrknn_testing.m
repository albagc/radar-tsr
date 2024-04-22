function test_output = radartsrknn_testing(X, A, ...
    acompmode, Mtrue, mtrue, xscale, alpha, maxit, told, tola, spelimcal, t2limcal, ...
    minpctge, n_min_cluster)
%--------------------------------------------------------------------------
% radartsrknn - RADAR-TSR algorithm for missing data imputation and outlier 
% detection
%--------------------------------------------------------------------------
%
% Inputs:
%   X: Data matrix with observations in rows and variables in columns.
%   A: Number of principal components.
%   acompmode: Mode for selecting the number of principal components.
%   Mtrue: Indicator matrix for known missing values (optional, default is all false).
%   mtrue: Indicator vector for known row outliers (optional, default is all false).
%   xscale: Logical value indicating whether to perform feature scaling (optional, default is false).
%   alpha: Significance level for outlier detection (optional, default is 0.01).
%   maxit: Maximum number of iterations for iterative PCA (optional, default is 50).
%   told: Tolerance for convergence in iterative PCA (optional, default is 1e-6).
%   tola: Tolerance for convergence in TSR-PCA (optional, default is 5e-3).
%   spelimcal: Type of SPE limit calculation (optional, default is 'nomikmac').
%   t2limcal: Type of T2 limit calculation (optional, default is 'theor').
%   minpctge: Minimum percentage of variance explained by each principal component (optional, default is 0.2).
%   n_min_cluster: Minimum number of observations in a cluster for radar-based outlier detection (optional, default is 5).
%
% Outputs:
%   Ximp: Imputed data matrix.
%   cwpred: Cell-wise outlier indicator matrix.
%   rwpred: Row-wise outlier indicator vector.
%   cltagE: Cluster tags based on radar-based outlier detection.
%   modeloutput: Additional model output and diagnostics.
%
% Example:
%   X = load('data.mat');
%   A = 10; % Number of principal components
%   acompmode = "auto"; % Mode for selecting the number of principal components
%   Mtrue = false(size(X)); % Indicator matrix for known missing values
%   mtrue = false(size(X,1),1); % Indicator vector for known row outliers
%   xscale = false; % Perform feature scaling
%   alpha = 0.01; % Significance level for outlier detection
%   maxit = 50; % Maximum number of iterations for iterative PCA
%   told = 1e-6; % Tolerance for convergence in iterative PCA
%   tola = 5e-3; % Tolerance for convergence in TSR-PCA
%   spelimcal = 'nomikmac'; % Type of SPE limit calculation
%   t2limcal = 'theor'; % Type of T2 limit calculation
%   minpctge = 0.2; % Minimum percentage of variance explained by each principal component
%   n_min_cluster = 5; % Minimum number of observations in a cluster for radar-based outlier detection
%   [Ximp, cwpred, rwpred, cltagE, modeloutput] = radartsrknn(X, A, acompmode, Mtrue, mtrue, xscale, alpha, maxit, told, tola, spelimcal, t2limcal, minpctge, n_min_cluster);
%
% Note: The X input should be a data matrix with observations in rows and variables in columns.

arguments
    X double
    A {double,char}
    acompmode {string,char} = "auto"
    Mtrue = false(size(X))
    mtrue = false(size(X,1),1)
    xscale (1,1)logical = false
    alpha (1,1)double = 0.01
    maxit (1,1)double = 50
    told (1,1)double = 1e-6
    tola (1,1)double = 5e-3
    spelimcal {string,char} = 'nomikmac'
    t2limcal {string,char} = 'theor'
    minpctge (1,1)double = 0.2
    n_min_cluster (1,1)double = 5

end
% Error palette %%%%
err1 = " missing for all observations - deleted from the loop";
err2 = " ill-conditioned covariance created at iteration 2.";
%%%%%%%%%%%%%%%%%%%%
[n,k] = size(X);
Mmd = isnan(X);
X0 = X;
% Compute scales
if xscale
    sx0 = scalemask(X,Mmd,"rob");
    prepro = "autosc";
else
    sx0 = ones(1,k);
    prepro = "cent";
end
%
% Standardization 
Xnai = X./repmat(sx0, n, 1);
%%
colmis = find(or(sum(isnan(Xnai))>=(size(Xnai,1)-2), sum(isinf(Xnai))>0));
if ~isempty(colmis)
    remvar = true;
    X(:,colmis) = [];
    Xnai(:,colmis) = [];
    Mmd(:,colmis) = [];
    sx0(:,colmis) = [];
    disp(strcat("Variable(s): ", strjoin(string(colmis), ", "), err1)) 
else
    remvar = false;
end
%% Step 1: Big univariate outliers (potential cellwise outliers) - - - - - - - - - -
[cwb, uvfilt] = univarfilter(Xnai, alpha, 'rob', Mmd, acompmode);
%% Step 2: Rows with too many univariate outliers stored in Irw0
[m0, rw0filt] = rowfiltcw(cwb, Mmd, alpha, "rob");
Xnai_wocw = Xnai;
Xnai_wocw(cwb) = nan;
%% Step 3: TSR for iterative PCA MB + iterative detection of outliers
Mall = logical(Mmd + cwb);
% - - - - - - - - - - - - - - - - - - - - -
%% TEST -- influence flunction
[mbloop, SPE, T2, pca0, test_output] = pcambloop_test(Xnai_wocw(~m0,:), A, Mmd(~m0,:), Mall(~m0,:), find(~m0), alpha, ...
    'rob', tola, told, maxit, true, spelimcal, t2limcal, false, minpctge, acompmode, mtrue(~m0), Xnai(~m0,:));
% Store testing results for convergence and for outliers' detection
test_output.Mcw = cwb;
test_output.Mrwpre = m0;
test_output.hP = mbloop.h.maxp;
test_output.hD = mbloop.h.diff;
test_output.hNobs = mbloop.h.nobs;
test_output.hNcells = mbloop.h.ncells;
test_output.c_hP = tola;
test_output.c_hD = told;
test_output.c_hIt = maxit;

Mcw_loop = false(size(Xnai_wocw(:,mbloop.varsid)));
Mcw_loop(~m0,:) = mbloop.cwloop;
m_speloop = false(n,1);
m_speloop(~m0) = mbloop.rwloop;

%% Step 4: NA-imputation with TSR for PCA-ME, projection and detection 
Xna_imp = tsrme(Xnai(:, mbloop.varsid),pca0);
[Xrec_na_i, ~, ~ , SPE_na_i, T2_na_i] = pcame(Xna_imp, pca0);
R_na_i = Xnai(:, mbloop.varsid) - Xrec_na_i;
Mcw_na = univarfilter(R_na_i, alpha, 'rob', Mmd(:, mbloop.varsid), acompmode);
Mcwfull = or(or(cwb(:, mbloop.varsid), Mcw_na), Mcw_loop);
m_cw_na = (sum(Mcwfull,2)/k) > mbloop.pcwfilt.ccwpre;
[~, m_naspe] = rw_detect(SPE_na_i, T2_na_i, pca0.limspe, pca0.limt2, 2.5, 2.5, alpha);

% Rows left out with high na-imputed SPE and a high p_cw + rows detected as cell-imputed SPE outliers in the PCA-MB loop
m_na = or(and(m_cw_na, m_naspe), m_speloop); % 
Mcw_full = or(cwb(:, mbloop.varsid), Mcw_na);
test_output.Pstep4 = pca0.P;

%% Step 5: Estimation with TSR for PCA-MB without any rowwise nor cellwise outliers
Xcell = Xnai(:, mbloop.varsid);
Xcell(Mcw_full) = nan; % Make sure there are enough observations to refit the model
Xcell(Mcw_full) = nan; % Make sure there are enough observations to refit the model
if any(sum(Mcw_full)==sum(~m_na))
    Xcell_imp = tsrme(Xcell,pca0);
    [~,~,~,~,pca0] = pcamb(Xcell_imp(~m_na,:), pca0.ncomp, alpha, 'cent', "rob", spelimcal, t2limcal);
else
    [~, ~, ~,  ~,  ~,  ~, pca0] = pcambtsr_diffP(Xcell(~m_na,:), pca0.ncomp, false, maxit, tola, alpha);
end
test_output.Pfinal = pca0.P;
% Cell-imputation with TSR for PCA-ME, projection and detection
X_cellimp = tsrme(Xcell, pca0);
mbloop.varsid = 1:size(Xnai_wocw(:, mbloop.varsid),2);
[X_cellimp_rec, ~, ~, SPEci, T2ci] = pcame(X_cellimp, pca0); 
m_cell_i = rw_detect(SPEci, T2ci, pca0.limspe, pca0.limt2, 2.5, 2.5, alpha); %m_cell_i_2
R_cell_i = Xnai(:,mbloop.varsid) - X_cellimp_rec;
Mcw = univarfilter(R_cell_i, alpha, 'rob', Mmd(:, mbloop.varsid), acompmode);
% New detection based on joining NA imputed rows and Cell imputed rows
m_cell = or(m_na, m_cell_i);

% Last imputation and projection according to the detection of outliers
Mcw(m_cell,:) = false;
Xcell = Xnai(:,mbloop.varsid);
Xcell(Mcw) = nan;
% Impute, project and detect 
X_cellimp = tsrme(Xcell, pca0);
[X_cellimp_rec, Tci, Eci, SPEci, T2ci] = pcame(X_cellimp, pca0);
[Mrw, rw_spe, rw_t2] = rw_detect(SPEci, T2ci, pca0.limspe, pca0.limt2, 2.5, 2.5, alpha);
[~, rw_cwfilt] = rowfiltcw(Mcw(~Mrw,:), Mmd(~Mrw,:), alpha);
Mcw = univarfilter(R_cell_i, alpha, 'rob', Mmd(:, mbloop.varsid), acompmode);
Mcw(Mrw,:) = false;
%% Return clean model terms 
% Rowwise outliers classification
rw_extr = and(~rw_spe, rw_t2);

% Update cut-offs for rowwise outliers:
rwCellImp.cspe = 2.5*pca0.limspe;
rwCellImp.ct2 = 2.5*pca0.limt2_ph2;
rwCellImp.limits_t = pca0.limits_t;
rwCellImp.pcw = rw_cwfilt;

% Undo initial scaling
pca0.prepro = prepro;
pca0.s = sx0(mbloop.varsid);
pca0.m = pca0.m.*sx0(mbloop.varsid);
mbloop = rmfield(mbloop, "pcamodel");
mbloop.cl_0.pcamodel = pca0;
X_cellimp = X_cellimp.*repmat(sx0(mbloop.varsid),size(X_cellimp,1),1);
X_cellimp_rec = X_cellimp_rec.*repmat(sx0(mbloop.varsid),size(X_cellimp_rec,1),1);
X_fullimp_rec_cluster = X_cellimp_rec;
X_fullimp = X_cellimp;

%% Radar step for orthogonal outlying patterns (detection of clusters by scores (radar on scores))
do_cluster = false;
if and(any(m_na), sum(rw_spe) > n_min_cluster)
    [~,~,~,~,pcaE] = pcamb(Eci(rw_spe,:), 0, alpha, 'cent', "rob", spelimcal, t2limcal, 0.1, acompmode, "eachpc");
    if any(pcaE.vexp > 0.1)
        if pcaE.ncomp == 1
            y_max = 1;
        else
            y_max = min([sum(rw_spe),floor(sum(rw_spe)/n_min_cluster),10]);
        end
        [~,Teci] = pcame(Eci(rw_spe,:),pcaE);
        [cltagK, cluster_models] = scoreclusters(Teci, y_max, "knn", acompmode, pcaE.ncomp, acompmode);
        cltagE = double(rw_spe);
        cltagE(rw_spe) = cltagK;
        cltagE(isnan(cltagE)) = 0;
        if cluster_models.optimal_num_clusters > 0
            do_cluster = true;
            for k = 1:size(cluster_models.cluster_centers,1)
                vars_clk = 1:size(Xnai,2);
                vars_clk(sum(isnan(Xnai))==size(Xnai,1)) = [];
                X_imp_vblesref = Xnai(cltagE==k,vars_clk);
                X_rec_vblesref = Xnai(cltagE==k,vars_clk);
                [X_k_imp, ~, ~,  ~,  ~,  ~, pcamodelxcl] = pcambtsr_diffP(X_imp_vblesref, 0, false, maxit, tola, alpha, acompmode);
                [X_rec_cluster ,~,~,SPEk,T2k] = pcame(X_k_imp, pcamodelxcl);
                % Detect rowwise
                [~, rw_spe_k] = rw_detect(SPEk, T2k, pcamodelxcl.limspe, pcamodelxcl.limt2, 2.5, 2.5, alpha);
                ind_k = find(cltagE==k);
                ind_kin = ind_k;
                ind_kin(rw_spe_k) = [];
                ind_k(~rw_spe_k) = [];
                cltagE(ind_k) = 0;
                X_imp_vblesref(:,vars_clk) = X_k_imp.*repmat(sx0(vars_clk),size(X_k_imp,1),1);
                X_fullimp(ind_kin,ismember(mbloop.varsid, vars_clk)) = X_imp_vblesref(~rw_spe_k, ismember(vars_clk, mbloop.varsid));
                X_rec_vblesref(:,vars_clk) = X_rec_cluster.*repmat(sx0(vars_clk),size(X_k_imp,1),1);
                X_fullimp_rec_cluster(ind_kin,ismember(mbloop.varsid, vars_clk)) = X_rec_vblesref(~rw_spe_k, ismember(vars_clk, mbloop.varsid));
                pcamodelxcl.prepro = prepro;
                pcamodelxcl.s = sx0(vars_clk);
                pcamodelxcl.m = pcamodelxcl.m.*sx0(vars_clk);
                mbloop.(strcat("cl_", string(k))).pcamodel = pcamodelxcl;
            end            
        end
    else
        cltagE = zeros(size(rw_spe));
    end
else
    X_fullimp = X_cellimp;
    cltagE = zeros(n,1);
end
[~,~,cltag_uns] = unique(cltagE,'rows');
cltags = sortclass(cltag_uns);

n_clusters.all = unique(cltags);
n_clusters.rw_E = unique(cltagE);

cl_tags.all = cltags;
cl_tags.rw_E = cltagE;
cl_tags.rw_T = rw_extr;

X_naimp = X(:, mbloop.varsid);
X_naimp(Mmd(:, mbloop.varsid)) = X_fullimp(Mmd(:, mbloop.varsid));

%% Get projection metrics with reference model
[Xrec_fi,Tfi,Efi,SPEfi,T2fi,T2matrixfi] = pcame(X_fullimp, pca0);
[Xrec_ci,Tci,Eci,SPEci,T2ci,T2matrixci] = pcame(X_cellimp, pca0);
[Xrec_na_i,Tnai,Enai,SPE_na_i,T2_na_i,T2matrixnai] = pcame(X_naimp, pca0);
[~, erfilter] = univarfilter(Eci, alpha, 'rob', Mmd(:, mbloop.varsid), acompmode, false, true);
rwpred = Mrw;
cwpred = Mcw;
%cwpred_full = cwfull;
modeloutput = struct();

% X mdimp only in not outlying rows
Xmdimp = X(:, mbloop.varsid);
Xmdimp(~rwpred,:) = X_naimp(~rwpred,:);

Xfullimp_rec_cluster = X_fullimp_rec_cluster;
Xcellimp_rec_cluster = X_fullimp_rec_cluster;
Xcellimp_rec_cluster(cltagE == 0,:) = Xrec_ci(cltagE == 0,:);

Ximp = struct();
Ximp.Xclrows_md = Xmdimp;
Ximp.Xci = X_cellimp; % NAs in all, CW in not oytlying rows
Ximp.Xnaimp = X_naimp; % only NAs
Ximp.Xfi = X_fullimp; % NAs and CW in all rows (even outlying cells in RW)

modeloutput.Ximp = Ximp;
modeloutput.remvar.bool = remvar;
vini = 1:size(X0,2);
vini(colmis) = [];
vini = vini(mbloop.varsid);

Ximp.XBigOutput = struct();
Ximp.XBigOutput.Xclrows_md = X0;
Ximp.XBigOutput.Xclrows_md(:,vini) = Ximp.Xclrows_md;
Ximp.XBigOutput.Xci = X0;
Ximp.XBigOutput.Xci(:,vini) = Ximp.Xci;
Ximp.XBigOutput.Xnai = X0;
Ximp.XBigOutput.Xnai(:,vini) = Ximp.Xnaimp;
Ximp.XBigOutput.Xfi = X0;
Ximp.XBigOutput.Xfi(:,vini) = Ximp.Xfi;

modeloutput.remvar.colmis = setdiff(1:size(X0,2), vini);
modeloutput.remvar.colin = vini;
modeloutput.PCA_Model = pca0;
if do_cluster
    modeloutput.FP_E_clusters = cluster_models.cluster_centers;
    modeloutput.cluster_models = cluster_models;
end
% Parameters of the preprocessning (for all observations)

modeloutput.PCA_Model_allclusters = mbloop;
modeloutput.NAimp.Xnairec = Xrec_na_i;
modeloutput.NACWimp.Xcirec = Xrec_ci;
modeloutput.FULLimp.Xfirec = Xrec_fi;
modeloutput.Xfirec_allclusters = Xfullimp_rec_cluster;
modeloutput.Xcirec_allclusters = Xcellimp_rec_cluster;

% FUll imp
modeloutput.FULLimp.Tfi = Tfi;
modeloutput.FULLimp.Efi = Efi;
modeloutput.FULLimp.Efi_na = Efi;
modeloutput.FULLimp.Efi_na(Mmd(:,mbloop.varsid)) = nan;
modeloutput.FULLimp.SPEfi = SPEfi;
modeloutput.FULLimp.T2fi = T2fi;
modeloutput.FULLimp.T2matrixfi = T2matrixfi;

% CW + NA imp
modeloutput.NACWimp.Tci = Tci;
modeloutput.NACWimp.Eci = Eci;
modeloutput.NACWimp.Eci_na = Eci;
modeloutput.NACWimp.Eci_na(Mmd(:,mbloop.varsid)) = nan;
modeloutput.NACWimp.SPEci = SPEci;
modeloutput.NACWimp.T2ci = T2ci;
modeloutput.NACWimp.T2matrixci = T2matrixci;

% NA imp
modeloutput.NAimp.Tnai = Tnai;
modeloutput.NAimp.Enai = Enai;
modeloutput.NAimp.SPEnai = SPE_na_i;
modeloutput.NAimp.T2nai = T2_na_i;
modeloutput.NAimp.T2matrixnai = T2matrixnai;

% Filters
modeloutput.cwfilt.uv = uvfilt;
modeloutput.cwfilt.erec = erfilter;
modeloutput.rw0filt = rw0filt;
modeloutput.rwfilt = rwCellImp;

modeloutput.nacells = Mmd;
modeloutput.nclusters = sum(n_clusters.all>0);
modeloutput.cluster_tag = cl_tags;
modeloutput.rwout = rwpred;
modeloutput.cwout = cwpred;
%modeloutput.cwout_fullimp = cwpred_full;
modeloutput.convergence.diffloop = mbloop.loopdiff;
modeloutput.convergence.minangle = mbloop.loopanglep;
modeloutput.convergence.niter = mbloop.it;
modeloutput.convergence.h = mbloop.h;
modeloutput.A = pca0.ncomp;
end
