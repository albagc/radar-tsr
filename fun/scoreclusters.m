function [cltagE, cl_model] = scoreclusters(T, max_clusters, model, acompmode, A, kcompmode, kdist)
%--------------------------------------------------------------------------
% scoreclusters - Perform clustering analysis on scores matrix
%--------------------------------------------------------------------------
%
% Inputs:
% T: Scores matrix
% max_clusters: Maximum number of clusters to consider (default: 10)
% model: Clustering model to use ('gmm' for Gaussian Mixture Model,
% 'knn' for K-means clustering) (default: 'gmm')
% acompmode: Automatic computation mode for the number of components in
% GMM ('auto' or 'manual') (default: 'auto')
% A: Number of components to use in GMM (default: min([10, size(E,2)]))
% kcompmode: Computation mode for the number of clusters in K-means
% ('auto' or 'manual') (default: 'auto')
% kdist: Distance metric used in K-means clustering (default: 'sqeuclidean')
%
% Outputs:
% cltagE: Cluster tags assigned to each observation
% cl_model: Cluster model containing information about the clusters
%
% Example:
% T = [1, 2; 3, 4; 5, 6]; % Scores matrix
% max_clusters = 5; % Maximum number of clusters
% model = 'gmm'; % Clustering model
% acompmode = 'auto'; % Automatic computation mode for GMM
% A = 3; % Number of components in GMM
% kcompmode = 'manual'; % Computation mode for K-means
% kdist = 'sqeuclidean'; % Distance metric for K-means
% [cltagE, cl_model] = scoreclusters(T, max_clusters, model, acompmode, A, 
% kcompmode, kdist); % Perform clustering analysis
arguments
    T double
    max_clusters double = 10
    model {string,char} = "gmm"
    acompmode {string,char} = "auto"
    A double = min([10, size(E,2)])
    kcompmode {string,char} = "auto"
    kdist = "sqeuclidean"
end
    % Error palette %%%%
    err2 = " ill-conditioned covariance created at iteration 2.";
    % Perform k-means clustering on scores matrix and optimize number of clusters
    % X is already scaled at the beginnning of the RadarTSR algorithm
    
    % Initialize variables
    scores_matrix = T;
    a_wcss = zeros(max_clusters, 1);
    wcss = zeros(max_clusters, 1);
    if strcmp(model, "knn")
        if max_clusters > 1
            if strcmp(kcompmode,"manual")
                figure,plotmatrix(T(:,1:A)), title("K-means step - Residual scores distributions ")
                for k = 1:max_clusters
                    % Perform k-means clustering
                    [cluster_labels, cluster_centers] = kmeans(scores_matrix, k,"Distance","sqeuclidean",'Replicates',5);
                    % Compute within-cluster sum of squares
                    a_wcss(k) = mean(sum((scores_matrix - cluster_centers(cluster_labels, :)).^2,2));
                    wcss(k) = sum((scores_matrix - cluster_centers(cluster_labels, :)).^2,"all");
                end
                figure, plot(wcss, 'o-','DisplayName',"Total within group sums of squares"),
                hold on, grid on, xlabel("K (number of clusters)"), title("WCSSQ vs. K")
                %plot(a_wcss,"sq-",'DisplayName',"Average within group sums of squares")
                %legend
                ylabel("within-cluster sums of squares")
                optimal_num_clusters = input("Enter optimal K: ");
            else
                % Find the optimal number of clusters using the elbow method
                eval_num_clusters = evalclusters(scores_matrix, 'kmeans', 'CalinskiHarabasz', 'KList', 1:max_clusters);
                if eval_num_clusters.OptimalK >= size(scores_matrix,1)
                    optimal_num_clusters = 0;
                else
                    optimal_num_clusters = eval_num_clusters.OptimalK;
                end
            end
        else
            optimal_num_clusters =1;
        end

        if optimal_num_clusters > 0
             % Perform k-means clustering with optimal number of clusters
            [cluster_labels, cluster_centers] = kmeans(scores_matrix, optimal_num_clusters);
            u_K = unique(cluster_labels);
            mDlimit = nan(length(u_K),1);
            for k = 1:length(u_K)
                mD = sum(pdist2(scores_matrix(cluster_labels==u_K(k),:), cluster_centers(k,:)).^2,2);
                mDlimit(k) = max(mD);
            end
            cl_model.mahal_lim = mDlimit;
        else
            cluster_labels = 0;
        end
    elseif strcmp(model,"gmm")
        if A==0
            A = size(T,2);

            mahalDist = nan(1, A);

            for k = 1:max_clusters
                try
                    gm = fitgmdist(scores_matrix, k,"CovarianceType","diagonal","SharedCovariance",false,"Start","plus");
                    cluster_labels = cluster(gm,scores_matrix);
                    mD = mahal(gm,scores_matrix);
                    ind = sub2ind(size(mD),(1:size(mD,1))',cluster_labels);
                    mahalDist(k) = sum(mD(ind).^2);
                    cl_model.gmm_model = gm;
                catch ME
                end
            end
            if strcmp(kcompmode,"manual")
                figure,plotmatrix(T), title("GMM step - Scores distributions ")
                figure,plotmatrix(T(:,1:A)), title("GMM step - Residual scores distributions ")
                figure,plot(mahalDist), title("Mahalanobis distance vs. number of components for GMM")
                optimal_num_clusters = input("Enter optimal K (number of components in the mixed gaussian model): ");
            elseif strcmp(kcompmode,"auto")
                [~,optimal_num_clusters] = min(mahalDist,[],"omitmissing");
            else
                optimal_num_clusters = A;
            end
        else
            optimal_num_clusters = A;
        end

        % Perform GMM clustering with optimal number of clusters
        cluster_labels = zeros(size(scores_matrix,1),1);
        
        try
            gm = fitgmdist(scores_matrix, optimal_num_clusters);
            cluster_labels = cluster(gm,scores_matrix);
            cl_model.gmm_model = gm;
            u_K = unique(cluster_labels);
            mDlimit = nan(length(u_K), 1);
            for k = 1:length(u_K)
                mD = mahal(gm,scores_matrix);
                mDlimit(k) = max(mD(cluster_labels==u_K(k)));
            end
            cl_model.mahal_lim = mDlimit;
        catch ME
            if strcmp(ME.identifier, err2)
               optimal_num_clusters = input(strcat("Enter different K (<", string(optimal_num_clusters),"): "));
               gm = fitgmdist(scores_matrix, optimal_num_clusters);
               cluster_labels = cluster(gm,scores_matrix);
               cl_model.gmm_model = gm;
               u_K = unique(cluster_labels);
               mDlimit = nan(length(u_K), 1);
               for k = 1:length(u_K)
                   mD = mahal(gm,scores_matrix);
                   mDlimit(k) = max(mD(cluster_labels==u_K(k)));
               end
               cl_model.mahal_lim = mDlimit;
            end
        end
    end

    if numel(cluster_labels) > 1
        [~,~,cltag_uns] = unique(cluster_labels,'rows');
        cltagE = sortclass(cltag_uns) + 1;
        cltagE = checkminobs(cltagE, 5) + 1;
    
        u_cltag = unique(cltagE(~isnan(cltagE)));
    
        cluster_centers = zeros(length(u_cltag), size(scores_matrix,2));
        for k = 1:length(u_cltag)
            cluster_centers(k,:) = mean(scores_matrix(cltagE == u_cltag(k)));
        end
    
        cl_model.cluster_centers = cluster_centers;
        cl_model.optimal_num_clusters = optimal_num_clusters;
        cltagE(isnan(cltagE)) = 0;
    else
        cltagE = zeros(size(scores_matrix,1),1);
        cl_model.optimal_num_clusters = optimal_num_clusters;
    end
end
