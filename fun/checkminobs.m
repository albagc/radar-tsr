function s_j_cltag = checkminobs(cltags, minobs)
%--------------------------------------------------------------------------
% checkminobs - Check and assign cluster tags based on minimum observations
%--------------------------------------------------------------------------
%   s_j_cltag = checkminobs(cltags, minobs) checks the cluster tags in
%   'cltags' and assigns a new set of cluster tags based on the minimum
%   number of observations 'minobs' for each cluster.
%
%   Input:
%   - cltags: Cluster tags. An array or vector representing cluster
%             assignments.
%   - minobs: Minimum number of observations required for a cluster to be
%             considered valid. Clusters with fewer observations will be
%             removed. Default is 3.
%
%   Output:
%   - s_j_cltag: Updated cluster tags with clusters reassigned based on
%                the minimum number of observations.
%
%   Example:
%   cltags = [1 1 2 2 3 3 3 4 4 4 4]; % Cluster tags
%   minobs = 3;
%   s_j_cltag = checkminobs(cltags, minobs);
%
%--------------------------------------------------------------------------

% Find unique cluster tags
utags = unique(cltags);

% Sort unique cluster tags in ascending order
u_jstep = sort(utags, 'ascend');

% Initialize matrix to store cluster tags and their respective counts
nx_cltags = [u_jstep, zeros(length(u_jstep), 1)];

% Count the number of observations for each cluster tag
for jntag = 1:size(nx_cltags, 1)
    nx_cltags(jntag, 2) = sum(cltags == nx_cltags(jntag, 1));
end

% Sort cluster tags based on the count of observations in descending order
[~, istags] = sort(nx_cltags(:, 2), 'descend');
nx_cltags = nx_cltags(istags, :);

% Delete cluster tags with fewer observations than the specified minimum
del_tag = nx_cltags(:, 2) <= minobs;
nx_cltags(del_tag, :) = [];

% Assign new cluster tags based on the updated order and count
s_j_cltag = nan(size(cltags));
for jntag = 1:size(nx_cltags, 1)
    s_j_cltag(cltags == nx_cltags(jntag, 1)) = jntag;
end

% Sort the new cluster tags in ascending order
s_j_cltag(~isnan(s_j_cltag)) = sortclass(s_j_cltag(~isnan(s_j_cltag)));

end
