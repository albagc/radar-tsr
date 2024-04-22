function stag = sortclass(tags)
%--------------------------------------------------------------------------
% sortclass - Sort cluster tags in descending order based on tag frequency
%--------------------------------------------------------------------------
%
% Inputs:
%   tags: Cluster tags assigned to each observation
%
% Output:
%   stag: Sorted cluster tags with values ranging from 0 to (number of tags - 1)
%
% Example:
%   tags = [1; 2; 2; 3; 1; 2]; % Cluster tags
%   stag = sortclass(tags); % Sort the cluster tags

    u_cltag_uns = sort(unique(tags), 'ascend');
    nx_cltags = [u_cltag_uns, zeros(length(u_cltag_uns), 1)];
    
    % Count the frequency of each tag
    for jntag = 1:size(nx_cltags, 1)
        nx_cltags(jntag, 2) = sum(tags == nx_cltags(jntag, 1));
    end
    
    % Sort the tags based on their frequency in descending order
    [~, istags] = sort(nx_cltags(:, 2), 'descend');
    nx_cltags = nx_cltags(istags, :);
    
    s_j_cltag = nan(size(tags));
    
    % Assign sorted tag values to observations
    for jntag = 1:size(nx_cltags, 1)
        s_j_cltag(tags == nx_cltags(jntag, 1)) = jntag;
    end
    
    % Adjust tag values to range from 0 to (number of tags - 1)
    stag = s_j_cltag - 1;

end
