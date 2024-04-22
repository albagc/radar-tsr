function [svec] = scalemask(X, M, est)
%--------------------------------------------------------------------------
% scalemask - Calculate scaled standard deviation using a mask
%--------------------------------------------------------------------------
%
% Inputs:
% X: Data matrix
% M: Mask matrix specifying elements to exclude from the calculation
% (default: isnan(X))
% est: Estimation method for scaled standard deviation calculation
% ('ls' for least squares estimation, 'rob' for robust estimation)
% (default: 'ls')
%
% Outputs:
% svec: Scaled standard deviation vector calculated from the available data
%
% Example:
% X = [1, NaN, 3; 4, 5, 6; NaN, 8, 9];
% M = isnan(X);
% est = 'ls';
% svec = scalemask(X, M, est);
arguments
   X double
   M {double,logical} = isnan(X)
   est {string,char} = 'ls'
end
% Calculate mean only using mask elements set to 1
% Imputs mean of available data
Z = X;
[n,p] = size(X);
if strcmp(est,'ls')
    Z(M) = 0;
    mvec = locmask(Z,M,est);
    difmed = Z-mvec;
    difmed(M) = 0;
    svec = sqrt(sum(difmed.^2)./(n-1-sum(M)));
elseif strcmp(est,'rob')
    svec = zeros(1,size(Z,2));
    for j = 1:p
        z1 = Z(M(:,j)==0,j);
        svec(j) = robscale(z1, 2.5, 0.844472, 1e-12);
    end
end
end
