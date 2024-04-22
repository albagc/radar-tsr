function [mse] = erroreval(Xtrue, Ximp, rindex, cindex)
%--------------------------------------------------------------------------
% erroreval - Evaluate mean squared error between true and imputed data
%--------------------------------------------------------------------------
%   [mse] = erroreval(Xtrue, Ximp, rindex, cindex) evaluates the mean
%   squared error (MSE) between the true data Xtrue and the imputed data
%   Ximp, considering only the specified rows and columns.
%
%   Input:
%   - Xtrue: True data matrix.
%   - Ximp: Imputed data matrix.
%   - rindex: Index or logical mask indicating which rows to consider for
%             MSE calculation. Default is true(size(Ximp,1),1) (all rows).
%   - cindex: Index or logical mask indicating which columns to consider
%             for MSE calculation. Default is true(size(Ximp)) (all
%             columns).
%
%   Output:
%   - mse: Mean squared error between the true and imputed data.
%
%   Example:
%   Xtrue = [1, 2, 3; 4, 5, 6];
%   Ximp = [1, NaN, 3; NaN, 5, 6];
%   rindex = [true, false];
%   cindex = [true, true, false];
%   mse = erroreval(Xtrue, Ximp, rindex, cindex);
%
%--------------------------------------------------------------------------

arguments
   Xtrue double
   Ximp double
   rindex {double, logical} = true(size(Ximp, 1), 1)
   cindex {double, logical} = true(size(Ximp))
end

M_rin = cindex;

if length(rindex) ~= size(cindex, 1)
    M_rin = cindex(rindex, :);
end

if sum(~rindex) > 0
    M_rin(~rindex, :) = false(size(cindex(~rindex, :)));
end

Ximp_rin = Ximp(M_rin);
Xtrue_rin = Xtrue(M_rin);
mse = mean((Ximp_rin - Xtrue_rin).^2);

end
