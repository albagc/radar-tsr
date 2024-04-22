function [remX, xinfo] = checkds(X, epsilon, nmaxcat, numboxplot)
%--------------------------------------------------------------------------
% checkds - Check and preprocess dataset
%--------------------------------------------------------------------------
%   [remX, xinfo] = checkds(X, epsilon, nmaxcat, numboxplot) performs
%   various checks and preprocessing steps on the input dataset X.
%
%   Input:
%   - X: Input dataset. It can be a matrix or a structure. If X is a matrix,
%        it will be converted to a table. If X is a structure, it will be
%        converted to a table.
%   - epsilon: Threshold for zero variance. Variables with variance less
%              than epsilon will be removed. Default is 0.
%   - nmaxcat: Maximum number of categories for categorical variables.
%              Variables with more than nmaxcat categories will be removed.
%              Default is 10.
%   - numboxplot: Flag indicating whether to plot box plots for numerical
%                 variables. If numboxplot is true, box plots will be
%                 displayed. Default is false.
%
%   Output:
%   - remX: Processed dataset after removing variables with zero variance
%           and NaN values.
%   - xinfo: Information about the dataset.
%     - skewnum: Skewness values for numerical variables.
%     - numvars: Names of numerical variables.
%     - catvars: Names of categorical variables after converting them to
%                dummy variables.
%     - delallNA: Logical vector indicating variables with NaN values.
%     - del0var: Logical vector indicating variables with zero variance.
%     - delvars_tommanycat: Variables removed due to excessive categories.
%     - dummy_original_corresp: Original categorical variables and their
%                               corresponding dummy variables.
%     - ind_numvars: Indices of numerical variables in remX.
%     - ind_catvars: Indices of categorical variables in remX.
%
%   Example:
%   X = rand(100, 5); % Random numerical dataset
%   [remX, xinfo] = checkds(X, 0.1, 5, true);
%
%   X = struct('A', [1; 2; 3], 'B', {'Red'; 'Green'; 'Blue'}); % Struct dataset
%   [remX, xinfo] = checkds(X);
%
%--------------------------------------------------------------------------

% Convert matrix to table if X is a matrix and not a table
if ismatrix(X) && ~istable(X)
    X = array2table(X);
    
% Convert structure to table if X is a structure
elseif isstruct(X)
    X = struct2table(X);
end

% Check variable types and convert text variables to appropriate types
[X, delvars_tommanycat] = checktxtvars(X);

% Convert categorical variables to dummy variables
[Xdummies, ind_catvars] = makedummy(X, nmaxcat);

% Remove original categorical variables from the table
dumX = X;
dumX(:,ind_catvars.OriginalNames) = [];

% Plot box charts for numerical variables if numboxplot is true
if numboxplot
    px = floor(sqrt(size(dumX,2)));
    py = ceil(size(dumX,2)/px);
    figure,
    for i = 1:size(dumX,2)
        subplot(px, py, i),
        boxchart(dumX{:,i},'Orientation','horizontal'),
        title(dumX.Properties.VariableNames{i}),
    end
end

% Calculate skewness and store variable names
xinfo.skewnum = skewness(dumX{:,:});
xinfo.numvars = string(dumX.Properties.VariableNames);
xinfo.catvars = string(Xdummies.Properties.VariableNames);

% Append dummy variables to the table
dumX = [dumX, Xdummies];

% Remove variables with NaN values based on a threshold
[n,~] = size(dumX);
xinfo.delallNA = sum(isnan(dumX{:,:})) > (n - 3);

% Remove variables with zero variance based on epsilon threshold
xinfo.del0var = var(dumX{:,:}, 0, "omitnan") < epsilon;

% Store additional information in xinfo
xinfo.delvars_tommanycat = delvars_tommanycat;
xinfo.dummy_original_corresp = ind_catvars;

% Remove variables with NaN values and zero variance from the dataset
remX = dumX(:,~logical(xinfo.delallNA + xinfo.del0var));

% Store indices of numerical and categorical variables
xinfo.ind_numvars = 1 : length(xinfo.numvars);
xinfo.ind_catvars = (length(xinfo.numvars) + 1) : (size(dumX,2));

end
