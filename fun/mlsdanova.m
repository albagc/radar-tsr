function [lsd_tbl] = mlsdanova(data, vary, fac1, fac2, rndfac, disptbl)
%--------------------------------------------------------------------------
% mlsdanova - Perform a multi-level single difference (MLSD) analysis of
% variance (ANOVA)
%--------------------------------------------------------------------------
%   [lsd_tbl] = mlsdanova(data, vary, fac1, fac2, rndfac, disptbl) performs
%   a multi-level single difference (MLSD) analysis of variance (ANOVA) on
%   the specified data. It calculates the group means and LSD (least
%   significant difference) intervals based on the ANOVA results.
%
%   Inputs:
%   - data: Data table or array containing the observations.
%   - vary: Name or index of the dependent variable in the data.
%   - fac1: Name or index of the first fixed factor in the data.
%   - fac2: Name or index of the second fixed factor in the data.
%   - rndfac: Name or index of the random factor in the data.
%   - disptbl: Display table option. 'off' to hide the ANOVA table
%              (optional, default is 'off').
%
%   Output:
%   - lsd_tbl: Table containing the group means and LSD intervals based on
%              the MLSD ANOVA.
%
%   Example:
%   data = readtable('data.csv');
%   vary = 'Response';
%   fac1 = 'Factor1';
%   fac2 = 'Factor2';
%   rndfac = 'RandomFactor';
%   disptbl = 'on';
%   lsd_tbl = mlsdanova(data, vary, fac1, fac2, rndfac, disptbl);
%
%--------------------------------------------------------------------------

if nargin == 5
    disptbl = 'off';
end

[~, ~, stats] = anovan(data{:, vary}, {data{:, fac1}, data{:, fac2}, data{:, rndfac}}, ...
    'model', [1 0 0; 0 1 0; 0 0 1; 1 1 0], 'varnames', [fac1, fac2, rndfac], ...
    'nested', [0 0 0; 0 0 0; 0 1 0], 'random', 3, 'display', disptbl);

[~, i_md] = sort(data{:, fac2}, 'ascend');
data = data(i_md, :);
[u_meth_md, ~, i_meth_md] = unique(data{:,[fac1, fac2]}, 'rows', 'stable');
mlsd_groups = nan(length(u_meth_md), 2);

% Group means
for kg = 1:length(u_meth_md)
    mlsd_groups(kg, 1) = mean(data{i_meth_md == kg, vary}, 'omitnan');
end

nrep = length(unique(data{:,"Rep"}));

% LSD intervals according to ANOVA stored in stats
mlsd_groups(:, 2) = sqrt(stats.mse * 2 / nrep) * (tinv(1 - 0.025, stats.dfe));  

lsd_tbl = array2table(u_meth_md);
lsd_tbl.Properties.VariableNames = [fac1, fac2];
lsd_tbl.m = mlsd_groups(:, 1);
lsd_tbl.wlsd = mlsd_groups(:, 2);

end
