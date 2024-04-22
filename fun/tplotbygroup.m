function tplotbygroup(T, t1, t2, limt1, limt2, g, outnames, T2, SPE, limT2, limspe, obstags, textsizelabs, diffltype)
%--------------------------------------------------------------------------
% tplotbygroup - Scatter plot of data points colored by group membership
%--------------------------------------------------------------------------
%
% Inputs:
%   T: Data matrix
%   t1: Index of the first variable for the scatter plot
%   t2: Index of the second variable for the scatter plot
%   limt1: Limits for the first variable axis
%   limt2: Limits for the second variable axis
%   g: Group membership vector
%   outnames: Flag to display observation names (default: false)
%   T2: Second data matrix (default: empty)
%   SPE: SPE statistic vector (default: empty)
%   limT2: Limits for the second data matrix axis (default: empty)
%   limspe: Limit for the SPE statistic axis (default: empty)
%   obstags: Observation tags (default: 1:size(T,1))
%   textsizelabs: Font size for observation tags (default: 10)
%   diffltype: Flag to differentiate line types (default: false)
%
% Example:
%   T = rand(100, 3); % Data matrix
%   t1 = 1; % Index of the first variable
%   t2 = 2; % Index of the second variable
%   limt1 = [0, 1]; % Limits for the first variable axis
%   limt2 = [0, 1]; % Limits for the second variable axis
%   g = randi(3, 100, 1); % Group membership vector
%   outnames = true; % Display observation names
%   tplotbygroup(T, t1, t2, limt1, limt2, g, outnames); % Scatter plot with group coloring and observation names

arguments
    T
    t1
    t2
    limt1
    limt2
    g
    outnames logical = false
    T2 double = []
    SPE double = []
    limT2 double = []
    limspe double = []
    obstags double = 1:size(T,1)
    textsizelabs double = 10
    diffltype = false
end

ug = unique(g);
m_rw = rw_detect(SPE, T2, limspe, limT2);
id_single_rw = and(m_rw, g == min(ug));
y_labels = string(g);
y_labels(id_single_rw) = "0 sRW";
y_labels(g > 0) = strcat("gRW - ", y_labels(g > 0));

u_labels = unique(y_labels);
col_scale = [0, 0, 0; 1, 0, 0; 0, 0, 1; 0, 185/256, 0; 1, 0.8275, 0; 1, 0.1020, 0.7255;...
    0, 0.1725, 0; 0.5176, 0.5176, 1];
col_scale = col_scale(1:length(ug), :);

% Adjust colors and shapes based on single rowwise outliers
if any(id_single_rw)
    y_colors = [0, 0, 0; col_scale];
else
    y_colors = col_scale;
end

% Determine line types and widths based on different group types
if diffltype
    ltype = ["o", "^", "sq", "d", "+", "*", "x"]';
    lw = [0.5, 0.5, 1.2, 1.2, 0.5, 1.2];
    la = [0, 0, 0, 0, 1, 1];
    y_shapes = ltype(1:length(u_labels));
else
    y_shapes = ["o", repmat(["^"], 1, length(u_labels) - 1)];
    lw = ones(length(y_shapes), 1);
    la = zeros(length(y_shapes), 1);
end

[xyellips] = doellips(limt1, limt2, 100, 0, 0);
plot(xyellips(:, 1), xyellips(:, 2), 'r--', 'HandleVisibility', 'off', 'LineWidth', 1.2), 
hold on, grid on

% Plot data points colored by group membership
for jg = 1:length(u_labels)
    idobs = ismember(y_labels, u_labels(jg));
    scatter(T(idobs, t1), T(idobs, t2), 'filled', ...
        "Marker", y_shapes(jg), "MarkerEdgeColor", y_colors(jg, :), ...
        "MarkerFaceColor", y_colors(jg, :), 'MarkerFaceAlpha', 0.5, ...
        'DisplayName', u_labels(jg), "linewidth", lw(jg), ...
        "MarkerEdgeAlpha", la(jg)), 
end

xlabel(strcat("t_{", string(t1), "}")),
ylabel(strcat("t_{", string(t2), "}")),
title('Data Cortex Nuclear'),
xline(0, 'k', 'HandleVisibility', 'off')
yline(0, 'k', 'HandleVisibility', 'off')

lg = legend('location', 'bestoutside', "Orientation", "horizontal");
title(lg, "Cluster tag")

% Display observation names if specified
if outnames
    id_rw = find(m_rw);
    text(1.01 * T(id_rw, t1), 1.01 * T(id_rw, t2), string(obstags(id_rw)), ...
        'FontSize', textsizelabs)
end

end
