function dplotbygroup(T2, SPE, limt2, limspe, g, outnames, obstags, textsizelabs, diffltype)
%--------------------------------------------------------------------------
% dplotbygroup - Plot data points grouped by clusters
%--------------------------------------------------------------------------
%   dplotbygroup(T2, SPE, limt2, limspe, g, outnames, obstags, textsizelabs, diffltype)
%   plots the data points grouped by clusters and labeled accordingly.
%
%   Input:
%   - T2: T^2 values for each data point.
%   - SPE: SPE (Squared Prediction Error) values for each data point.
%   - limt2: Threshold for T^2 values for highlighting.
%   - limspe: Threshold for SPE values for highlighting.
%   - g: Cluster tags for each data point.
%   - outnames: Flag indicating whether to display observation tags.
%               Default is false.
%   - obstags: Observation tags or labels associated with each data point.
%              Default is 1:size(T2, 1).
%   - textsizelabs: Font size for observation tags. Default is 10.
%   - diffltype: Flag indicating whether to use different line types for
%                different clusters. Default is false.
%
%   Example:
%   T2 = [1.2, 0.8, 1.5, 0.9, 1.1];
%   SPE = [0.6, 0.9, 0.7, 1.0, 0.8];
%   limt2 = 1.0;
%   limspe = 0.8;
%   g = [1, 1, 2, 2, 3];
%   outnames = true;
%   obstags = ["A", "B", "C", "D", "E"];
%   textsizelabs = 12;
%   diffltype = true;
%   dplotbygroup(T2, SPE, limt2, limspe, g, outnames, obstags, textsizelabs, diffltype);
%
%--------------------------------------------------------------------------

arguments
    T2
    SPE
    limt2
    limspe
    g
    outnames logical = false
    obstags double = 1:size(T2, 1)
    textsizelabs double = 10
    diffltype = false
end

ug = unique(g);
m_rw = rw_detect(SPE, T2, limspe, limt2);
id_single_rw = and(m_rw, g == min(ug));
y_labels = string(g);
y_labels(id_single_rw) = "0 sRW";
y_labels(g > 0) = strcat("gRW - ", y_labels(g > 0));

u_labels = unique(y_labels);
col_scale = [0, 0, 0; 1, 0, 0; 0, 0, 1; 0, 185/256, 0; 1, 0.8275, 0; 1, 0.1020, 0.7255;...
    0, 0.1725, 0; 0.5176, 0.5176, 1];
col_scale = col_scale(1:length(ug), :);
if any(id_single_rw)
    y_colors = [0, 0, 0; col_scale];
else
    y_colors = col_scale;
end

if diffltype
    ltype = ["o", "^", "sq", "d", "+", "*", "x"]';
    lw = [0.5, 0.5, 1.2, 1.2, 0.5, 1.2];
    la = [0, 0, 0, 0, 1, 1];
    y_shapes = ltype(1:length(u_labels));
else
    y_shapes = ["o", repmat(["^"], 1, length(u_labels) - 1)];
    lw = zeros(length(y_shapes), 1);
    la = zeros(length(y_shapes), 1);
end

% Plot clusters
for jg = 1:length(u_labels)
    idobs = ismember(y_labels, u_labels(jg));
    scatter(T2(idobs), SPE(idobs), 'filled', ...
        "Marker", y_shapes(jg), "MarkerEdgeColor", y_colors(jg, :), ...
        "MarkerFaceColor", y_colors(jg, :), 'MarkerFaceAlpha', 0.5, ...
        'DisplayName', u_labels(jg), "linewidth", lw(jg), ...
        "MarkerEdgeAlpha", la(jg)), 
    
    hold on, grid on
end

xlabel('T^2'), ylabel('SPE'),
xline(limt2, 'r--', 'DisplayName', 'UCL_{SPE, T^2}', 'LineWidth', 1.2),
yline(limspe, 'r--', 'HandleVisibility', 'off', 'LineWidth', 1.2),
xline(2.5*limt2, 'r-', 'DisplayName', 'c_{SPE, T^2}', 'LineWidth', 1.2),
yline(2.5*limspe, 'r-', 'HandleVisibility', 'off', 'LineWidth', 1.2),
if outnames
    id_rw = find(m_rw);
    text(1.01*T2(id_rw), 1.01*SPE(id_rw), string(obstags(id_rw)), "FontSize", textsizelabs)
end

lg = legend('location', 'bestoutside', "Orientation", "horizontal");
title(lg, "Cluster tag")
end
