function errormap(x, xscale, xname, ethresh, rnames, ylabels, vnames, xlabels)
%--------------------------------------------------------------------------
% errormap - Plot error map of a data matrix
%--------------------------------------------------------------------------
%   errormap(x, xscale, xname, ethresh, rnames, ylabels, vnames, xlabels)
%   plots the error map of a data matrix x, highlighting the magnitude of
%   the errors.
%
%   Input:
%   - x: Data matrix.
%   - xscale: Flag indicating whether to perform autoscaling of the data.
%             Default is false.
%   - xname: Name of the data.
%   - ethresh: Error threshold for coloring. Values exceeding ethresh will
%              be capped. Default is 0.
%   - rnames: Flag indicating whether to display row labels. Default is false.
%   - ylabels: Labels for the y-axis (observations). Default is string(1:size(x, 1)).
%   - vnames: Flag indicating whether to display column labels. Default is false.
%   - xlabels: Labels for the x-axis (variables). Default is string(1:size(x, 2)).
%
%   Example:
%   x = [1, 2, 3; 4, 5, 6];
%   xscale = true;
%   xname = "Data X";
%   ethresh = 1;
%   rnames = true;
%   ylabels = ["Obs1", "Obs2"];
%   vnames = true;
%   xlabels = ["Var1", "Var2", "Var3"];
%   errormap(x, xscale, xname, ethresh, rnames, ylabels, vnames, xlabels);
%
%--------------------------------------------------------------------------

arguments
   x
   xscale
   xname
   ethresh
   rnames = false
   ylabels = string(1:size(x, 1));
   vnames = false
   xlabels = string(1:size(x, 2));
end

if xscale
    x = (x - robloc(x)) ./ repmat(robscale(x), size(x, 1), 1);
    ttext = strcat(xname, " (autoscaled)");
    x_sc = x;
    x(x_sc > ethresh) = ethresh;
    x(x_sc < -ethresh) = -ethresh;
    M = ethresh;
else
    ttext = xname;
    M = max(abs(x(:)));
end

rwb_colscale = rbscolors(10);

tiledlayout(1, 10 * round(log(size(x, 2))/log(10)))
nexttile([1, 10 * round(log(size(x, 2))/log(10)) - 1])
imagesc(x),
colormap(rwb_colscale), clim([-M, M]),
if ~rnames
    yticks([1, size(x, 1)]), yticklabels(["1", "n"])
else
    yticks(1:size(x, 1)), yticklabels(ylabels)
end
if ~vnames
    xticks([1, size(x, 2)]), xticklabels(["1", "d"]), xtickangle(90)
else
    xticks(1:size(x, 2)), xticklabels(xlabels), xtickangle(90)
end
title(ttext), xlabel('Variables'), ylabel('Observations')
a = get(gca, 'XTicklabel');
set(gca, 'XTickLabel', a, 'fontsize', 6)
ax = nexttile([1, 1]);
imagesc(sum(x.^2, 2)), colormap(ax, flipud(gray))
yticks([1, size(x, 1)]), yticklabels(["1", "n"])
title('SE')
a = get(gca, 'XTicklabel');
set(gca, 'XTickLabel', a, 'fontsize', 6)

end
