function [errparams] = errormap_rowavge_new(x, standod, cod_factor, xscale, estnorm, xname, ndiv, ethresh, nrowsblock, rnames, ylabels, vnames, xlabels)
%--------------------------------------------------------------------------
% errormap_rowavge_new - Plot error map with row averaging
%--------------------------------------------------------------------------
%   [errparams] = errormap_rowavge_new(x, xscale, estnorm, xname, ndiv, ethresh, nrowsblock, rnames, ylabels, vnames, xlabels)
%   plots the error map of a data matrix x with row averaging, highlighting
%   the magnitude of the errors.
%
%   Input:
%   - x: Data matrix.
%   - xscale: Flag indicating whether to perform autoscaling of the data.
%             Default is false.
%   - estnorm: Estimation norm used for scaling. Default is not specified.
%   - xname: Name of the data.
%   - ndiv: Number of divisions for tiled layout.
%   - ethresh: Error threshold for coloring. Values exceeding ethresh will
%              be capped. Default is 0.
%   - nrowsblock: Number of rows per block for row averaging. Can be either
%                 a scalar or a 2-element vector. Default is not specified.
%   - rnames: Flag indicating whether to display row labels. Default is false.
%   - ylabels: Labels for the y-axis (observations). Default is string(1:size(x, 1)).
%   - vnames: Flag indicating whether to display column labels. Default is false.
%   - xlabels: Labels for the x-axis (variables). Default is string(1:size(x, 2)).
%
%   Output:
%   - errparams: Error parameters struct, containing m (mean) and s (scale)
%                if xscale is true, and Esc (error-scaled data matrix).
%
%   Example:
%   x = [1, 2, NaN; NaN, 5, 6];
%   xscale = true;
%   estnorm = 2;
%   xname = "Data X";
%   ndiv = 2;
%   ethresh = 1;
%   nrowsblock = [1, 2];
%   rnames = true;
%   ylabels = ["Obs1", "Obs2"];
%   vnames = true;
%   xlabels = ["Var1", "Var2", "Var3"];
%   errormap_rowavge_new(x, xscale, estnorm, xname, ndiv, ethresh, nrowsblock, rnames, ylabels, vnames, xlabels);
%
%--------------------------------------------------------------------------

arguments
    x
    standod
    cod_factor
    xscale
    estnorm
    xname
    ndiv
    ethresh
    nrowsblock
    rnames = false
    ylabels = string(1:size(x,1));
    vnames = false
    xlabels = string(1:size(x,2));
end

if xscale
    mE = locmask(x, isnan(x), estnorm);
    sE = scalemask(x, isnan(x), estnorm);
    x = (x - mE) ./ repmat(sE, size(x, 1), 1);
    errparams.m = mE;
    errparams.s = mE;
    errparams.xscale = xscale;
end
errparams.Esc = x;

ttext = strcat(xname);

M = 3;
x1dec = round(x, 1);
ncols = length(round(ethresh, 1):0.1:round(M, 1));
[r2o, yel, whi, p2b] = rybscolors(ncols);
n_levin = length(round(-ethresh):0.1:round(ethresh));

xind = x1dec;
xind(x > M) = M;
xind(x < -M) = M;
xind(isnan(x)) = 0;
xind(abs(x) < ethresh) = x(abs(x) < ethresh);

oy = mkgrad(r2o(end, :), yel, round(n_levin / 4));
yw = flipud(mkgrad(whi, yel, round(n_levin / 2)));
yw = yw(2:round(0.8 * n_levin / 2), :);
yp = mkgrad(yel, p2b(1, :), round(n_levin / 4));
oywyp = [oy; yw(2:end-1, :); flipud(yw); yp];

cgrad = flipud([r2o; oywyp; p2b]);
SE_all = standod;%sum(xind.^2, 2, 'omitnan');

x(isnan(x)) = 0;
if length(nrowsblock) == 1
    n_blocks = round(size(x, 1) / nrowsblock);
    m_block_x = nan(n_blocks, size(x, 2));
    m_block_x_RGB = nan(n_blocks, size(x, 2), 3);
    for jb = 1:n_blocks
        m_block_x(jb, :) = mean(x(nrowsblock*(jb-1)+1:nrowsblock*jb, :));
        m_block_x_RGB(jb, :, 1) = mean(x(nrowsblock*(jb-1)+1:nrowsblock*jb, :, 1));
        m_block_x_RGB(jb, :, 2) = mean(x(nrowsblock*(jb-1)+1:nrowsblock*jb, :, 2));
        m_block_x_RGB(jb, :, 3) = mean(x(nrowsblock*(jb-1)+1:nrowsblock*jb, :, 3));
    end
elseif length(nrowsblock) == 2
    m_block_x = blockproc(xind, [nrowsblock(1), nrowsblock(2)], @(x)mean(x.data, 'all'));
    se_block_x = blockproc(SE_all, [nrowsblock(1), nrowsblock(2)], @(x)mean(x.data, 'all'));
end

% Generate layout of plots
tiledlayout(1, ndiv)
% Error map
nexttile([1, ndiv-1])
% Plot squares 
imagesc(m_block_x), hold on
colormap(cgrad), clim([-M, M]),
% Plot lines between squares
for xcellmap = 1:size(m_block_x, 2)
    line([xcellmap-0.5, xcellmap-0.5], [0.5, size(m_block_x, 1)+0.5], 'Color', 'w', 'HandleVisibility', 'off')
end
for ycellmap = 1:size(m_block_x, 1)
    line([0.5, size(m_block_x, 2)+0.5], [ycellmap-0.5, ycellmap-0.5], 'Color', 'w', 'HandleVisibility', 'off')
end
% Observation names as y-axis ticks
if ~rnames
    yticks([1, size(m_block_x, 1)]), yticklabels(["1", "n"])
else
    yticks(1:size(m_block_x, 1)), yticklabels(ylabels)
end
% Variable names as x-axis ticks
if ~vnames
    xticks([1, size(m_block_x, 2)]), xticklabels(["1", "d"]), xtickangle(90)
else
    xticks(1:size(m_block_x, 2)), xticklabels(xlabels), xtickangle(90)
end
% Error map title
title(ttext), xlabel('Variables'), ylabel('Observations')

% Distance colorbar
%cod_factor;
limL = 0; %1/cod_factor;
limH = cod_factor; %1;
se_block_x(se_block_x > cod_factor) = limH;
stand_se_block_x = ((se_block_x - limL)/(limH - limL));
ax = nexttile([1, 1]);
imagesc(stand_se_block_x), colormap(ax, flipud(gray))
yticks([1, size(m_block_x, 1)]), yticklabels(["1", "n"]), clim([0, 1])
title('SE')
end
