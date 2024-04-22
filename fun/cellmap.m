function cellmapimg = cellmap(R, ucl_od, indcells, indrows, standOD, showVals, D, ...
    rowlabels, columnlabels, mTitle, rowtitle, columntitle, ndiv, showrows, ...
    showcolumns, nrowsinblock, ncolumnsinblock, manualrowblocksizes, ...
    manualcolumnblocksizes, rowblocklabels, columnblocklabels, autolabel, ...
	columnangle, colContrast, outlyingGrad, darkestColor, drawCircles) 
    arguments
        R
        ucl_od
        indcells = []
        indrows = []
        standOD = []
        showVals = []
        D = []
        rowlabels = []
        columnlabels = []
        mTitle = "cell map"
        rowtitle = "cases"
        columntitle = "variables"
        ndiv = 15
        showrows = []
        showcolumns = []
        nrowsinblock = []
        ncolumnsinblock = []
        manualrowblocksizes = []
        manualcolumnblocksizes = []
        rowblocklabels = []
        columnblocklabels = []
        autolabel = true
        columnangle = 90
        colContrast = 1
        outlyingGrad = true
        darkestColor = ucl_od
        drawCircles = true
    end
    numColors = 20;
    type = 'cell';
    n = size(R, 1);
    d = size(R, 2);
    
    if ~isempty(showVals)
        if ~any(strcmp(showVals, {'D', 'R'}))
            error('Invalid "showVals" argument. Should be one of: [], "D", "R"');
        end
    end
    
    if isempty(showVals)
        D = R;
    elseif strcmp(showVals, 'D') && isempty(D)
        error('When showVals="D", you must input the argument D');
    end
    
    if isempty(D)
        D = R;
    elseif ~isequal(size(D), size(R))
        error('The dimensions of D and R must match');
    end
    
    if isempty(indcells)
        indcells = find(abs(R) > norminv(1 - 0.01/2));
    end
    
    if ~isempty(rowlabels)
        if length(rowlabels) ~= n
            error(['Number of rowlabels does not match n = ', num2str(n)]);
        end
    else
        if istable(R) 
            if isempty(R.Properties.RowNames)
                rowlabels = 1:n;
            else
                rowlabels = R.Properties.RowNames;
            end
        else 
            rowlabels = 1:n;
        end
    end
    
    if ~isempty(columnlabels)
        if length(columnlabels) ~= d
            error(['Number of columnlabels does not match d = ', num2str(d)]);
        end
    else
        if istable(R) 
            if isempty(R.Properties.VariableNames)
                columnlabels = 1:d;
            else
                columnlabels = R.Properties.VariableNames;
            end
        else 
            columnlabels = 1:d;
        end
    end
    
    if ~isempty(showcolumns) || ~isempty(showrows)
        if isempty(showrows)
            showrows = 1:n;
        else
            if ~all(ismember(showrows, 1:n))
                error('showrows goes out of bounds');
            end
        end
    
        if isempty(showcolumns)
            showcolumns = 1:d;
        else
            if ~all(ismember(showcolumns, 1:d))
                error('showcolumns goes out of bounds');
            end
        end
    
        tempMat = zeros(n, d);
        tempMat(indcells) = 1;
        tempMat = tempMat(showrows, showcolumns);
        indcells = find(tempMat == 1);
    
        tempVec = zeros(1, n);
        tempVec(indrows) = 1;
        tempVec = tempVec(showrows);
        indrows = find(tempVec == 1);
        clear tempMat tempVec;
    
        R = R(showrows, showcolumns);
        D = D(showrows, showcolumns);
        rowlabels = rowlabels(showrows);
        columnlabels = columnlabels(showcolumns);
    
        n = size(R, 1);
        d = size(R, 2);
    
        if ~isempty(standOD)
            standOD = standOD(showrows);
        end
    end
    
    blockRows = false;
    blockColumns = false;
    
    if ~isempty(manualrowblocksizes)
        if ~isempty(nrowsinblock)
            disp('Input argument manualrowblocksizes has overruled argument nrowsinblock.');
            blockRows = true;
        end
    
        if ~isvector(manualrowblocksizes)
            error('manualrowblocksizes should be a vector with strictly positive integers, adding up to at most n');
        end
    
        if ~isnumeric(manualrowblocksizes)
            error('manualrowblocksizes should be a vector with strictly positive integers, adding up to at most n');
        end
    
        if sum(manualrowblocksizes < 1) > 0
            error('manualrowblocksizes should be a vector with strictly positive integers, adding up to at most n');
        end
    
        if sum(manualrowblocksizes) > n
            error('manualrowblocksizes should be a vector with strictly positive integers, adding up to at most n');
        end
    
        if sum(manualrowblocksizes ~= 1) == 0
            error('All manualrowblocksizes are 1');
        end
    
        blockRows = true;
    else
        if ~isempty(nrowsinblock)
            if nrowsinblock > 1
                if nrowsinblock > n
                    error(['Input argument nrowsinblock cannot be more than n = ', num2str(n)]);
                end
                blockRows = true;
            end
        end
    end
    
    if ~isempty(manualcolumnblocksizes)
        if ~isempty(ncolumnsinblock)
            disp('Input argument manualcolumnblocksizes has overruled argument ncolumnsinblock.');
        end
    
        if ~isvector(manualcolumnblocksizes)
            error('manualcolumnblocksizes should be a vector with strictly positive integers, adding up to at most d');
        end
    
        if ~isnumeric(manualcolumnblocksizes)
            error('manualcolumnblocksizes should be a vector with strictly positive integers, adding up to at most d');
        end
    
        if sum(manualcolumnblocksizes < 1) > 0
            error('manualcolumnblocksizes should be a vector with strictly positive integers, adding up to at most d');
        end
    
        if sum(manualcolumnblocksizes) > d
            error('manualcolumnblocksizes should be a vector with strictly positive integers, adding up to at most d');
        end
    
        if sum(manualcolumnblocksizes ~= 1) == 0
            error('All manualcolumnblocksizes are 1');
        end
    
        blockColumns = true;
    else
        if ~isempty(ncolumnsinblock)
            if ncolumnsinblock > 1
                if ncolumnsinblock > d
                    error(['Input argument ncolumnsinblock cannot be more than d = ', num2str(d)]);
                end
                blockColumns = true;
            end
        end
    end
    
    
    if (blockRows || blockColumns) && ~isempty(showVals)
        warning('The option showVals="D" or showVals="R" cannot be combined with blocking rows and/or columns, so showVals is set to NULL here.');
        showVals = [];
    end
    
    if strcmp(type, 'residual')
        outlyingGrad = 1;
    end
    
    X = zeros(n, d);
    Xrow = zeros(n, 1);
    Xrow(indrows) = 3;
    
    if strcmp(type, 'cell') || blockRows || blockColumns
        indcells_mat = zeros(size(R));
        indcells_mat(indcells) = 1;
        pcells = find(indcells_mat > 0 & R >= 0);
        ncells = find(indcells_mat > 0 & R < 0);
    else
        pcells = find(R >= 0);
        ncells = find(R < 0);
    end
    
    X(ncells) = 1;
    X(pcells) = 2;
    X(isnan(R)) = 4;
    
    if (blockRows || blockColumns)
        rowblocksizes = ones(1, n);
        
        if blockRows
            if ~isempty(manualrowblocksizes)
                rowblocksizes = manualrowblocksizes;
                n = length(rowblocksizes);
            elseif ~isempty(nrowsinblock) && nrowsinblock > 1
                n = floor(n / nrowsinblock);
                rowblocksizes = repmat(nrowsinblock, 1, n);
            end
        end
    
        columnblocksizes = ones(1, d);
        
        if blockColumns
            if ~isempty(manualcolumnblocksizes)
                columnblocksizes = manualcolumnblocksizes;
                d = length(columnblocksizes);
            elseif ~isempty(ncolumnsinblock) && ncolumnsinblock > 1
                d = floor(d / ncolumnsinblock);
                columnblocksizes = repmat(ncolumnsinblock, 1, d);
            end
        end
    
        [Xcol, Xgrad]  = colorBlocks(X, rowblocksizes, columnblocksizes, colContrast);

        if drawCircles
            [Xcol_sd, Xgrad_sd] = colorBlocks(Xrow, rowblocksizes, ones(1, 1), colContrast);
            Xrowgrad = Xgrad_sd;
            Xrowgrad(Xcol_sd == 0) = 0;
        end
    
        if blockRows
            if isempty(rowblocklabels)
                fprintf('No rowblocklabels were given, so they are constructed automatically.\n');
                laby = rowlabels;
                rowlabels = zeros(1, n);
                rind = cumsum([0 rowblocksizes]);
                
                for i = 1:n
                    if rowblocksizes(i) == 1
                        rowlabels(i) = laby(rind(i) + 1);
                    else
                        rowlabels(i) = string(strcat(laby(rind(i) + 1), '-', laby(rind(i + 1))));
                    end
                end
            else
                if length(rowblocklabels) ~= n
                    error(['The number of rowblocklabels is ', num2str(length(rowblocklabels)), ...
                        ' but there are ', num2str(n), ' row blocks.']);
                end
                rowlabels = rowblocklabels;
            end
        end
    
        if blockColumns
            if isempty(columnblocklabels)
                fprintf('No columnblocklabels were given, so they are constructed automatically.\n');
                labx = columnlabels;
                columnlabels = zeros(1, d);
                cind = cumsum([0 columnblocksizes]);
                
                for j = 1:d
                    if columnblocksizes(j) == 1
                        columnlabels(j) = labx(cind(j) + 1);
                    else
                        columnlabels(j) = string(strcat(labx(cind(j) + 1), '-', labx(cind(j + 1))));
                    end
                end
            else
                if length(columnblocklabels) ~= d
                    error(['The number of columnblocklabels is ', num2str(length(columnblocklabels)), ...
                        ' but there are ', num2str(d), ' column blocks.']);
                end
                columnlabels = columnblocklabels;
            end
        end
        
        Xdf = array2table([transpose(1:n), Xcol]);
        Xdf.Properties.VariableNames = [{'rownr'}, string(1:d)];
        Xdf.Properties.RowNames = {};
        
        [~, new_order] = ismember(categorical(Xdf.rownr), categorical(1:n));
        Xdf.rownr = Xdf.rownr(new_order);
        mX = table(reshape(Xdf{:,2:end}, numel(Xdf{:,2:end}), 1));
        mX.Properties.VariableNames = "CatNr";
        mX.rownr = repmat(Xdf.rownr, size(Xdf{:,2:end}, 2),1);

        Xgraddf = array2table([transpose(1:n), Xgrad]);
        Xgraddf.Properties.VariableNames = [{'rownr'}, string(1:d)];
        Xgraddf.Properties.RowNames = {};
        
        [~, new_order_grad] = ismember(categorical(Xgraddf.rownr), categorical(1:n));
        Xgraddf.rownr = Xgraddf.rownr(new_order_grad);
        mXgrad = table(reshape(Xgraddf{:,2:end}, numel(Xgraddf{:,2:end}), 1));
        mXgrad.Properties.VariableNames = "grad";
        mXgrad.rownr = repmat(Xgraddf.rownr, size(Xgraddf{:,2:end}, 2),1);

        
        mX.grad = mXgrad.grad;
        mX.rescaleoffset = mXgrad.grad + 10 * mX.CatNr;
    
        if drawCircles
            mXrow = table(transpose(1:size(Xrowgrad, 1)), Xrowgrad + 10 * 3, ...
                'VariableNames', {'rownr', 'rescaleoffset'});
        end
    
        [r2o, yel, whi, p2b] = rybscolors(numColors);
        mX.contscale = mX.grad;
        mX.contscale = mX.grad + mX.CatNr - 1;
        uColors = unique(mX.CatNr);
        rneg_rpos = ismember([1, 2], uColors);
        if strcmp(type, 'cell')
            if all(rneg_rpos)
                colorends = [mkgrad(yel, [0.8164, 0.6094, 0.9961], 10);...
                                mkgrad([0.8164, 0.6094, 0.9961], [0, 0, 0.7812], 10);...
                                flipud(mkgrad([0.9961, 0.5234, 0.0703], yel, 10));...
                                flipud(mkgrad([0.7812, 0, 0], [0.9961, 0.5234, 0.0703], 10))];
            elseif and(rneg_rpos(1), ~rneg_rpos(2))
                colorends = [mkgrad(yel, [0.8164, 0.6094, 0.9961], 10);...
                                mkgrad([0.8164, 0.6094, 0.9961], [0, 0, 0.7812], 10)];
            elseif and(rneg_rpos(2), ~rneg_rpos(1))
                colorends = [flipud(mkgrad([0.9961, 0.5234, 0.0703], yel, 10));...
                                flipud(mkgrad([0.7812, 0, 0], [0.9961, 0.5234, 0.0703], 10))];
            end
            if ismember(3,uColors)
                colorends = [colorends; flipud(gray(20))];
            elseif ismember(4,uColors)
                colorends = [colorends; flipud(gray(20)); mkgrad(yel, whi, 20)];
            end
        elseif strcmp(type, 'residual')
            if all(rneg_rpos)
                colorends = [mkgrad(whi, [0.7812, 0, 0], 20);...
                    mkgrad(whi, [0, 0, 0.7812], 20)];
            elseif and(rneg_rpos(1), ~rneg_rpos(2))
                colorends = [mkgrad(whi, [0.7812, 0, 0], 20)];
            elseif and(rneg_rpos(2), ~rneg_rpos(1))
                colorends = [mkgrad(whi, [0, 0, 0.7812], 20)];
            end
            if ismember(3,uColors)
                colorends = [colorends; flipud(gray(20))];
            elseif ismember(4,uColors)
                colorends = [colorends; flipud(gray(20)); mkgrad(whi, whi, 20)];
            end
        end
    else
    
        Ddf = table((1:n)', D);
        Ddf.Properties.VariableNames = {'rownr', 'CatNr', 'Value'};
        Ddf.Properties.RowNames = {};
        
        mD = table(Ddf.rownr, Ddf.Value);
        mD.Properties.VariableNames = {'rownr', 'data'};
        mD.Properties.RowNames = {};
        
        Rdf = table((1:n)', R);
        Rdf.Properties.VariableNames = {'rownr', 'CatNr', 'Value'};
        Rdf.Properties.RowNames = {};
        
        mR = table(Rdf.rownr, Rdf.Value);
        mR.Properties.VariableNames = {'rownr', 'data'};
        mR.Properties.RowNames = {};
        
        Xdf = table((1:n)', X);
        Xdf.Properties.VariableNames = {'rownr', 'CatNr', 'Value'};
        Xdf.Properties.RowNames = {};
        
        mX = table(Xdf.rownr, Xdf.Value);
        mX.Properties.VariableNames = {'rownr', 'data'};
        mX.Properties.RowNames = {};
        
        if ~isempty(showVals)
            if strcmp(showVals, 'D')
                mX.data = mD.data;
            elseif strcmp(showVals, 'R')
                mX.data = mR.data;
            end
        end
        
        if ~outlyingGrad
            mX.rescaleoffset = 10 * mX.CatNr;
            if strcmp(type, 'cell')
                colorends = [linspace(0, 1, numColors); linspace(0, 1, numColors); linspace(1, 0, numColors)]'; 
            elseif strcmp(type, 'residual')
                colorends = [linspace(0, 1, numColors); linspace(0, 1, numColors); linspace(1, 1, numColors)]'; 
            end
        else
            Xgrad = nan(size(R));
            if strcmp(type, 'cell')
                Xgrad(indcells) = abs(R(indcells));
                limL = norminv(1 - 0.01/2);
            else
                Xgrad = abs(R);
                limL = 0;
            end
        
            limH =  norminv(1 - 0.01/2);
            Xgrad(Xgrad > limH) = limH;
            Xgrad = ((Xgrad - limL) / (limH - limL)).^colContrast;
            Xgrad(isnan(Xgrad)) = 0;
        
            Xgraddf = table((1:n)', Xgrad);
            Xgraddf.Properties.VariableNames = {'rownr', 'grad'};
            Xgraddf.Properties.RowNames = {};
        
            mXgrad = table(Xgraddf.rownr, Xgraddf.grad);
            mXgrad.Properties.VariableNames = {'rownr', 'grad'};
            mXgrad.Properties.RowNames = {};
        
            mX.grad = mXgrad.grad;
            mX.rescaleoffset = mXgrad.grad + 10 * mX.CatNr;
            [r2o, yel, whi, p2b] = rybscolors(ncols);
            mX.contscale = mX.grad;
            mX.contscale = mX.grad + mX.CatNr - 1;
            uColors = unique(mX.CatNr);
            rneg_rpos = ismember([1, 2], uColors);
            if strcmp(type, 'cell')
                if all(rneg_rpos)
                    colorends = [mkgrad(yel, [0.8164, 0.6094, 0.9961], 10);...
                                    mkgrad([0.8164, 0.6094, 0.9961], [0, 0, 0.7812], 10);...
                                    flipud(mkgrad([0.9961, 0.5234, 0.0703], yel, 10));...
                                    flipud(mkgrad([0.7812, 0, 0], [0.9961, 0.5234, 0.0703], 10))];
                elseif and(rneg_rpos(1), ~rneg_rpos(2))
                    colorends = [mkgrad(yel, [0.8164, 0.6094, 0.9961], 10);...
                                    mkgrad([0.8164, 0.6094, 0.9961], [0, 0, 0.7812], 10)];
                elseif and(rneg_rpos(2), ~rneg_rpos(1))
                    colorends = [flipud(mkgrad([0.9961, 0.5234, 0.0703], yel, 10));...
                                    flipud(mkgrad([0.7812, 0, 0], [0.9961, 0.5234, 0.0703], 10))];
                end
                if ismember(3,uColors)
                    colorends = [colorends; flipud(gray(20))];
                elseif ismember(4,uColors)
                    colorends = [colorends; flipud(gray(20)); mkgrad(yel, whi, 20)];
                end
            elseif strcmp(type, 'residual')
                if all(rneg_rpos)
                    colorends = [mkgrad(whi, [0.7812, 0, 0], 20);...
                        mkgrad(whi, [0, 0, 0.7812], 20)];
                elseif and(rneg_rpos(1), ~rneg_rpos(2))
                    colorends = [mkgrad(whi, [0.7812, 0, 0], 20)];
                elseif and(rneg_rpos(2), ~rneg_rpos(1))
                    colorends = [mkgrad(whi, [0, 0, 0.7812], 20)];
                end
                if ismember(3,uColors)
                    colorends = [colorends; flipud(gray(20))];
                elseif ismember(4,uColors)
                    colorends = [colorends; flipud(gray(20)); mkgrad(whi, whi, 20)];
                end
            end
        end
    end

    % Reverse the row labels
    rowlabels = flip(rowlabels);
    
    % Create the base size
    base_size = 10;
    limsupcolor = max(uColors);
    % Create the MATLAB figure
    tiledlayout(1, ndiv);
    % Error map
    nexttile([1, ndiv-1])
    imagesc(reshape(mX.contscale, size(Xgraddf, 1), size(Xgraddf, 2)-1)), colormap(colorends), clim([0,limsupcolor])
    % Plot lines between squares
    for xcellmap = 1:size(Xgraddf, 2)-1
        line([xcellmap-0.5, xcellmap-0.5], [0.5, size(Xgraddf, 1)+0.5], 'Color', 'w', 'HandleVisibility', 'off')
    end
    for ycellmap = 1:size(Xgraddf, 1)
        line([0.5, size(Xgraddf, 2)+0.5], [ycellmap-0.5, ycellmap-0.5], 'Color', 'w', 'HandleVisibility', 'off')
    end
    % Observation names as y-axis ticks
    if isempty(rowlabels)
        yticks([1, size(Xgraddf, 1)]), yticklabels(["1", "n"])
    else
        yticks(1:size(Xgraddf, 1)), yticklabels(rowlabels)
    end
    % Variable names as x-axis ticks
    if isempty(columnlabels)
        xticks([1, size(Xgraddf, 2)-1]), xticklabels(["1", "d"]), xtickangle(90)
    else
        xticks(1:size(Xgraddf, 2)-1), xticklabels(columnlabels), xtickangle(90)
    end
    % Error map title
    title(mTitle), xlabel(rowtitle), ylabel(columntitle)
    %%
     % Set X-axis properties
    xticks(1:size(mX, 2));
    xticklabels(columnlabels);
    xlabel(columntitle);
    
    % Set Y-axis properties
    yticks(1:size(mX, 1));
    yticklabels(rowlabels);
    ylabel(rowtitle);
    
    % Set title and axis labels
    title(mTitle);
    
    % Remove axis ticks and lines
    box off;
    ax0 = gca;
    ax0.XAxis.TickValues = [];
    ax0.YAxis.TickValues = [];

    % Add polygons if drawCircles is true
    if drawCircles
        hold on;
        % Error map
        ax = nexttile([1, 1]);
        imagesc(mXrow.rescaleoffset), colormap(ax, flipud(gray(20))), clim(ax, [30, 31])
        for ycellmap = 1:size(mXrow, 1)
            line([0.5, size(mXrow, 2)+0.5], [ycellmap-0.5, ycellmap-0.5], 'Color', 'w', 'HandleVisibility', 'off')
        end
        hold off;
    end

    % Remove axis ticks and lines
    box off;
    ax.XAxis.TickValues = [];
    ax.YAxis.TickValues = [];

    
    % Add text labels if showVals is not empty
    if ~isempty(showVals)
        txtcol = mX;
        txtcol(txtcol == 0) = "black";
        txtcol(txtcol == 1) = "white";
        txtcol(txtcol == 2) = "white";
        if strcmp(type, 'residual')
            txtcol = "black";
            txtcol(mXgrad > 0.5) = "white";
        end
        txtcol(txtcol == 4) = "black";
        [rows, cols] = size(mX);
        for row = 1:rows
            for col = 1:cols
                val = mX(row, col);
                if isnan(val)
                    txt = '';
                else
                    txt = sprintf('%1.0f', val);
                end
                text(col, row, txt, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', txtcol(row, col), 'FontSize', base_size * 0.5);
            end
        end
    end
        
       cellmapimg = gcf;
end


    


