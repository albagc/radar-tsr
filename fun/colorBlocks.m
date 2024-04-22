function [Xblock, Xblockgrad] = colorBlocks(Xin, rowblocksizes, columnblocksizes, colContrast)
    n = length(rowblocksizes);
    d = length(columnblocksizes);
    Xblock = zeros(n, d);
    Xblockgrad = zeros(n, d);
    rind = cumsum([0, rowblocksizes]);
    cind = cumsum([0, columnblocksizes]);
    
    for i = 1:n
        for j = 1:d
            Xsel = Xin((rind(i) + 1):rind(i + 1), (cind(j) + 1):cind(j + 1));
            seltable = histcounts(Xsel, 1:5);
            
            if sum(seltable) > 0
                [cntmax, indmax] = max(seltable);
                ncells = rowblocksizes(i) * columnblocksizes(j);
                gradmax = (cntmax / ncells).^ (1 / colContrast);
            else
                indmax = 0;
                gradmax = 1;
            end
            
            Xblock(i, j) = indmax;
            Xblockgrad(i, j) = gradmax;
        end
    end
    
    Xblock = Xblock; % Transpose to match R's column-major format
    Xblockgrad = Xblockgrad; % Transpose to match R's column-major format
    
    % Return as cell array (equivalent to R list)
    % Xblock = {Xblock, Xblockgrad};
end
