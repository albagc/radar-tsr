function metrics = evalresults(yobs, ypred, mode)
%--------------------------------------------------------------------------
% evalresults - Evaluate classification or prediction results
%--------------------------------------------------------------------------
%   metrics = evalresults(yobs, ypred, mode) evaluates the results of
%   classification or prediction and returns a table of evaluation metrics.
%
%   Input:
%   - yobs: Observed values.
%   - ypred: Predicted values.
%   - mode: Mode of evaluation. 'class' for classification, 'multiclass'
%           for multiclass classification, 'pred' for prediction.
%
%   Output:
%   - metrics: Table of evaluation metrics.
%
%   Example:
%   yobs = [1, 0, 1, 1];
%   ypred = [1, 1, 0, 1];
%   mode = 'class';
%   metrics = evalresults(yobs, ypred, mode);
%
%--------------------------------------------------------------------------

arguments
    yobs {logical, double}
    ypred {logical, double}
    mode char
end

switch mode
    case 'class'
        if islogical(yobs)
            yobs = double(yobs);
        end
        if islogical(ypred)
            ypred = double(ypred);
        end
        % True Positives
        tp = sum(yobs .* ypred);
        % True negatives
        tn = sum((yobs == 0) .* (ypred == 0));
        % False Positives
        fp = sum((yobs == 0) .* (ypred == 1));
        % False Negatives
        fn = sum((yobs == 1) .* (ypred == 0));
        % Accuracy (only recommended with balanced classes)
        acc = (tp + tn) / (tp + tn + fp + fn);
        % Precision (should we trust a positive result?)
        prec = tp / (tp + fp);
        % Recall (sensitivity/true positive rate)
        recall = tp / (tp + fn);
        % Specifitity (1 - false positive rate)
        fpr = fp / (fp + tn);
        spec = 1 - fpr;
        
        metrics = struct2table(struct('TruePositives', tp, ...
            'TrueNegatives', tn, 'FalsePositives', fp, ...
            'FalseNegatives', fn, 'Accuracy', acc, ...
            'Precision', prec, ...
            'Recall', recall, 'Specificity', spec, ...
            'FalsePositiveRate', fpr));
        
    case 'multiclass'
        yobs_dum = dummyvar(categorical(yobs));
        ypred_dum = dummyvar(categorical(ypred));
        cat_value = categories(categorical(yobs));
        metrics = array2table(nan(length(cat_value), 9));
        metrics.Properties.VariableNames = cellstr({'TruePositives', ...
            'TrueNegatives', 'FalsePositives', 'FalseNegatives', 'Accuracy', ...
            'Precision', 'Recall', 'Specificity', 'FalsePositiveRate'});
        
        for j = 1:size(yobs_dum, 2)
            % True Positives
            tp = sum(yobs_dum(:, j) .* ypred_dum(:, j));
            % True negatives
            tn = sum((yobs_dum(:, j) == 0) .* (ypred_dum(:, j) == 0));
            % False Positives
            fp = sum((yobs_dum(:, j) == 0) .* (ypred_dum(:, j) == 1));
            % False Negatives
            fn = sum((yobs_dum(:, j) == 1) .* (ypred_dum(:, j) == 0));
            % Accuracy (only recommended with balanced classes)
            acc = (tp + tn) / (tp + tn + fp + fn);
            % Precision (should we trust a positive result?)
            prec = tp / (tp + fp);
            % Recall (sensitivity/true positive rate)
            recall = tp / (tp + fn);
            % Specifitity (1 - false positive rate)
            fpr = fp / (fp + tn);
            spec = 1 - fpr;
            
            metrics.TruePositives(j) = tp;
            metrics.TrueNegatives(j) = tn;
            metrics.FalsePositives(j) = fp;
            metrics.FalseNegatives(j) = fn;
            metrics.Precision(j) = prec;
            metrics.Recall(j) = recall;
            metrics.Specificity(j) = spec;
            metrics.FalsePositiveRate(j) = fpr;
        end      
        
        metrics.Categories = cat_value;
        
    case 'pred'
        % MSE
        mse = mean((yobs - ypred).^2);
        % MAE
        mae = median(abs(yobs - ypred));
        % Rel error
        [minerr, imin] = min(abs(yobs - ypred));
        [maxerr, imax] = max(abs(yobs - ypred));
        error_range = [minerr, maxerr];
        relerror_range = [minerr / yobs(imin), maxerr / yobs(imax)];
        
        metrics = struct2table(struct('MSE', mse, ...
            'MAE', mae, 'ErrorRange', error_range, ...
            'RelErrorRange', relerror_range));
end

end
