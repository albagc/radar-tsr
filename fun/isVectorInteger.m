function result = isVectorInteger(x, tol)
    if nargin < 2
        tol = sqrt(eps);
    end

    y = abs(x - round(x)) < tol;
    result = sum(1 - y) == 0;
end
