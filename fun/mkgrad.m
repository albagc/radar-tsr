function cgrad = mkgrad(clow, cup, nlevs)
%--------------------------------------------------------------------------
% mkgrad - Create a color gradient
%--------------------------------------------------------------------------
%   cgrad = mkgrad(clow, cup, nlevs) creates a color gradient from the
%   starting color clow to the ending color cup with a specified number of
%   levels nlevs.
%
%   Input:
%   - clow: Starting color specified as an RGB triplet [R, G, B].
%   - cup: Ending color specified as an RGB triplet [R, G, B].
%   - nlevs: Number of levels in the color gradient.
%
%   Output:
%   - cgrad: Color gradient matrix with size [nlevs, 3].
%
%   Example:
%   clow = [0, 0, 0];
%   cup = [1, 1, 1];
%   nlevs = 10;
%   cgrad = mkgrad(clow, cup, nlevs);
%
%--------------------------------------------------------------------------

cgrad = [linspace(clow(1), cup(1), nlevs)', ...
         linspace(clow(2), cup(2), nlevs)', ...
         linspace(clow(3), cup(3), nlevs)'];

end
