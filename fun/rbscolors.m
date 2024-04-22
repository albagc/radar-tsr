function colscale_red2blue = rbscolors(nlevels)
%--------------------------------------------------------------------------
% rbscolors - Generate a red-to-blue color scale
%--------------------------------------------------------------------------
%
% Inputs:
%   nlevels: Number of color levels in the scale.
%
% Output:
%   colscale_red2blue: Color scale matrix ranging from red to blue.
%
% Example:
%   nlevels = 10; % Number of color levels
%   colscale_red2blue = rbscolors(nlevels);
%
% Note: The output colscale_red2blue is a matrix with nlevels rows and 3 columns representing RGB values.
%--------------------------------------------------------------------------

colscale_red2white = [linspace(200, 255, nlevels)', linspace(0, 255, nlevels)', linspace(0, 255, nlevels)'];
colscale_white2blue = [linspace(255, 0, nlevels)', linspace(255, 0, nlevels)', linspace(255, 200, nlevels)'];
colscale_red2blue = flipud([colscale_red2white; colscale_white2blue] ./ 255);
end
