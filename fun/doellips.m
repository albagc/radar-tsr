function [xyellips, xy_int] = doellips(xr, yr, n, xc, yc)
%--------------------------------------------------------------------------
% doellips - Generate ellipse coordinates
%--------------------------------------------------------------------------
%   [xyellips, xy_int] = doellips(xr, yr, n, xc, yc) generates coordinates
%   for an ellipse based on the given parameters.
%
%   Input:
%   - xr: X-axis radius of the ellipse.
%   - yr: Y-axis radius of the ellipse.
%   - n: Number of intervals for the intermediate points of the ellipse.
%        Default is 1.
%   - xc: X-coordinate of the center of the ellipse. Default is 0.
%   - yc: Y-coordinate of the center of the ellipse. Default is 0.
%
%   Output:
%   - xyellips: Coordinates of points along the ellipse boundary.
%   - xy_int: Intermediate coordinates within the ellipse.
%
%   Example:
%   xr = 5; % X-axis radius
%   yr = 3; % Y-axis radius
%   n = 2;  % Number of intervals
%   xc = 1; % X-coordinate of center
%   yc = 2; % Y-coordinate of center
%   [xyellips, xy_int] = doellips(xr, yr, n, xc, yc);
%
%--------------------------------------------------------------------------

arguments
    xr
    yr
    n = 1
    xc = 0
    yc = 0
end

theta = 0 : 0.01 : 2*pi;
x = xr * cos(theta) + xc;
y = yr * sin(theta) + yc;
xyellips = [x', y'];

theta_int = [linspace(0, pi/2, n+2), linspace(pi/2, pi, n+2)];
x_int = xr * cos(theta_int) + xc;
y_int = yr * sin(theta_int) + yc;
xy_int = [x_int', y_int'];
xy_int = [xy_int; -xy_int];

end
