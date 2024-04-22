function [r2o, yel, whi, p2b, cgrad] = rybscolors(nlevels)
%--------------------------------------------------------------------------
% rybscolors - Generate a colormap transitioning from red to yellow to blue
%--------------------------------------------------------------------------
%
% Inputs:
% nlevels: Number of levels in the colormap
%
% Outputs:
% r2o: Colormap transitioning from red to orange
% yel: Yellow color
% whi: White color
% p2b: Colormap transitioning from purple to blue
% cgrad: Complete colormap transitioning from red to yellow to blue
%
% Example:
% nlevels = 256;
% [r2o, yel, whi, p2b, cgrad] = rybscolors(nlevels);

r2o = [linspace(200,255,nlevels)',linspace(0,134,nlevels)',linspace(0,18,nlevels)']./256;
% colscale_yellow2whyellow = [linspace(255,255,ceil(nlevels/4))',linspace(255,255,ceil(nlevels/4))',linspace(25,150,ceil(nlevels/4))'];
yel = [255, 223, 0]./256;
whi = [1, 1, 1];
p2b = [linspace(209,0,nlevels)',linspace(156,0,nlevels)',linspace(255,200,nlevels)']./256;
% Intermediates
nlevhalf = round(nlevels/4);
o2y = [linspace(255,yel(1),nlevhalf)',linspace(134,yel(2),nlevhalf)',linspace(18,yel(3),nlevhalf)']./256;
y2w = [linspace(yel(1),256,nlevhalf)',linspace(yel(2),256,nlevhalf)',linspace(yel(3),256,nlevhalf)']./256;
w2y = flipud(y2w);
y2p = [linspace(yel(1),209,nlevhalf)',linspace(yel(2),156,nlevhalf)',linspace(yel(3),255,nlevhalf)']./256;
cgrad = flipud([r2o;o2y;y2w;w2y;y2p;p2b]);

end