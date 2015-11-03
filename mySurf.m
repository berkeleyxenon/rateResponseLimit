function h=mySurf(xdata,ydata,xbins,ybins);
%
%
% 120807 pfs - wrapper for surf

[xmesh,ymesh]=meshgrid(xbins(1:end-1),ybins(1:end-1));
clear ncts;
for xi=1:length(xbins)-1
for yi=1:length(ybins)-1
    ncts(yi,xi) = sum( (xdata>=xbins(xi)) & (xdata<xbins(xi+1)) & (ydata>=ybins(yi)) & (ydata<ybins(yi+1)));
	dis('xi=%d , yi=%d',xi,yi);
end
end


h=surf(xmesh,ymesh,-ncts); % by default, plot inverse, so I can plot lines over the 0,90 view
set(h,'EdgeColor','none'...
    ,'LineStyle','none' ...
    );
colormap('bone');

%shading interp;
view(0,90);

%    hsv        - Hue-saturation-value color map.
%    hot        - Black-red-yellow-white color map.
%    gray       - Linear gray-scale color map.
%    bone       - Gray-scale with tinge of blue color map.
%    copper     - Linear copper-tone color map.
%    pink       - Pastel shades of pink color map.
%    white      - All white color map.
%    flag       - Alternating red, white, blue, and black color map.
%    lines      - Color map with the line colors.
%    colorcube  - Enhanced color-cube color map.
%    vga        - Windows colormap for 16 colors.
%    jet        - Variant of HSV.
%    prism      - Prism color map.
%    cool       - Shades of cyan and magenta color map.
%    autumn     - Shades of red and yellow color map.
%    spring     - Shades of magenta and yellow color map.
%    winter     - Shades of blue and green color map.
%    summer     - Shades of green and yellow color map.
 
