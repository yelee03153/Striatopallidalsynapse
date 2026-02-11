function frac=frac(xp,yp,d,loc_x,loc_y)
% calculate fraction of the circle with center at loc_x,loc_y that is inside of polygon
%% input
% - xp,yp: coordinates of the polygon 
% - d: radius of the circle 
% - loc_x,loc_y: center of the circle
%% output: frac -fraction of the circle that is inside of polygon 
%Requires the arclength function by John D'Errico
%   http://www.mathworks.com/matlabcentral/fileexchange/34871-arclength
%%
% set coordinates of the circle
dt=d-.05;
th=0:0.0001:2*pi;
xc = dt.*cos(th)+loc_x; 
yc = dt.*sin(th)+loc_y;

in1 =find(inpolygon(xc, yc,xp, yp)); % find intersenction of the polygon and circle

if isempty(in1)
in1=[1 2]; % if it's empty fill it with the single point
end

if d>0 
xi=xc(in1); % take the coordinates of the circel that are inside the polygon 
yi=yc(in1);
arclen = arclength(xi,yi); % calculate its length 
frac= (2*pi*d)/arclen; % take fracntion 
else
    frac=0;
end

