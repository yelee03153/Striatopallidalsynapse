function [w]=edge_corr(DIST,xp,yp,loc_x,loc_y)
% calculate edge correction
%%input
% - DIST: distances between points
% - xp, yp: coordinates of the window (eg AZ)
% - loc_x, loc_y : locations of the points
%% output: w- edge correciton factor


[n,m]=size(DIST);

for i=1:n
   for j=1:m
      w(i,j)=frac(xp,yp,DIST(i,j),loc_x(i),loc_y(i));   
   end
end