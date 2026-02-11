% random point generation 
function [keepX,keepY]=points2(xp, min_d)
%%input
% - xp: box
% - min_d:minimum allowed distance between the points 
%%output: keepX,keepY: coordinates of randomly generated points


%%for circular area
%  theta = rand(1,5000)*(2*pi);
%  r = sqrt(rand(1,5000))*175;
%  x =  r.*cos(theta);
%  y = r.*sin(theta);

%% for the box
x=xp(1)+(xp(2)-xp(1))*rand(1,10000);
y=xp(3)+(xp(4)-xp(3))*rand(1,10000);


keepX = x(1);
keepY = y(1);

counter = 2;
for k = 2 :10000
thisX = x(k);
	thisY = y(k);
	% See how far is is away from existing keeper points
	distances = sqrt((thisX-keepX).^2 + (thisY - keepY).^2);
	minDistance = min(distances);
	if (minDistance >=min_d) % check in points are not further than min_d
		keepX(counter) = thisX;
		keepY(counter) = thisY;
		counter = counter + 1;
	end
end
end
