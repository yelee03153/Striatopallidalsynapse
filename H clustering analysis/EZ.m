% Maria Reva 2016.05

function [xu_c,yu_c,keeperX1,keeperY1]=EZ(xp,yp,area,ez,ves_num)
% generation of the EZ 
%% input
% - xp,yp: coordinates of the AZ polygon
% - area: area of the EZ
% - ez: size of the EZ (nm)
%% outup
% - xu_c,yu_c: coordinates of the EZ
% - keeperX1,keeperY1: coordinates of the EZ centers

% generate points inside the smallest shape aroun the EZ
x=min(xp)+(max(xp)-min(xp))*rand(1,1000);
y=min(yp)+(max(yp)-min(yp))*rand(1,1000);

% for the circular AZ
% rad =175; %radius of the circulaer AZ
%  theta = rand(1,5000)*(2*pi);
%  r = rad*sqrt(rand(1,500)); %
%  x = rad+ r.*cos(theta);
%  y = rad+ r.*sin(theta);

[in,~] = inpolygon(x,y,xp,yp); % keep only those that are inside the AZ
x=x(in);
y=y(in);

minAllowableDistance =40; % min distance between vesicle centers

% Initialize first EZ central point
keeperX1 = x(1);
keeperY1 = y(1);
% dropping down more points
counter = 2;
for k = 2 : length(x)
	% Get a trial point
	thisX = x(k);
	thisY = y(k);
	% See how far is is away from existing keeper points
	distances = sqrt((thisX-keeperX1).^2 + (thisY - keeperY1).^2);
	minDistance = min(distances);
	if (minDistance >= minAllowableDistance) 
		keeperX1(counter) = thisX;
		keeperY1(counter) = thisY;
		counter = counter + 1;
	end
end

% draw vesciels around these poits
th=0:0.01:2*pi;
for i =1 : length(keeperX1)
    xunit(i,:)=20* cos(th) + keeperX1(i);
    yunit(i,:)=20* sin(th) + keeperY1(i);
end

% draw EZ arounf these points 
xunit1=zeros(20, length(th));
yunit1=zeros(20, length(th));
for i =1 : length(keeperX1)
    xunit1(i,:)=ez* cos(th) + keeperX1(i);
    yunit1(i,:)=ez* sin(th) + keeperY1(i);
end

x_ez={};
y_ez={};
x_ez{1}=xunit1(1,:);
y_ez{1}=yunit1(1,:);

% create intersections of the EZs
for j=2:20
    [x_ez{j},y_ez{j}]=polybool('union', (x_ez{1,j-1}),(y_ez{1,j-1}), xunit1(j,:),yunit1(j,:));
end 

% keep number of vesicles that is proportionla to the size of the AZ
n_ves=round(ves_num*area/9e4);
if n_ves>1
     ind=min(20,n_ves + round(-1+2*rand));
else
    ind= max(1,n_ves+round(rand));
end

xu_c=x_ez{1,ind};
yu_c=y_ez{1,ind};

end
