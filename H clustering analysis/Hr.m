function H = Hr_g(locs,dist,box)
% Calculate H(r) function for the point pattern 
%%input 
% - locs: location of the points
% - dist: range of distances for calculaing H(r) (r=dist)
% - area: area of the polygon (eg AZ)
% - DIST: distances between the points 
% - w: edge correction factor
%% output: H normilized Ripley function H(r)

%initialization
locs=locs';
[N,~] = size(locs);
dbox = min([locs(:,1)'-box(1);box(2)-locs(:,1)';locs(:,2)'-box(3); box(4)-locs(:,2)'] ); 
%distance to a box
%distances among the points
diX = repmat(locs(:,1),1,N)-repmat(locs(:,1)',N,1);
diY = repmat(locs(:,2),1,N)-repmat(locs(:,2)',N,1);
DIST = sqrt(diX.^2+diY.^2);
DIST = sort(DIST);
K = zeros(length(dist),1);
for k=1:length(K)
    K(k) = sum(sum(DIST(2:end,:)<dist(k)))/N;
end
%Normalization
lambda = N/((box(2)-box(1))*(box(4)-box(3)));
K = K/lambda;
%H = diff(K')./(2*pi*pi*dist(1:end-1));
%H=sqrt(K/pi);
H=sqrt(K/pi)-dist';