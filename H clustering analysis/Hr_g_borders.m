function H = Hr_g_borders(locs,dist,area,xp,yp )
% Calculate H(r) function for the point pattern 
%%input 
% - locs: location of the points
% - dist: range of distances for calculaing H(r) (r=dist)
% - area: area of the polygon (eg AZ)
% - xp,yp: coordinates of the boundaries 

%% output: H normilized Ripley function H(r)

%initialization
locs=locs';
[N,~] = size(locs);
K = zeros(length(dist),1);

%distances among the points
an=pdist(locs);
DIST=squareform(an);
K = zeros(length(dist),1);
%calcualte weigths 
w=edge_corr(DIST,xp,yp,locs(:,1),locs(:,2));
for k=1:length(K)
    b=(DIST(1:end,:)<dist(k));
    K(k) = sum(sum(w(b)))/N;
end
%Normalization
lambda = N/area;
K = K/lambda;
%H = diff(K')./(2*pi*pi*dist(1:end-1));
%H=sqrt(K/pi);
H=sqrt(K/pi)-dist';