function [H,H_upper, H_lower, H_all]=Hr_fun_sim(loc,area,b,sims,poly,min_d,ce)
% calculation of the H(r) function and CE for null model
%% input:
% - loc: coordinates of two points
% - area: AZ area 
% - b: 
% - sims: number of simulations 
% - poly: coordinates of the AZ 
% - min_d : minimum allowed distance between points
% - ce: confidence level 
%% output
% - H : H(r) function for the data
% - H_upper : upper CE
% - H_lower : lower CE
% - H_all : H(r) function for all the randomly generated samples
%%

[n,~] = size(loc');
%% H(r) function for the data points
H = Hr(loc,b,poly);
%% CE generation 
for i = 1:sims   
    %generation of random points in the AZ
    [xr,yr]=points2(poly,min_d);
    %[in,~]=inpolygon(xr,yr,poly(:,1),poly(:,2));
    %x_temp=xr(in);
    %y_temp=yr(in);
    rand_p(1,:)=xr(1:n);
    rand_p(2,:)=yr(1:n);   
    % H(r) function for the random pattern
    H_all(i,:) = Hr(rand_p,b,poly);   
end   
%Build envelopes
H_rank=sort(H_all);
ce_n=round(sims/ce);
H_upper =H_rank(ce_n,:);
H_lower = H_rank(end-ce_n,:);
