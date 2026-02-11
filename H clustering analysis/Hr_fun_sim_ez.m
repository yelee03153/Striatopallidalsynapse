function [H,H_upper, H_lower, H_all]=Hr_fun_sim_ez(loc,area,b,sims,ez_sims,poly,min_d,ce,ez_radius,num_ves)
% calculation of the H(r) function and CE for EZ model
%% input:
% - loc: coordinates of two points
% - area: AZ area 
% - b: 
% - sims: number of simulations 
% - ez_sims: number of EZ simulations 
% - poly: coordinates of the AZ 
% - min_d : minimum allowed distance between points
% - ce: confidence level 
% - ez_radius: radius of the EZ
% - num_ves: number of vesicles
%% output
% - H : H(r) function for the data
% - H_upper : upper CE
% - H_lower : lower CE
% - H_all : H(r) function for all the randomly generated samples
%%

[n,~] = size(loc);
H_all=[];
%% H(r) function for the data points
H = Hr(loc,b,area,poly(:,1),poly(:,2));
%% CE generation 
for j = 1:ez_sims   
    %generation of EZ model in the AZ
     [xu_c,yu_c,~,~]=EZ(poly(:,1),poly(:,2),area,ez_radius,num_ves);
     k=1;
     while k<100
        [xr,yr]=points2(poly(:,1),poly(:,2),min_d);
        [in,~]=inpolygon(xr,yr,poly(:,1),poly(:,2));
        xr_i=xr(in);
        yr_i=yr(in);
        [in1,~]=inpolygon(xr_i,yr_i,xu_c,yu_c);
        xr_i1=xr_i(~in1);
        yr_i1=yr_i(~in1);
        if length(yr_i1)>=n
            rand_p(1,:)=xr_i1(1:n);
            rand_p(2,:)=yr_i1(1:n);
             k=k+1;     
            % H(r) function for the random pattern
             H_all =[H_all Hr(rand_p,b,area,poly(:,1),poly(:,2))];
        end
     end
end   
%Build envelopes
H_rank=sort(H_all);
ce_n=round(sims*ez_sims/ce);
H_upper =H_rank(ce_n,:);
H_lower = H_rank(end-ce_n,:);
