clear all
close all

% Univariate H(r) anlaysis
d = dir('C:\Users\USER\Downloads\Univariate_Analysis_JH\Univariate_Analysis_JH\Bassoon Gephyrin\Bassoon Gephyrin ACC');

%initialization 
sims=100; % number of simulations
dist=0:1:50; % spatial scale to calculate H(r) (in pixel units)
min_d=1;

H_all=[];
d(1:2)=[];
ncols = length(d);

figure, 
% loop trough all the files in the folder
for i=1:ncols
    %read coordinates from ascii file, no normalization
    A = dlmread(d(i).name);
    x=A(:,1)';
    y=A(:,2)';
    
    
    %draw an artifical minimum contour for the pattern, defined be max and
    %min coordinates of the pattern
    xp(1)=min(x)+1;
    xp(2)=max(x)+1;
    xp(3)=min(y)+1;
    xp(4)=max(y)+1;    
    
    are(i)=831744; % area of the box
    
 
    num_ca(i)=length(x); %number of points
    %% Calculate H(r) for data, CE for randomly sampled patterns
    [H_data(i,:), H_upper(i,:), H_lower(i,:), H_rand_temp]=Hr_fun_sim([x ;y],are(i),dist,sims,xp,min_d,95);
    %% Calculate H(r) for data, CE for EZ model
    %ez_sim=100; % number of EZ patterns generated
    %ez_rad=50; % radius of EZ
    %num_ves=4; % number of vesicels
    %[H_data(i,:), H_upper(i,:), H_lower(i,:), H_rand_temp]=Hr_fun_sim_ez([ca_x ca_y],are(i),dist,sims,ez_sim,[xp yp],min_d,99,ez_rad,num_ves);
    %% plot individual functions
    hold on, plot(dist,H_data(i,:),'b')
    hold on, plot(dist,H_upper(i,:),'--b' )
    hold on, plot(dist,H_lower(i,:),'--b' )

    % store all the H(r) generated for CE
    H_all=[H_all; H_rand_temp];
    % MAD test for the current pattern
    [h_pattern(i,:),pval_pattern(i)]=mad_test1(H_rand_temp,sims,H_data(i,:),95);
end

%pooled H(r)
dens=num_ca./are;
[n_H,~]=size(H_data);
[n_Hall,~]=size(H_all);
H_ncols=pool_it(H_data,n_H,dens);
H_rank=sort(H_all);
ce_n=round(n_Hall/95);
H_upper = H_rank(ce_n,:);
H_lower = H_rank(end-ce_n,:);
% MAD test for the population
[H_tot,pval_tot]=mad_test1(H_all,ncols*sims,H_ncols,95);
% plot
figure, plot(dist, H_ncols,'r')
hold on, plot(dist, H_upper,'b')
hold on, plot(dist, H_lower,'b')
hold on, plot(dist, mean(H_all),'g')
title('Pooled')






