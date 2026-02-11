close all
clear all

d = dir('C:\Users\USER\Downloads\Univariate_Analysis_JH\Univariate_Analysis_JH\Larger_Scale\OTu');
d(1:2)=[];
ncols = length(d);

for i=1:ncols
     A = dlmread(d(i).name);
    x=A(:,1);
    y=A(:,2);
    
    
    %draw an artifical minimum contour for the pattern, defined be max and
    %min coordinates of the pattern
    xp(1)=min(x)+1;
    xp(2)=max(x)+1;
    xp(3)=min(y)+1;
    xp(4)=max(y)+1;    
    
    are(i)=(xp(2)-xp(1))*(xp(4)-xp(3)); % area of the box
    num_ca(i)=length(x); %number of "cav" points
    min_d=1;%minimum allowed distance between points
    
    %% generate null model
   for j=1:500
         [xr,yr]=points2(xp,min_d);
          xr_i_le=xr(1:num_ca(i));
          yr_i_le=yr(1:num_ca(i));
         % find nnd of the randomchannels
          NND_rand=nnds(num_ca(i),xr_i_le,yr_i_le);
          Mean_NND_rand(i,j)=mean(NND_rand);
   end    
    
    %% generate EZ model
%   num_ves=8;
%   ez_rad=50;
%    for j=1:100
%         [xu_c,yu_c,keeperX1,keeperY1]=EZ(xp,yp,area,ez_rad,num_ves);
%         k=1;
%         while k<100
%               [xr,yr]=points2(xp,yp,min_d);
%               [in,on]=inpolygon(xr,yr,xp,yp);
%                xr_i=xr(in);
%                yr_i=yr(in);
%                [in1,on1]=inpolygon(xr_i,yr_i,xu_c,yu_c);
%                xr_i1=xr_i(~in1);
%                yr_i1=yr_i(~in1);
%                if length(yr_i1)>=num_ca(i)
%                    xr_i_le=xr_i1(1:num_ca(i));
%                    yr_i_le=yr_i1(1:num_ca(i));
%                    k=k+1;
%                    % find nnd of the randomchannels
%                    NND_ez=nnds(num_ca(i),xr_i_le,yr_i_le);
%                    temp_m(j,k)=mean(NND_ez);
%                end
%         end
%    Mean_NND_EZ(i,j)=mean(temp_m(j,:));
%    end

    % mean NND for the data
    NND_exp=nnds(num_ca(i),x',y');
    Mean_NND_ex(i)=mean(NND_exp);
end 
% MW test
[h_mean_NND,p_mean_NND]=ranksum(Mean_NND_ex,mean(Mean_NND_rand'))
% plotting
figure, cdfplot(Mean_NND_ex)
hold on, cdfplot(mean(Mean_NND_rand'))
