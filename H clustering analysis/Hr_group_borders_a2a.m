%% generate H(r) functions for particular model 
function [H,dens]=Hr_group_borders(folder_name,dist,folder_name_borders)

%% input: 
%folder_name: path the folder where patterns on the group of interest
%are located
% dist: diastenced used to compute H(r) functions 

d = dir(folder_name);

d_borders=dir(folder_name_borders);

d(1:2)=[];
d_borders(1:2)=[];
ncols = length(d);
num_=0;

for i=1:ncols
    cd(folder_name);
    A = dlmread(d(i).name);
    x=A(:,1)';
    y=A(:,2)';
 
    % read the borders 
    cd(folder_name_borders);
    Borders = dlmread(d_borders(i).name);
    xp=Borders(:,1)';
    yp=Borders(:,2)';

    num(i)=length(x); % nu,ber of points in the pattern
    are(i)=polyarea(xp,yp); % area of the box
    
    if num(i)>3 % consider patterns that havle only more than 10 points
       num_=num_+1;
       dens(num_)=num(i)./are(i);   
    %% Calculate H(r) for data, CE for randomly sampled patterns
       [H(num_,:)]=Hr_g_borders([x ;y],dist,are(i),xp,yp);
    end
    clear A
    clear borders
end