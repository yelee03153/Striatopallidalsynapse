%% generate H(r) functions for particular model 
function [H,dens]=Hr_group(folder_name,dist)

%% input: 
%folder_name: path the folder where patterns on the group of interest
%are located
% dist: diastenced used to compute H(r) functions 

d = dir(folder_name);
d(1:2)=[];
ncols = length(d);
num_=0;

cd(folder_name);

for i=1:ncols
A = dlmread(d(i).name);
    x=A(:,1)';
    y=A(:,2)';
 
    %draw an artifical minimum contour for the pattern, defined be max and
    %min coordinates of the pattern
    xp(1)=min(x)+1;
    xp(2)=max(x)+1;
    xp(3)=min(y)+1;
    xp(4)=max(y)+1;  

    num(i)=length(x); % nu,ber of points in the pattern
    are(i)=(xp(2)-xp(1))*(xp(4)-xp(3)); % area of the box
    
    if num(i)>10 % consider patterns that havle only more than 10 points
       num_=num_+1;
       dens(num_)=num(i)./are(i);   
    %% Calculate H(r) for data, CE for randomly sampled patterns
       [H(num_,:)]=Hr([x ;y],dist,xp);
    end
end