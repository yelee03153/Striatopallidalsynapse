%% create new pattern by dividing the given one in two
function [H,dens]=Hr_half_group(folder_name,dist)

%% input: 
%folder_name: path the folder where patterns on the group of interest
%are located
% dist: diastenced used to compute H(r) functions 

d = dir(folder_name);
d(1:2)=[];
ncols = length(d);
num_=0;

for i=1:ncols
A = dlmread(d(i).name);
    x=A(:,1)';
    y=A(:,2)';
 
    %draw an artifical minimum contour for the pattern, defined be max and
    %min coordinates of the pattern
     xp(1)=min(x)+1;
     xp(2)=min(x)+(max(x)-min(x))/2;
     xp(3)=min(y)+1;
     xp(4)=max(y);    

    are(i)=(xp(2)-xp(1))*(xp(4)-xp(3)); % area of the box
    ind_h=x<(xp(2));
    ind_l=x>(xp(2));

    num(i)=length(x(ind_h)); % nu,ber of points in the pattern

    if num(i)>10 % consider patterns that havle only more than 10 points
       num_=num_+1;
       dens_(num_)=num(i)./are(i);   
    %% Calculate H(r) for data, CE for randomly sampled patterns
       [H_1(num_,:)]=Hr([x(ind_h) ;y(ind_h)],dist,xp);
       [H_2(num_,:)]=Hr([x(ind_l) ;y(ind_l)],dist,xp);
    end
    
end
H=[H_1; H_2];
dens=[dens_, dens_];