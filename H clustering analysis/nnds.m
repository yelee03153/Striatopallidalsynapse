function [nnd_le]= nnds(l,x,y)
% find nnd
%%input
% - l: number of the poitns
% - x,y: coordinates
%% output- nnd_le: nnd vector

for i=1:l
   [k,DIX]=knnsearch([x(i); y(i)]', [x; y]', 'k',1);
   DIX=sort(DIX);
   nnd_le(i)=DIX(2);
    
end