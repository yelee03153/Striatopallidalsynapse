%mad test
function [H0,p]=mad_test1(arr,n,data,ce)
% pefrorm MAD test
%% input 
% - arr: H(r) function for the randomly generated patterns
% - n: number of simulations 
% - data: H(r) for the data pattern 
% - ce : CE level 
%% output: H0: null hypothesis, p- p value 

iso_mean=mean(arr);
H0=0;

% calculae differnces between H(r) for the data pattern and its mean

for i=1:n
    delta_t_sim(i)= max(abs((arr(i,:)-iso_mean)));
end

del_t_data=max(abs(data-iso_mean)); % diffrence between H(r) data and mean of H(r) for simulated patterns 
del_t_sim=sort(delta_t_sim);
l=ceil(ce*(n)/100); 
del_t_hi=del_t_sim(l); % take the difference value that corresponds to the CE level
if (del_t_data> del_t_hi) % MAD test
    H0=1;
end
num_t=sum(del_t_sim>del_t_data);
% p value 
p=(1+num_t)/(n+1);