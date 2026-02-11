%mad test
function [H0,p, del_t_hi, T_obs]=dclf_test(arr,n,data,ce,r)
% pefrorm DCLF test
%% input 
% - arr: H(r) function for the randomly generated patterns
% - n: number of simulations 
% - data: H(r) for the data pattern 
% - ce : CE level 
% - r : the distaces at which DCLF is performed
%% output: H0: null hypothesis, p- p value 

iso_mean=mean(arr); % null model 
H0=0; % intialization of the null hypothesis 

% calculae differnces between H(r) of each pattern of the null model and
% the  null model
for i=1:max(r)
    for j=1:n
    T_sim_temp(i,j)= mean((arr(j,r(1:i))-iso_mean(r(1:i))).^2);
    end
end

% take the differnce between H(r) for observed pattern and null model for
% each distacne in r
for i=1:max(r)
    T_obs(i)= mean((data(r(1:i))-iso_mean(r(1:i))).^2);
end
l=ceil(ce*n/100); % take the index that corresponds to 95% of the number of patterns in the null model 

% for each distance in r perpform test 
for i=1:max(r)
    del_t_sim(i,:)=sort(T_sim_temp(i,:)); % sort the differnces
    del_t_hi(i)=del_t_sim(i,l); % take lth biggest difference (corresponds to the 95% of all difference)
    if (T_obs(i)> del_t_hi(i)) 
        H0(i)=1;
    else
        H0(i)=0;
    end
    num_t(i)=sum(del_t_sim(i,:)>T_obs(i));
    p(i)=(num_t(i))/(n);
end

