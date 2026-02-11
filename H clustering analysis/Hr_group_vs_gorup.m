% Univariate H(r) anlaysis for group vs group

dist=1:1:4000; % spatial scale to calculate H(r) (in pixel units)
h_pattern_di_h=[];
h_pattern_mi_h=[];
pval_pattern_di_h=[];
pval_pattern_mi_h=[];


[H_group_a,dens_a]=Hr_group('C:\Users\USER\Downloads\Univariate_Analysis_JH\Univariate_Analysis_JH\Larger_Scale\OTu',dist);
% in case you want to increase number of patterns by dividing each in two use this 
%[H_group_a,dens_a]=Hr_half_group('/Users/maria/Documents/IP/Final_IDN/codes/Univariate_Analysis_JH/test/A',dist);

[H_group_b,dens_b]=Hr_group('C:\Users\USER\Downloads\Univariate_Analysis_JH\Univariate_Analysis_JH\Larger_Scale\NAs',dist);

null_group= H_group_a; %group that will a null model
null_dens=dens_a;

test_group=H_group_b; %group that will be compared to then null model
test_dens=dens_b;
[n_H,~]=size(test_group);
[n_Hall,~]=size(null_group);
% loop trough all the H(r) of the test group
for i=1:n_H
    
    %% plot individual functions
      
    %% MAD & DCLF tests for the current pattern
    [h_pattern_di_h(i,:),pval_pattern_di_h(i,:)]=dclf_test(null_group,n_Hall,test_group(i,:),95, dist);

    [h_pattern_mi_h(i,:),pval_pattern_mi_h(i,:)]=mad_test1(null_group,n_Hall,test_group(i,:),95);
    
    end

%%  % of patterns that differ from the null_droup
%for DCLF test 
su=sum(h_pattern_di_h');
res_dclf=sum(su>0)/n_H*100

%for MAD test
res_mad=sum(h_pattern_mi_h)/n_H*100

%pooled H(r)
 H_ncols=pool_it(test_group,n_H,test_dens);
 H_rank=sort(null_group);
%  ce_n=n_Hall;
 
ce_n=ceil(n_Hall/95);
H_upper = H_rank(ce_n,:);
H_lower = H_rank(end-ce_n,:);
H_ncols_null_group=pool_it(null_group,n_Hall,null_dens);
 
%  MAD test for the population
 [H_mad,pval_mad]=mad_test1(null_group,n_Hall,H_ncols,95);
% DCLF test for the population
 [H_dclf,pval_dclf]=dclf_test(null_group,n_Hall,H_ncols,95,dist);


      hold on,plot(dist,H_upper,'--b')
      hold on,plot(dist,H_lower,'--b')
     for i=1:n_H
    hold on, plot(dist,test_group(i,:),'r' )
     end
     hold on, plot(dist,H_ncols_null_group,'b' )





