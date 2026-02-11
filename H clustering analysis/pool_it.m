function  [pooled,varr,y_all,X]=pool_it(arr,n,weights)
%pool
    all_pf=[];
    y_all=[];
    sum_dens=sum(weights(1:n));
    
   for i=1:n
        all_pf(i,:)=arr(i,:).*weights(i);
    end 
    y_all=sum(all_pf);
    pooled=y_all./(sum_dens);
    
    %var
    m=length(sum_dens);
    x_s=sum(y_all);
    y_s=sum(weights);
%    construct a temp vec X
    for i =1:m
       X1(i)=n*y_all(i);
       Y1(i)=n*sum_dens(i);
    end
   


   
