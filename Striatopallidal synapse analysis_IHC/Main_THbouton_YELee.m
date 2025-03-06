
clear all
close all
%% Main
% part1: extracting pre signal
% part2: extracting post signal
% part3: extracting synapse
% part4: analysis 
%% initialize 
% pre Vmat2 488
params.deno = 0.1;
params.hard_lev1 = 40000;
params.otsu = 50; 
params.deno_px_number = 20;
params.dist = 2;

% body TH 647
params.deno2 = 0.1;
params.hard_lev2 = 32000;
params.otsu2 = 50; 
params.deno_px_number2 = 20;

% post RFP 594
params.deno3 = 0.1;
params.hard_lev3 = 25000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

 for v = 1
%% Data loading - pre vmat2

for ci = 1:3
    data1(:,:,ci) = imread('THbouton_DAPIVmat2RFPTH_63X3_L1_LD_Out',ci);
end

% data_pre = data(1:300,1:300,3:3:end);
% data_pre = data1(1:300,1:300,3:3:end);
data_pre = data1(:,:,2);

%input_pre = data_pre(1:500,1:500,13);
input_pre = data_pre(:,:,v);

% input = input(351:450,231:330);




%% Data loading -body TH

for ci = 1:3
    data2(:,:,ci) = imread('THbouton_DAPIVmat2RFPTH_63X3_L1_LD_Out',ci);
end

data_body = data2(:,:,3);
%  data_body = data2(1:300,1:300,1:3:end);

%input_body = data_body(1:500,1:500,13);
input_body = data_body(:,:,v);

% input = input(351:450,231:330);
% input = imread('cell_sample01.tif');
%% Data loading - post RFP
% input = imread('cell_sample01.tif');

for ci = 1:3
    data(:,:,ci) = imread('THbouton_DAPIVmat2RFPTH_63X3_L1_LD_Out',ci);
end

% data_post = data(1:300,1:300,2:3:end);
data_post = data(:,:,1);


%input_post = data_post(1:500,1:500,13);
input_post = data_post(:,:,v);
% input = input(351:450,231:330);
%% part1
% function - pre clustering= TH, Vmat2 colocal


res_clustering_pre = extpre(input_pre,input_body ,params);
%% part1-2
% function - TH area

THarea = extTH (input_body,params);

THnum = nnz(THarea);

%% part2
% function - post clustering= RFP area

res_clustering_post = extpostyelee (input_post,params);

RFPnum = nnz(res_clustering_post);
%% part3 
% function -  synapse extracting= TH, Vmat2, RFP colocal

d =1;
[res_synapse] = extsynapseyelee (res_clustering_pre,res_clustering_post);

% temp3 = tril(dmat_prepost);
% temp4 = sort(temp3(:));
% res_sort2_prepost = temp4(temp4~=0);
% figure, histogram(res_sort2_prepost)
% xlabel('distance')
% ylabel('N')
% title('pre2post')
% 
% if v==1
%     res_sort_prepost1 = res_sort2_prepost;
%     dmat_prepost1 = dmat_prepost;
% end


   

%% part4 

%% extract pre-pre graph
clear Distmat;
clear cluster_mat;
Nc = max(res_clustering_pre(:));
% cluster_mat = zeros(100,2,Nc);
for nc = 1:Nc
    
    [row ,col] = find(res_clustering_pre==nc) ;
temp = [row, col];

[t1 ,t2] = size(temp);
cluster_mat(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
end

%% Pairwise distance


dmat = zeros(Nc,Nc);
for c1 = 1:Nc
    temp1 = nonzeros(cluster_mat(:,:,c1));
    [m1,m2] = size(temp1);
    
    cmat_nonzero1 = reshape(temp1,[m1/2,2]);
    for c2 = 1:Nc
        
        if c1<c2
            dmat(c1,c2)=9999;
        else
            
            temp2 = nonzeros(cluster_mat(:,:,c2));
            [m1,m2] = size(temp2);
            
            cmat_nonzero2 = reshape(temp2,[m1/2,2]); 
            Distmat = pdist2(cmat_nonzero1,cmat_nonzero2);
            dmat(c1,c2) = min(Distmat(:));
            
        end
    end
    
end

dmat_prepre = dmat;

%%

temp3 = tril(dmat_prepre);
temp4 = sort(temp3(:));
res_sort_pre = temp4(temp4~=0);
figure, histogram(res_sort_pre)
xlabel('distance')
ylabel('N')
title('pre2pre')

if v==1
    res_sort_pre1 = res_sort_pre;
    dmat_prepre1 = dmat_prepre;
    
    for  ni= 1:max(res_clustering_pre(:))
    index_num_pre1(ni,1) =  sum(nnz(res_clustering_pre==ni));
    end
end


res_center= regionprops(res_clustering_pre,'centroid');
    cent_mat = cat(1, res_center.Centroid);

    figure, imagesc(res_clustering_pre), colormap gray
    hold on 
    plot(cent_mat(:,1), cent_mat(:,2), 'r*') , title('outline & centeroids');
    
    Xfilename = sprintf('M3_L1_LD_Analysis.ascii',v);
    
     dlmwrite(Xfilename,cent_mat);

clear Distmat;
clear cluster_mat;
%% extracting post-post graph
% Nc = max(res_clustering_post(:));
% % cluster_mat = zeros(100,2,Nc);
% for nc = 1:Nc
%     
%     [row ,col] = find(res_clustering_post==nc) ;
% temp = [row, col];
% 
% [t1 ,t2] = size(temp);
% cluster_mat(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
% end

%% Pairwise distance


% dmat = zeros(Nc,Nc);
% for c1 = 1:Nc
%     temp1 = nonzeros(cluster_mat(:,:,c1));
%     [m1,m2] = size(temp1);
%     
%     cmat_nonzero1 = reshape(temp1,[m1/2,2]);
%     for c2 = 1:Nc
%         
%         if c1<c2
%             dmat(c1,c2)=9999;
%         else
%             
%             temp2 = nonzeros(cluster_mat(:,:,c2));
%             [m1,m2] = size(temp2);
%             
%             cmat_nonzero2 = reshape(temp2,[m1/2,2]); 
%             Distmat = pdist2(cmat_nonzero1,cmat_nonzero2);
%             dmat(c1,c2) = min(Distmat(:));
%             
%         end
%     end
%     
% end
% 
% dmat_postpost = dmat;

%%

% temp3 = tril(dmat_postpost);
% temp4 = sort(temp3(:));
% res_sort_post = temp4(temp4~=0);
% figure, histogram(res_sort_post)
% xlabel('distance')
% ylabel('N')
% title('post2post')
% 
% if v==1
%     res_sort_post1 = res_sort_post;
%     dmat_postpost1 = dmat_postpost;
%     
%     for  ni= 1:max(res_clustering_post(:))
%     index_num_post1(ni,1) =  sum(nnz(res_clustering_post==ni));
%     end
%     
% end


clear Distmat;
clear cluster_mat;
%% Extracting distances among clusters

Nc = max(res_synapse(:));
% cluster_mat = zeros(100,2,Nc);
for nc = 1:Nc
    
    [row ,col] = find(res_synapse==nc) ;
temp = [row, col];

[t1 ,t2] = size(temp);
cluster_mat(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
end

%% Pairwise distance


dmat = zeros(Nc,Nc);
for c1 = 1:Nc
    temp1 = nonzeros(cluster_mat(:,:,c1));
    [m1,m2] = size(temp1);
    
    cmat_nonzero1 = reshape(temp1,[m1/2,2]);
    for c2 = 1:Nc
        
        if c1<c2
            dmat(c1,c2)=9999;
        else
            
            temp2 = nonzeros(cluster_mat(:,:,c2));
            [m1,m2] = size(temp2);
            
            cmat_nonzero2 = reshape(temp2,[m1/2,2]); 
            Distmat = pdist2(cmat_nonzero1,cmat_nonzero2);
            dmat(c1,c2) = min(Distmat(:));
            
        end
    end
    
end

dmat_synapse = dmat; 

%%

temp3 = tril(dmat_synapse);
temp4 = sort(temp3(:));
res_sort_synapse = temp4(temp4~=0);
figure, histogram(res_sort_synapse)
xlabel('distance')
ylabel('N')
title('synapse')




if v==1
    res_sort_synapse1 = res_sort_synapse;
    dmat_synapse1 = dmat_synapse;
    temp_dmat1 = dmat_synapse1;
    
       for c1 = 1:Nc
        for c2 = 1:Nc
            if c1<c2
                temp_dmat1(c1,c2) = temp_dmat1(c2,c1);
            end
            if c1==c2
                temp_dmat1(c1,c2) = 9999;
            end
        end
    end

    for e = 1:Nc
        min_dmat1(e,1) = min(temp_dmat1(:,e));
    end
    
    figure, histogram(min_dmat1)
    xlabel('distance')
    ylabel('N')
    title('min_synapse')
    
    for  ni= 1:max(res_synapse(:))
    index_num_synapse1(ni,1) =  sum(nnz(res_synapse==ni));
    end
    
    
    
end


    res_center= regionprops(res_synapse,'centroid');
    cent_mat = cat(1, res_center.Centroid);

    figure, imagesc(res_synapse), colormap gray
    hold on 
    plot(cent_mat(:,1), cent_mat(:,2), 'r*') , title('outline & centeroids');
    
    Xfilename = sprintf('M3_L1_LD_Analysis.ascii',v);
    
    dlmwrite(Xfilename,cent_mat);

clear Distmat;
clear cluster_mat;

%% confirm
figure, imagesc(res_clustering_pre), colormap colorcube, title('pre')
figure, imagesc(res_clustering_post), colormap colorcube, title('post')
figure, imagesc(res_synapse), colormap colorcube, title('synapse')



% test6 = zeros(size(res_synapse));
% 
% for  ni= 1:max(res_synapse(:))
%     test6 = zeros(size(res_synapse));
%     
%     idx5 =  find(res_synapse==ni);
%     
%     test6(idx5) = ni;
%     figure(123), imagesc(test6), colormap colorcube
%     name = sprintf('Lable number is %d', ni);
%     title(name)
%     % ['Acquisition number: ',num2str(ni)])
%      pause(1)
% 
% end

 end