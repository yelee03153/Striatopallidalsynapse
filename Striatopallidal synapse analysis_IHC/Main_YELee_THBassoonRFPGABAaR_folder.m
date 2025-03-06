
clear all
close all
%% Main
% part1: extracting pre signal
% part2: extracting post signal
% part3: extracting synapse
% part4: analysis 
%% initialize 

% pre 488 Bassoon
params.deno = 0.1;
params.hard_lev1 = 35000;
params.otsu = 50; 
params.deno_px_number = 5;
params.dist = 2;

% body 405 TH
params.deno2 = 0.1;
params.hard_lev2 = 30000;
params.otsu2 = 50; 
params.deno_px_number2 = 20;

% postpost 647 GABAaR
params.deno3 = 0.1;
params.hard_lev3 = 38000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

% postbody 594 RFP
params.deno4 = 0.1;
params.hard_lev4 = 40000;
params.otsu4 = 25;
params.deno_px_number4 = 20;

folder_name = "E:\sample"; 
image_list = extractFileList_imagenameincluding(folder_name); 


variable_list = {'colocarea', 'prenum', 'postcolocarea', 'postnum', 'synapsearea', 'synapsenum', 'THarea', 'THnum', 'rfparea'};
variable_list2 = {'min_dmat1'}
variable_list3 = {'min_dmatPost'}
variable_list4 = {'min_dmatPre'}

filename = 'test1.xlsx'; 

writematrix(image_list,filename,'Sheet',1,'Range','B4'); 
writecell(variable_list,filename,'Sheet',1,'Range','C3:K3');

writematrix(image_list,filename,'Sheet',2,'Range','B2'); 
writecell(variable_list2,filename,'Sheet',2,'Range','A1');

writematrix(image_list,filename,'Sheet',3,'Range','B2'); 
writecell(variable_list3,filename,'Sheet',3,'Range','A1');

writematrix(image_list,filename,'Sheet',4,'Range','B2'); 
writecell(variable_list4,filename,'Sheet',4,'Range','A1');


for iteration = 1:length(image_list) 


 for v = 1
%% Data loading - pre
% 1: body, 2: post, 3: pre
% input = imread('cell_sample01.tif');

for ci = 1:4
    data1(:,:,ci) = imread(image_list(iteration),ci);
end

% data_pre = data(1:300,1:300,3:3:end);
% data_pre = data1(1:300,1:300,3:3:end);
data_pre = data1(:,:,3);

%input_pre = data_pre(1:500,1:500,13);
input_pre = data_pre(:,:,v);

% input = input(351:450,231:330);




%% Data loading -body
% input = imread('cell_sample01.tif');

for ci = 1:4
    data2(:,:,ci) = imread(image_list(iteration),ci);
end

data_body = data2(:,:,1);
%  data_body = data2(1:300,1:300,1:3:end);

%input_body = data_body(1:500,1:500,13);
input_body = data_body(:,:,v);

% input = input(351:450,231:330);
% input = imread('cell_sample01.tif');
%% Data loading - postpost
% input = imread('cell_sample01.tif');

for ci = 1:4
    data(:,:,ci) = imread(image_list(iteration),ci);
end

% data_post = data(1:300,1:300,2:3:end);
data_postpost = data(:,:,4);


%input_post = data_post(1:500,1:500,13);
input_postpost = data_postpost(:,:,v);
% input = input(351:450,231:330);
%% Data loading - postbody
% input = imread('cell_sample01.tif');

for ci = 1:4
    data(:,:,ci) = imread(image_list(iteration),ci);
end

% data_post = data(1:300,1:300,2:3:end);
data_postbody = data(:,:,2);


%input_post = data_post(1:500,1:500,13);
input_postbody = data_postbody(:,:,v);
% input = input(351:450,231:330);
%% part1
% function - pre clustering
res_clustering_pre = extpresynaptophysin(input_pre,input_body ,params);
colocarea = nnz(res_clustering_pre);

res_clustering_TH = extTH(input_body, params);
THarea = nnz(res_clustering_TH);
THnum = max(max(res_clustering_TH));
prenum = max(max(res_clustering_pre));
%% part2
% function - post clustering

res_clustering_post = extpostcoloc(input_postpost,input_postbody,params);
postcolocarea = nnz(res_clustering_post);
postnum = max(max(res_clustering_post));

res_clustering_rfp = extRFPyelee(input_postbody, params);
rfparea = nnz(res_clustering_rfp);

%% part3 
% function -  synapse extracting

d =2;
[res_synapse,dmat_prepost] = extsynapv1 (res_clustering_pre,res_clustering_post,d);

temp3 = tril(dmat_prepost);
temp4 = sort(temp3(:));
res_sort2_prepost = temp4(temp4~=0);
figure, histogram(res_sort2_prepost)
xlabel('distance')
ylabel('N')
title('pre2post')

if v==1
    res_sort_prepost1 = res_sort2_prepost;
    dmat_prepost1 = dmat_prepost;
end
if v==2
    res_sort_prepost2 = res_sort2_prepost;
    dmat_prepost2 = dmat_prepost;
end

if v==3
    res_sort_prepost3 = res_sort2_prepost;
    dmat_prepost3 = dmat_prepost; 
end

synapsearea = nnz(res_synapse);
synapsenum = max(max(res_synapse));

NPpost = max(max(res_clustering_post));

for f = 1:NPpost
        min_dmatPost(1,f) = min(dmat_prepost1(f,:));
end
    
NPpre = max(max(res_clustering_pre));

for g = 1:NPpre
        min_dmatPre(1,g) = min(dmat_prepost1(:,g));
end


   

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
if v==2
    res_sort_pre2 = res_sort_pre;
    dmat_prepre2 = dmat_prepre;
    
    for  ni= 1:max(res_clustering_pre(:))
    index_num_pre2(ni,1) =  sum(nnz(res_clustering_pre==ni));
    end
end

if v==3
    res_sort_pre3 = res_sort_pre;
    dmat_prepre3 = dmat_prepre; 
    
    for  ni= 1:max(res_clustering_pre(:))
    index_num_pre3(ni,1) =  sum(nnz(res_clustering_pre==ni));
    end
end

% res_center= regionprops(res_clustering_pre,'centroid');
%     cent_mat = cat(1, res_center.Centroid);
% 
%     figure, imagesc(res_clustering_pre), colormap gray
%     hold on 
%     plot(cent_mat(:,1), cent_mat(:,2), 'r*') , title('outline & centeroids');
%     
%     Xfilename = sprintf('M1S1MDL%dpre.ascii',v);
%     
%      dlmwrite(Xfilename,cent_mat);

clear Distmat;
clear cluster_mat;
%% extracting post-post graph
Nc = max(res_clustering_post(:));
% cluster_mat = zeros(100,2,Nc);
for nc = 1:Nc
    
    [row ,col] = find(res_clustering_post==nc) ;
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

dmat_postpost = dmat;

%%

temp3 = tril(dmat_postpost);
temp4 = sort(temp3(:));
res_sort_post = temp4(temp4~=0);
figure, histogram(res_sort_post)
xlabel('distance')
ylabel('N')
title('post2post')

if v==1
    res_sort_post1 = res_sort_post;
    dmat_postpost1 = dmat_postpost;
    
    for  ni= 1:max(res_clustering_post(:))
    index_num_post1(ni,1) =  sum(nnz(res_clustering_post==ni));
    end
    
end
if v==2
    res_sort_post2 = res_sort_post;
    dmat_postpost2 = dmat_postpost;
    
    for  ni= 1:max(res_clustering_post(:))
    index_num_post2(ni,1) =  sum(nnz(res_clustering_post==ni));
    end
end
if v==3
    res_sort_post3 = res_sort_post;
    dmat_postpost3 = dmat_postpost;
    
    for  ni= 1:max(res_clustering_post(:))
    index_num_post3(ni,1) =  sum(nnz(res_clustering_post==ni));
    end
end

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
        min_dmat1(1,e) = min(temp_dmat1(:,e));
    end
    
    figure, histogram(min_dmat1)
    xlabel('distance')
    ylabel('N')
    title('min_synapse')
    
    for  ni= 1:max(res_synapse(:))
    index_num_synapse1(ni,1) =  sum(nnz(res_synapse==ni));
    end
    
    
    
end
if v==2
    res_sort_synapse2 = res_sort_synapse;
    dmat_synapse2 = dmat_synapse;
     temp_dmat2 = dmat_synapse2;
    
       for c1 = 1:Nc
        for c2 = 1:Nc
            if c1<c2
                temp_dmat2(c1,c2) = temp_dmat2(c2,c1);
            end
            if c1==c2
                temp_dmat2(c1,c2) = 9999;
            end
        end
    end

    for e = 1:Nc
        min_dmat2(e,1) = min(temp_dmat2(:,e));
    end
    
    figure, histogram(min_dmat2)
    xlabel('distance')
    ylabel('N')
    title('min_synapse')
    
    for  ni= 1:max(res_synapse(:))
    index_num_synapse2(ni,1) =  sum(nnz(res_synapse==ni));
    end
    
end
if v==3
    res_sort_synapse3 = res_sort_synapse;
    dmat_synapse3 = dmat_synapse;
     temp_dmat3 = dmat_synapse3;
    
       for c1 = 1:Nc
        for c2 = 1:Nc
            if c1<c2
                temp_dmat3(c1,c2) = temp_dmat3(c2,c1);
            end
            if c1==c2
                temp_dmat3(c1,c2) = 9999;
            end
        end
    end

    for e = 1:Nc
        min_dmat3(e,1) = min(temp_dmat3(:,e));
    end
    
    figure, histogram(min_dmat3)
    xlabel('distance')
    ylabel('N')
    title('min_synapse')
    
    for  ni= 1:max(res_synapse(:))
    index_num_synapse3(ni,1) =  sum(nnz(res_synapse==ni));
    end
    
end

%     res_center= regionprops(res_synapse,'centroid');
%     cent_mat = cat(1, res_center.Centroid);
% 
%     figure, imagesc(res_synapse), colormap gray
%     hold on 
%     plot(cent_mat(:,1), cent_mat(:,2), 'r*') , title('outline & centeroids');
%     
%     Xfilename = sprintf('test.ascii',v);
%     
%     dlmwrite(Xfilename,cent_mat);

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
 %% Save data into Excel file
filename = 'test1.xlsx'; 
writematrix(colocarea,filename,'Sheet',1,'Range','C'+string(iteration+3)); 
writematrix(prenum,filename,'Sheet',1,'Range','D'+string(iteration+3));
writematrix(postcolocarea,filename,'Sheet',1,'Range','E'+string(iteration+3));
writematrix(postnum,filename,'Sheet',1,'Range','F'+string(iteration+3));
writematrix(synapsearea,filename,'Sheet',1,'Range','G'+string(iteration+3));
writematrix(synapsenum,filename,'Sheet',1,'Range','H'+string(iteration+3));
writematrix(THarea,filename,'Sheet',1,'Range','I'+string(iteration+3));
writematrix(THnum,filename,'Sheet',1,'Range','J'+string(iteration+3));
writematrix(rfparea,filename,'Sheet',1,'Range','K'+string(iteration+3));

writematrix(min_dmat1,filename,'Sheet',2,'Range','C'+string(iteration+1));
writematrix(min_dmatPost,filename,'Sheet',3,'Range','C'+string(iteration+1));
writematrix(min_dmatPre,filename,'Sheet',4,'Range','C'+string(iteration+1));

%% Re-initialization


clearvars -except iteration 
close all 


% pre 488 Bassoon
params.deno = 0.1;
params.hard_lev1 = 35000;
params.otsu = 50; 
params.deno_px_number = 5;
params.dist = 2;

% body 405 TH
params.deno2 = 0.1;
params.hard_lev2 = 30000;
params.otsu2 = 50; 
params.deno_px_number2 = 20;

% postpost 647 GABAaR
params.deno3 = 0.1;
params.hard_lev3 = 38000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

% postbody 594 RFP
params.deno4 = 0.1;
params.hard_lev4 = 40000;
params.otsu4 = 25;
params.deno_px_number4 = 20;

folder_name = "E:\sample";
image_list = extractFileList_imagenameincluding(folder_name);

end