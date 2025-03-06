
clear all
close all
%% Main
% part1: extracting pre signal
% part2: extracting post signal
% part3: extracting synapse
% part4: analysis 
%% initialize 
% pre 405 vmat2
params.deno = 0.1;
params.hard_lev1 = 40000;
params.otsu = 50; 
params.deno_px_number = 3;
params.dist = 2;

% body 647 TH
params.deno2 = 0.1;
params.hard_lev2 = 55000;
params.otsu2 = 50; 
params.deno_px_number2 = 200;

% post 594 RFP_Synaptophysin
params.deno3 = 0.1;
params.hard_lev3 = 24000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

% postbody 488 GFP
params.deno4 = 0.1;
params.hard_lev4 = 24000;
params.otsu4 = 25;
params.deno_px_number4 = 20;

folder_name = "E:\1) A2a-Cre 6-OHDA_Synaptophysin_vmat2 GFP RFP TH X63X3 stitching\6-OHDA DL"; % image ?ˆ?Š” ?´?” ?œ„ì¹?
image_list = extractFileList(folder_name); % image ?´ë¦? ì¶”ì¶œ

variable_list = {'colocarea', 'gfparea', 'prenum', 'rfparea', 'rfpnum', 'THarea', 'THnum', 'unionarea', 'unioncolocalvmat2THarea', 'unionTHarea'};
% ?•„?š”?•œ ë³??ˆ˜ ?´ë¦? ë¦¬ìŠ¤?Š¸
variable_list2 = {'min_dmatrr'}
variable_list3 = {'neighbor_RFP_pre'}
variable_list4 = {'neighbor_TH_RFP'}
filename = 'test1.xlsx'; % ?ƒ?„±?•œ ?—‘???ŒŒ?¼ ?´ë¦?(ë¯¸ë¦¬ ?ƒ?„± ?•´?†”?•¼?•¨), ?´ë¯¸ì? ?œ„ì¹˜ì? ?™?¼(?•„?‹ˆë©? ì£¼ì†Œê¹Œì? ? ê¸?)

writematrix(image_list,filename,'Sheet',1,'Range','B4:B'+string(length(image_list)+3)); % ?—‘?? ?ŒŒ?¼?— ?´ë¯¸ì? ?´ë¦„ë“¤ ???¥
writecell(variable_list,filename,'Sheet',1,'Range','C3:L3'); % ?—‘?? ?ŒŒ?¼?— ë³??ˆ˜ ?´ë¦„ë“¤ ???¥

writematrix(image_list.',filename,'Sheet',2,'Range','B2'); 
writecell(variable_list2,filename,'Sheet',2,'Range','A1');

writematrix(image_list.',filename,'Sheet',3,'Range','B2'); 
writecell(variable_list3,filename,'Sheet',3,'Range','A1');

writematrix(image_list.',filename,'Sheet',4,'Range','B2'); 
writecell(variable_list4,filename,'Sheet',4,'Range','A1');

for iteration = 1:length(image_list) % ?´ë¯¸ì??“¤ ë°˜ë³µ ?‹œ?‘?•˜?Š” êµ¬ë¬¸(?´ë¯¸ì? ë¦¬ìŠ¤?Š¸ ê°œìˆ˜ë§Œí¼ ë°˜ë³µ?•˜?‹œ?˜¤ ?¼?Š” ?œ»)
    
 for v = 1
%% Data loading - pre 405 vmat2
% 1: body, 2: post, 3: pre
% input = imread('cell_sample01.tif');

for ci = 1:4
    data1(:,:,ci) = imread(image_list(iteration),ci); 
end

% data_pre = data(1:300,1:300,3:3:end);
% data_pre = data1(1:300,1:300,3:3:end);
data_pre = data1(:,:,1);

%input_pre = data_pre(1:500,1:500,13);
input_pre = data_pre(:,:,v);

% input = input(351:450,231:330);

%% Data loading -body TH 647
% input = imread('cell_sample01.tif');

for ci = 1:4
    data2(:,:,ci) = imread(image_list(iteration),ci);
end

data_body = data2(:,:,4);
%  data_body = data2(1:300,1:300,1:3:end);

%input_body = data_body(1:500,1:500,13);
input_body = data_body(:,:,v);

% input = input(351:450,231:330);
% input = imread('cell_sample01.tif');
%% Data loading - post 594 RFP_synaptophysin
% input = imread('cell_sample01.tif');

for ci = 1:4
    data(:,:,ci) = imread(image_list(iteration),ci);
end

% data_post = data(1:300,1:300,2:3:end);
data_post = data(:,:,3);


%input_post = data_post(1:500,1:500,13);
input_post = data_post(:,:,v);
% input = input(351:450,231:330);
%% Data loading - postbody 488 GFP
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
% function - vmat2, TH, colocal area
res_clustering_pre = extpresynaptophysin(input_pre,input_body ,params);
colocarea = nnz(res_clustering_pre);

res_clustering_TH = extTH(input_body, params);
THarea = nnz(res_clustering_TH);
THnum = max(max(res_clustering_TH));
prenum = max(max(res_clustering_pre));
%% part2
% function - RFP, GFP, union(RFP U GFP) area
res_clustering_gfp = extpregfp(input_postbody, params);
gfparea = nnz(res_clustering_gfp);

res_clustering_rfp = extprerfp(input_post, params);
rfparea = nnz(res_clustering_rfp);
rfpnum = max(max(res_clustering_rfp));

union = unionfn (res_clustering_gfp, res_clustering_rfp);
unionarea = nnz(union);
%% part3
% function - TH area collocal with union
unioncolocalTH = colocalfnbinarize (union, res_clustering_TH);
unionTHarea = nnz(unioncolocalTH);

% function - vmat2, TH, colocal area collocal with union
unioncolocalvmat2TH = colocalfnbinarize (union, res_clustering_pre);
unioncolocalvmat2THarea = nnz(unioncolocalvmat2TH);
%% part4 
% nearest neighbor between other channel
% If you want make other factors' nearest neighbor, just define input1 and
% 2 without changing other sentences.

input1 = res_clustering_rfp;
input2 = unioncolocald2r;

Nc1 = max(input2(:));
% cluster_mat = zeros(100,2,Nc);
for nc = 1:Nc1
    
    [row ,col] = find(input2==nc) ;
temp = [row, col];

[t1 ,t2] = size(temp);
cluster_mat_post(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
end
Nc2 = max(input1(:));
% cluster_mat = zeros(100,2,Nc);
for nc = 1:Nc2
    
    [row ,col] = find(input1==nc) ;
temp = [row, col];

[t1 ,t2] = size(temp);
cluster_mat_pre(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
end
  

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

%

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
%  Making clusters before calculating nearest neighbor between RFP and RFP(without itself)

input1 = res_clustering_rfp;
Nc = max(input1(:));
% max= total image area, nc= search cluster group
% cluster_mat = zeros(100,2,Nc);
% input1= factor that I want to find out nearest neighbor 
for nc = 1:Nc
    
    [row ,col] = find(input1==nc) ;
temp = [row, col];

[t1 ,t2] = size(temp);
cluster_mat(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
end

%% Pairwise distance= nearest neighbor

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

min_dmat(e,1) = min(dmat(:,e));

%

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

    res_center= regionprops(res_synapse,'centroid');
    cent_mat = cat(1, res_center.Centroid);

    figure, imagesc(res_synapse), colormap gray
    hold on 
    plot(cent_mat(:,1), cent_mat(:,2), 'r*') , title('outline & centeroids');
    
    Xfilename = sprintf('test.ascii',v);
    
    dlmwrite(Xfilename,cent_mat);

clear Distmat;
clear cluster_mat;

%% confirm
figure, imagesc(res_clustering_pre), colormap colorcube, title('pre')
figure, imagesc(union), colormap colorcube, title('union')
figure, imagesc(unioncolocalvmat2TH), colormap colorcube, title('unioncolocalvmat2TH')



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

filename = 'test1.xlsx'; % ?ƒ?„±?•œ ?—‘???ŒŒ?¼ ?´ë¦?(ë¯¸ë¦¬ ?ƒ?„± ?•´?†”?•¼?•¨), ?•œë²? ?” ? ?? ?´?œ ?Š” clearvar ?•˜ë©? filename ë³??ˆ˜?„ ì´ˆê¸°?™” ?˜?‹ˆê¹? ê³„ì† ???¥?˜ê²? ?•¨

% ?‹œ?Š¸ ?´ë¦„ì— ?´ë¯¸ì? ?´ë¦? ?„£ê³? ê°ê° ?‹œ?Š¸?„ ?´ë¯¸ì? ?´ë¦„ìœ¼ë¡? ???¥?•  ?•Œ ??
% str1 = extractBefore(image_list(iteration),"_"); % ?´ë¯¸ì? ?´ë¦? ?„ˆë¬? ê¸¸ì–´?„œ ë¶?ë¶? ì¶”ì¶œ
% str2 = extractAfter(image_list(iteration),"3_"); % ?´ë¯¸ì? ?´ë¦? ?„ˆë¬? ê¸¸ì–´?„œ ë¶?ë¶? ì¶”ì¶œ
% sheet_name = str1+' '+str2; % ë¶?ë¶? ì¶”ì¶œ?•œ ?‘ê°?ì§? ê°? ?•©ì¹˜ê¸°
%sheet_list = ['sheet1', 'sheet2', 'sheet3']; % ?„?˜ ê°?

% ?—‘?? ?ŒŒ?¼?— ë³??ˆ˜?“¤ ???¥, iteration+3 ?•œ ?´?œ : 4ë²ˆì§¸ ?—´ë¶??„° ???¥?•˜? ¤ê³?
writematrix(colocarea,filename,'Sheet',1,'Range','C'+string(iteration+3)); 
writematrix(gfparea,filename,'Sheet',1,'Range','D'+string(iteration+3));
writematrix(prenum,filename,'Sheet',1,'Range','E'+string(iteration+3));
writematrix(rfparea,filename,'Sheet',1,'Range','F'+string(iteration+3));
writematrix(rfpnum,filename,'Sheet',1,'Range','G'+string(iteration+3));
writematrix(THarea,filename,'Sheet',1,'Range','H'+string(iteration+3));
writematrix(THnum,filename,'Sheet',1,'Range','I'+string(iteration+3));
writematrix(unionarea,filename,'Sheet',1,'Range','J'+string(iteration+3));
writematrix(unioncolocalvmat2THarea,filename,'Sheet',1,'Range','K'+string(iteration+3));
writematrix(unionTHarea,filename,'Sheet',1,'Range','L'+string(iteration+3));


output_matrix = alphabet_generator();
writematrix(min_dmatrr,filename,'Sheet',2,'Range',output_matrix(iteration+1)+string(3));
writematrix(neighbor_RFP_pre,filename,'Sheet',3,'Range',output_matrix(iteration+1)+string(3));
writematrix(neighbor_TH_RFP,filename,'Sheet',4,'Range',output_matrix(iteration+1)+string(3));

%% Re-initialization
% ?•?— ê°’ë“¤ ë°”ê¿”ì£¼ë©´ ?—¬ê¸°ë„ ë°”ê¿”?•¼?•¨.

clearvars -except iteration % iteration ë¹¼ê³ ?Š” ?‹¤ ?—†?• ê¸? ?œ„?•´?„œ
close all % figure ?‹¤ ?‹«?œ¼? ¤ê³?

% ?—†?•¤ ê°’ë“¤ ?‹¤?‹œ ?„¤? •

params.deno = 0.1;
params.hard_lev1 = 40000;
params.otsu = 50; 
params.deno_px_number = 3;
params.dist = 2;

% body 647 TH
params.deno2 = 0.1;
params.hard_lev2 = 55000;
params.otsu2 = 50; 
params.deno_px_number2 = 200;

% post 594 RFP_Synaptophysin

params.deno3 = 0.1;
params.hard_lev3 = 24000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

% postbody 488 GFP
params.deno4 = 0.1;
params.hard_lev4 = 24000;
params.otsu4 = 25;
params.deno_px_number4 = 20;


folder_name = "E:\1) A2a-Cre 6-OHDA_Synaptophysin_vmat2 GFP RFP TH X63X3 stitching\6-OHDA DL";
image_list = extractFileList(folder_name);

end