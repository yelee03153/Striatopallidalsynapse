
clear all
close all
%% Main
% part1: extracting pre signal
% part2: extracting post signal
% part3: extracting synapse
% part4: analysis 
%% initialize 
% pre (RFP)
params.deno = 0.1;
params.hard_lev1 = 30000;
params.otsu = 50; 
params.deno_px_number = 20;
params.dist = 2;

% body (GFP)
params.deno2 = 0.1;
params.hard_lev2 = 30000;
params.otsu2 = 50; 
params.deno_px_number2 = 20;

% post (D2R)
params.deno3 = 0.1;
params.hard_lev3 = 10000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

folder_name = "E:\3) A2a-Cre 6-OHDA_Synaptophysin_GFP RFP D2R X63X3 stack12\6-OHDA DL"; 
image_list = extractFileList(folder_name); 

variable_list = {'colocarea1', 'colocnum1', 'd2rarea1', 'd2rnum1', 'gfparea1', 'gfpd2rarea1', 'gfpd2rnum1', 'rfparea1','rfpnum1','rfpd2rarea1', 'rfpd2rnum1','unionarea1','unionnum1', 'uniond2rarea1', 'uniond2rnum1','colocarea6', 'colocnum6', 'd2rarea6', 'd2rnum6', 'gfparea6', 'gfpd2rarea6', 'gfpd2rnum6', 'rfparea6','rfpnum6','rfpd2rarea6', 'rfpd2rnum6','unionarea6','unionnum6', 'uniond2rarea6', 'uniond2rnum6', 'colocarea11', 'colocnum11', 'd2rarea11', 'd2rnum11', 'gfparea11', 'gfpd2rarea11', 'gfpd2rnum11', 'rfparea11','rfpnum11','rfpd2rarea11', 'rfpd2rnum11','unionarea11','unionnum11', 'uniond2rarea11', 'uniond2rnum11'};
variable_list2 = {'dmat_nnr1'}
variable_list3 = {'dmat_nnr6'}
variable_list4 = {'dmat_nnr11'}


filename = 'test1.xlsx'; 
writematrix(image_list,filename,'Sheet',1,'Range','B4'); 
writecell(variable_list,filename,'Sheet',1,'Range','C3:Z3,AA,AB,AC,AD,AE,AF,AG,AH,AI,AJ,AK,AL,AM,AN,AO,AP,AQ,AR,AS,AT,AU');

writematrix(image_list.',filename,'Sheet',2,'Range','B2'); 
writecell(variable_list2,filename,'Sheet',2,'Range','A1');

writematrix(image_list.',filename,'Sheet',3,'Range','B2'); 
writecell(variable_list3,filename,'Sheet',3,'Range','A1');

writematrix(image_list.',filename,'Sheet',4,'Range','B2'); 
writecell(variable_list4,filename,'Sheet',4,'Range','A1');

for iteration = 1:length(image_list) 


 for v = 1:5:12
%% Data loading - pre
% 1: body, 2: post, 3: pre
% input = imread('cell_sample01.tif');

for ci = 1:36
    data1(:,:,ci) = imread(image_list(iteration),ci);
end

% data_pre = data(1:300,1:300,3:3:end);
% data_pre = data1(1:300,1:300,3:3:end);
data_pre = data1(:,:,2:3:end);

%input_pre = data_pre(1:500,1:500,13);
input_pre = data_pre(:,:,v);

% input = input(351:450,231:330);




%% Data loading -body
% input = imread('cell_sample01.tif');

for ci = 1:36
    data2(:,:,ci) = imread('M1 4509 6-OHDA Left_GFP RFP D2R_X63X3_LD1_Out.tif',ci);
end

data_body = data2(:,:,1:3:end);
%  data_body = data2(1:300,1:300,1:3:end);

%input_body = data_body(1:500,1:500,13);
input_body = data_body(:,:,v);

% input = input(351:450,231:330);
% input = imread('cell_sample01.tif');
%% Data loading - post
% input = imread('cell_sample01.tif');

for ci = 1:36
    data(:,:,ci) = imread('M1 4509 6-OHDA Left_GFP RFP D2R_X63X3_LD1_Out.tif',ci);
end

% data_post = data(1:300,1:300,2:3:end);
data_post = data(:,:,3:3:end);


%input_post = data_post(1:500,1:500,13);
input_post = data_post(:,:,v);
% input = input(351:450,231:330);
%% part1
% function - gfp, rfp, colocal area/number


res_clustering_pre = extpre(input_pre,input_body ,params);

if v == 1
    colocarea1 = nnz(res_clustering_pre);
    colocnum1 = max(max(res_clustering_pre));
end
if v == 6
    colocarea6 = nnz(res_clustering_pre);
    colocnum6 = max(max(res_clustering_pre));
end
if v == 11
    colocarea11 = nnz(res_clustering_pre);
    colocnum11 = max(max(res_clustering_pre));
end

res_clustering_gfp = extpregfp(input_body, params);

if v == 1
    gfparea1 = nnz(res_clustering_gfp);
    gfpnum1 = max(max(res_clustering_gfp));
end
if v == 6
    gfparea6 = nnz(res_clustering_gfp);
end
if v == 11
    gfparea11 = nnz(res_clustering_gfp);
end

res_clustering_rfp = extprerfp(input_pre, params); 

if v == 1
    rfparea1 = nnz(res_clustering_rfp);
    rfpnum1 = max(max(res_clustering_rfp));
end
if v == 6
    rfparea6 = nnz(res_clustering_rfp);
     rfpnum6 = max(max(res_clustering_rfp));
end
if v == 11
    rfparea11 = nnz(res_clustering_rfp);
     rfpnum11 = max(max(res_clustering_rfp));
end
%% part2
% function - d2r, union area/number

res_clustering_post = extpostyelee (input_post,params);

if v == 1
    d2rarea1 = nnz(res_clustering_post);
    d2rnum1 = max(max(res_clustering_post));
end
if v == 6
    d2rarea6 = nnz(res_clustering_post);
    d2rnum6 = max(max(res_clustering_post));
end
if v == 11
    d2rarea11 = nnz(res_clustering_post);
    d2rnum11 = max(max(res_clustering_post));
end

union = unionfn (res_clustering_gfp, res_clustering_rfp);

if v == 1
    unionarea1 = nnz(union);
    unionnum1 = max(max(union));
end
if v == 6
    unionarea6 = nnz(union);
    unionnum6 = max(max(union));
end
if v == 11
    unionarea11 = nnz(union);
    unionnum11 = max(max(union));
end

%% part4
% gfp/rfp/union colocaled with d2r

rfpcolocald2r = colocalfnbinarize (res_clustering_rfp, res_clustering_post);

if v == 1
    rfpd2rarea1 = nnz(rfpcolocald2r);
    rfpd2rnum1 = max(max(rfpcolocald2r));
end
if v == 6
    rfpd2rarea6 = nnz(rfpcolocald2r);
    rfpd2rnum6 = max(max(rfpcolocald2r));
end
if v == 11
    rfpd2rarea11 = nnz(rfpcolocald2r);
    rfpd2rnum11 = max(max(rfpcolocald2r));
end

gfpcolocald2r = colocalfnbinarize (res_clustering_gfp, res_clustering_post);

if v == 1
    gfpd2rarea1 = nnz(gfpcolocald2r);
    gfpd2rnum1 = max(max(gfpcolocald2r));
end
if v == 6
    gfpd2rarea6 = nnz(gfpcolocald2r);
    gfpd2rnum6 = max(max(gfpcolocald2r));
end
if v == 11
    gfpd2rarea11 = nnz(gfpcolocald2r);
    gfpd2rnum11 = max(max(gfpcolocald2r));
end

unioncolocald2r = colocalfnbinarize (union, res_clustering_post);

if v == 1
    uniond2rarea1 = nnz(unioncolocald2r);
    uniond2rnum1 = max(max(unioncolocald2r));
end
if v == 6
    uniond2rarea6 = nnz(unioncolocald2r);
    uniond2rnum6 = max(max(unioncolocald2r));
end
if v == 11
    uniond2rarea11 = nnz(unioncolocald2r);
    uniond2rnum11 = max(max(unioncolocald2r));
end

%% part5
% nearest neighbor
% stack 'for' sentence goes from line 30
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

%% Pairwise distance


dmat = zeros(Nc1,Nc2);
for c1 = 1:Nc1
    temp1 = nonzeros(cluster_mat_post(:,:,c1));
    [m1,m2] = size(temp1);
    
    cmat_nonzero1 = reshape(temp1,[m1/2,2]);
    for c2 = 1:Nc2
        
%         if c1<c2
%             dmat(c1,c2)=9999;
%         else
            
            temp2 = nonzeros(cluster_mat_pre(:,:,c2));
            [m1,m2] = size(temp2);
            
            cmat_nonzero2 = reshape(temp2,[m1/2,2]); 
            Distmat = pdist2(cmat_nonzero1,cmat_nonzero2);
            dmat(c1,c2) = min(Distmat(:));
            
%         end
    end
    
end

temp_dmat = dmat;

 
    for e = 1:Nc2
        min_dmat(e,1) = min(temp_dmat(:,e));
    end


if v==1
    dmat_nnr1 = min_dmat;
end
if v==6
    dmat_nnr6 = min_dmat;
end
if v==11
    dmat_nnr11 = min_dmat; 
end


clear Distmat;
clear temp_dmat;
clear min_dmat;
clear cluster_mat_pre;
clear cluster_mat_post;
 end
 
 
 
 %% Save data into Excel file
filename = 'test1.xlsx'; 
writematrix(colocarea1,filename,'Sheet',1,'Range','C'+string(iteration+3)); 
writematrix(colocnum1,filename,'Sheet',1,'Range','D'+string(iteration+3));
writematrix(d2rarea1,filename,'Sheet',1,'Range','E'+string(iteration+3));
writematrix(d2rnum1,filename,'Sheet',1,'Range','F'+string(iteration+3));
writematrix(gfparea1,filename,'Sheet',1,'Range','G'+string(iteration+3));
writematrix(gfpd2rarea1,filename,'Sheet',1,'Range','H'+string(iteration+3));
writematrix(gfpd2rnum1,filename,'Sheet',1,'Range','I'+string(iteration+3));
writematrix(rfparea1,filename,'Sheet',1,'Range','J'+string(iteration+3));
writematrix(rfpnum1,filename,'Sheet',1,'Range','K'+string(iteration+3));
writematrix(rfpd2rarea1,filename,'Sheet',1,'Range','L'+string(iteration+3));
writematrix(rfpd2rnum1,filename,'Sheet',1,'Range','M'+string(iteration+3));
writematrix(unionarea1,filename,'Sheet',1,'Range','N'+string(iteration+3));
writematrix(unionnum1,filename,'Sheet',1,'Range','O'+string(iteration+3));
writematrix(uniond2rarea1,filename,'Sheet',1,'Range','P'+string(iteration+3));
writematrix(uniond2rnum1,filename,'Sheet',1,'Range','Q'+string(iteration+3));

writematrix(colocarea6,filename,'Sheet',1,'Range','R'+string(iteration+3)); 
writematrix(colocnum6,filename,'Sheet',1,'Range','S'+string(iteration+3));
writematrix(d2rarea6,filename,'Sheet',1,'Range','T'+string(iteration+3));
writematrix(d2rnum6,filename,'Sheet',1,'Range','U'+string(iteration+3));
writematrix(gfparea6,filename,'Sheet',1,'Range','V '+string(iteration+3));
writematrix(gfpd2rarea6,filename,'Sheet',1,'Range','W'+string(iteration+3));
writematrix(gfpd2rnum6,filename,'Sheet',1,'Range','X'+string(iteration+3));
writematrix(rfparea6,filename,'Sheet',1,'Range','Y'+string(iteration+3));
writematrix(rfpnum6,filename,'Sheet',1,'Range','Z'+string(iteration+3));
writematrix(rfpd2rarea6,filename,'Sheet',1,'Range','AA'+string(iteration+3));
writematrix(rfpd2rnum6,filename,'Sheet',1,'Range','AB'+string(iteration+3));
writematrix(unionarea6,filename,'Sheet',1,'Range','AC'+string(iteration+3));
writematrix(unionnum6,filename,'Sheet',1,'Range','AD'+string(iteration+3));
writematrix(uniond2rarea6,filename,'Sheet',1,'Range','AE'+string(iteration+3));
writematrix(uniond2rnum6,filename,'Sheet',1,'Range','AF'+string(iteration+3));

writematrix(colocarea11,filename,'Sheet',1,'Range','AG'+string(iteration+3)); 
writematrix(colocnum11,filename,'Sheet',1,'Range','AH'+string(iteration+3));
writematrix(d2rarea11,filename,'Sheet',1,'Range','AI'+string(iteration+3));
writematrix(d2rnum11,filename,'Sheet',1,'Range','AJ'+string(iteration+3));
writematrix(gfparea11,filename,'Sheet',1,'Range','AK'+string(iteration+3));
writematrix(gfpd2rarea11,filename,'Sheet',1,'Range','AL'+string(iteration+3));
writematrix(gfpd2rnum11,filename,'Sheet',1,'Range','AM'+string(iteration+3));
writematrix(rfparea11,filename,'Sheet',1,'Range','AN'+string(iteration+3));
writematrix(rfpnum11,filename,'Sheet',1,'Range','AO'+string(iteration+3));
writematrix(rfpd2rarea11,filename,'Sheet',1,'Range','AP'+string(iteration+3));
writematrix(rfpd2rnum11,filename,'Sheet',1,'Range','AQ'+string(iteration+3));
writematrix(unionarea11,filename,'Sheet',1,'Range','AR'+string(iteration+3));
writematrix(unionnum11,filename,'Sheet',1,'Range','AS'+string(iteration+3));
writematrix(uniond2rarea11,filename,'Sheet',1,'Range','AT'+string(iteration+3));
writematrix(uniond2rnum11,filename,'Sheet',1,'Range','AU'+string(iteration+3));


output_matrix = alphabet_generator();
writematrix(dmat_nnr1,filename,'Sheet',2,'Range',output_matrix(iteration+1)+string(3));
writematrix(dmat_nnr6,filename,'Sheet',3,'Range',output_matrix(iteration+1)+string(3));
writematrix(dmat_nnr11,filename,'Sheet',4,'Range',output_matrix(iteration+1)+string(3));

%% Re-initialization


clearvars -except iteration 
close all 
% pre (RFP)
params.deno = 0.1;
params.hard_lev1 = 30000;
params.otsu = 50; 
params.deno_px_number = 20;
params.dist = 2;

% body (GFP)
params.deno2 = 0.1;
params.hard_lev2 = 30000;
params.otsu2 = 50; 
params.deno_px_number2 = 20;

% post (D2R)
params.deno3 = 0.1;
params.hard_lev3 = 10000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

folder_name = "E:\3) A2a-Cre 6-OHDA_Synaptophysin_GFP RFP D2R X63X3 stack12\6-OHDA DL"; 
image_list = extractFileList(folder_name); 

end