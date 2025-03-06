
clear all
close all
%% Main
% part1: extracting pre signal
% part2: extracting post signal
% part3: extracting synapse
% part4: analysis 
%% initialize 
% pre VGAT 647
params.deno = 0.1;
params.hard_lev1 = 33000;
params.otsu = 50; 
params.deno_px_number = 20;
params.dist = 2;

% body RFP 594
params.deno2 = 0.1;
params.hard_lev2 = 30000;
params.otsu2 = 50; 


params.deno_px_number2 = 20;

% post Syt 488
params.deno3 = 0.1;
params.hard_lev3 = 33000;
params.otsu3 = 25;
params.deno_px_number3 = 20;
params.dist3 = 2;

 for v = 1
%% Data loading - pre VGAT

for ci = 1:3
    data1(:,:,ci) = imread('M1_Synaptotagmin7_synap7 RFP vGAT_63X3_R4_LV_Out.tif',ci);
end

% data_pre = data(1:300,1:300,3:3:end);
% data_pre = data1(1:300,1:300,3:3:end);
data_pre = data1(:,:,3);

%input_pre = data_pre(1:500,1:500,13);
input_pre = data_pre(:,:,v);

% input = input(351:450,231:330);




%% Data loading -body RFP

for ci = 1:3
    data2(:,:,ci) = imread('M1_Synaptotagmin7_synap7 RFP vGAT_63X3_R4_LV_Out.tif',ci);
end

data_body = data2(:,:,1);
%  data_body = data2(1:300,1:300,1:3:end);

%input_body = data_body(1:500,1:500,13);
input_body = data_body(:,:,v);

% input = input(351:450,231:330);
% input = imread('cell_sample01.tif');
%% Data loading - post Syt
% input = imread('cell_sample01.tif');

for ci = 1:3
    data(:,:,ci) = imread('M1_Synaptotagmin7_synap7 RFP vGAT_63X3_R4_LV_Out.tif',ci);
end

% data_post = data(1:300,1:300,2:3:end);
data_post = data(:,:,2);


%input_post = data_post(1:500,1:500,13);
input_post = data_post(:,:,v);
% input = input(351:450,231:330);
 %% part1
% function - pre clustering= TH, Vmat2 colocal


res_clustering_pre = extpreyelee (input_pre,input_body ,params);
%% part1-2
% function - RFP area

RFParea = extRFP (input_body,params);

RFPnum = nnz(RFParea);
%% part2
% function - post clustering= RFP area

res_clustering_post = extpostyelee (input_post,params);


%% part3 
% function -  synapse extracting= TH, Vmat2, RFP colocal

[res_synapse] = extsynapseyelee (res_clustering_pre,res_clustering_post);


Prenum = nnz(res_clustering_pre);

Synapsenum = nnz(res_synapse);

Postnum = nnz(res_clustering_post);


    
    for  ni= 1:max(res_clustering_pre(:))
    index_num_pre(ni,1) =  sum(nnz(res_clustering_pre==ni));
    end

    for  ni= 1:max(res_clustering_post(:))
    index_num_post(ni,1) =  sum(nnz(res_clustering_post==ni));
    end
    
    for  ni= 1:max(res_synapse(:))
    index_num_synapse(ni,1) =  sum(nnz(res_synapse==ni));
    end
    
    
    

clear Distmat;
clear cluster_mat;

%% confirm
figure, imagesc(res_clustering_pre), colormap colorcube, title('pre')
figure, imagesc(res_clustering_post), colormap colorcube, title('post')
figure, imagesc(res_synapse), colormap colorcube, title('synapse')




 end