
function res_clustering_rfp = extRFPyelee (input_body,params)


%% Main purpose: to find out distance among signals
% 2D
% outline & centroids, and distance & rearrangement 


%% procedure: 
% data loading (original data)
% Denosing
% Top hat
% Hard thresholding 
% Binarize 
% Opening
% Imfill
% Analysis


%% Denoising
res_denoising = imgaussfilt(input_body,params.deno4);
%  figure, imagesc(res_denoising), colormap gray
% title('denoising')


%% Top-hat (enhancement)
% opt1 = strel('disk',6);

opt1 = offsetstrel('ball', 100,100);
res_tophat = imtophat(res_denoising,opt1);

%  figure, imagesc(res_tophat), colormap gray
% title('top-hat')
% figure, bar3(input)
% 
% figure, bar3(res_tophat)
% figure, imshowpair(input,res_tophat)


%% Hard thresholding 

res_idx = find(res_tophat> params.hard_lev4);
res_thres = zeros(size(res_tophat));
res_thres(res_idx) = res_tophat(res_idx);

% figure, imagesc(input), colormap gray
% title('input')
% figure, imagesc(res_thres), colormap gray
% title('thres')

% figure, bar3(res_thres)

% figure, imshowpair(input,res_thres)



%% Binarization

check_threshold = graythresh(res_thres);
% param1 = 1.03; % RSpost 1.03

% res_bw = imbinarize(input, check_threshold* param1);
res_bw = imbinarize(res_thres, check_threshold* params.otsu4);

% figure, imagesc(res_bw), colormap gray, title('thresholding')

res_bw2 = bwareaopen(res_bw, params.deno_px_number4);
%  figure, imagesc(res_bw2), colormap gray, title('noise suppresion')



%% Opening 
opt2 = strel('square',1);

% res_open = imopen(res_bw2,opt2);
% res_erode = imerode(res_erode,opt2);
% 
% figure, imagesc(res_open), colormap gray
% title('opening')

%% Imfill
res_fill = imfill(res_bw2, 'holes');
% figure, imagesc(res_fill), colormap gray, title('binary image with filled holes');


%% Outline
%  

% SE_outline = false(3,3);
% SE_outline(2,2) = 1;

% SE_outline = strel('disk',1,8); % erosion size becomes smaller when the second parameter increases
% res_outline = bwperim(res_fill);
% % res_outline = res_fill- imerode(res_fill,SE_outline);
% res_outline = logical(res_outline);
% % figure, imagesc(res_outline),colormap gray, title('outline image');


% SE_outline = strel('disk',1,8); % erosion size becomes smaller when the second parameter increases
% res_outline = bwperim(res_fill);
% res_outline = logical(res_outline);



%% Lableing (cluster)

res_label_RFP = bwlabeln(res_fill);
res_clustering_rfp = res_label_RFP;
% 
% figure, imagesc(res_label_post), colormap colorcube
% title('labeling cluster')

%% Clustering

%% Extracting indices between cluster 

% Nc = max(res_label_TH(:));
% % cluster_mat = zeros(100,2,Nc);
% for nc = 1:Nc
%     
%     [row ,col] = find(res_label_TH==nc) ;
% temp = [row, col];
% 
% [t1 ,t2] = size(temp);
% cluster_mat(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
% end
% 
% % %% Measuring distance (single)
% % 
% % temp1 = nonzeros(cluster_mat(:,:,1));
% % [m1,m2] = size(temp1);
% % 
% % cmat_nonzero1 = reshape(temp1,[m1/2,2]);
% % 
% % 
% % temp2 = nonzeros(cluster_mat(:,:,2));
% % [m1,m2] = size(temp2);
% % 
% % cmat_nonzero2 = reshape(temp2,[m1/2,2]);
% % 
% % 
% % %%
% % Dtest = pdist2(cmat_nonzero1,cmat_nonzero2,'euclidean','Smallest',1);
% % Dtest2 = pdist2(cmat_nonzero1,cmat_nonzero2);
% % 
% 
% 
% %% Pairwise distance
% 
% 
% dmat = zeros(Nc,Nc);
% for c1 = 1:Nc
%     temp1 = nonzeros(cluster_mat(:,:,c1));
%     [m1,~] = size(temp1);
%     
%     cmat_nonzero1 = reshape(temp1,[m1/2,2]);
%     for c2 = 1:Nc
%         
% %         if c1<c2
% %             dmat(c1,c2)=9999;
% %         else
%             
%             temp2 = nonzeros(cluster_mat(:,:,c2));
%             [m1,~] = size(temp2);
%             
%             cmat_nonzero2 = reshape(temp2,[m1/2,2]); 
%             Distmat = pdist2(cmat_nonzero1,cmat_nonzero2);
%             dmat(c1,c2) = min(Distmat(:));
%             
% %         end
%     end
%     
% end
% 
% 
% %% Find out the minimum values
%  %dvec = dmat(:);
% %  ind = find(dvec <= d_criterion);
% ind2 = find(dmat <= params.dist3);
% 
% testmat = zeros(size(dmat));
% % D1mat = eye(size(testmat));
% testmat(ind2)=1;
% % testmat = testmat;%+D1mat;
% 
% %% 
% % save_labelv2 = testmat;
% % save_labelv2 = save_label;
% save_labelv2 = testmat;
% 
% [q1, q2] = size(save_labelv2);
% 
% 
% for w1 = 1:q1
%     for w2 = 1:q2
%         
%         if w1==w2
%         else
%            val1= sum(save_labelv2(:,w1).*save_labelv2(:,w2));
%             if val1 ==0
%                 
%             else
%                 save_labelv2(:,w1) = save_labelv2(:,w1)+save_labelv2(:,w2);
%                 save_labelv2(:,w2)=0;
%                 
%             end
%             
%             
%         end
%     end
%     
% end
% 
% %%
% 
% col2 = sum(save_labelv2,1);
% new_label = save_labelv2(:, find(col2));
% 
% 
% 
% %%
% THarea = zeros(size(res_label_TH));
% [~,n22] = size(new_label);
% 
% for ni= 1:n22
%     temp_idx = find( new_label(:,ni));
%     
%     for qi = 1: length(temp_idx)
% %         idx4 =  find(res_clustering==temp_idx(qi));
%         idx4 =  find(res_label_TH==temp_idx(qi));
% 
%         THarea(idx4) = ni;
%     end
%     
% end

    
    
end



