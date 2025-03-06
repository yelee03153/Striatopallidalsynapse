

function res_clustering_gfp = extpregfp(input_body,params)


%% Main purpose: colocalization

% 2D


%% procedure: 
% 1. body data
% data loading (original data)
% Denosing
% Top hat
% Hard thresholding 
% Binarize 
% Opening
% Imfill
% Analysis




%% Denoising
res_denoising1 = imgaussfilt(input_body,params.deno4);
%figure, imagesc(res_denoising2), colormap gray
%title('denoising')


%% Top-hat (enhancement)
% opt1 = strel('disk',6);

opt11 = offsetstrel('ball', 100,100);
res_tophat1 = imtophat(res_denoising1,opt11);

  figure, imagesc(res_tophat1), colormap gray
 title('top-hat')
%figure, bar3(input)
%
%figure, bar3(res_tophat2)
 %figure, imshowpair(input,res_tophat2)


%% Hard thresholding 

res_idx1 = find(res_tophat1> params.hard_lev4);
res_thres1 = zeros(size(res_tophat1));
res_thres1(res_idx1) = res_tophat1(res_idx1);

% figure, imagesc(input), colormap gray
% title('input')
% figure, imagesc(res_thres), colormap gray
% title('thres')

% figure, bar3(res_thres)

% figure, imshowpair(input,res_thres)



%% Binarization

check_threshold1 = graythresh(res_thres1);
% param1 = 1.03; % RSpost 1.03

% res_bw = imbinarize(input, check_threshold* param1);
res_bw11 = imbinarize(res_thres1, check_threshold1* params.otsu4);

% figure, imagesc(res_bw), colormap gray, title('thresholding')

res_bw22 = bwareaopen(res_bw11, params.deno_px_number4);
%  figure, imagesc(res_bw2), colormap gray, title('noise suppresion')



%% Opening 
% opt4 = strel('square',5);
% 
% res_open2 = imopen(res_bw22,opt4);
% res_erode = imerode(res_erode,opt2);
% 
% figure, imagesc(res_open), colormap gray
% title('opening')

%% Imfill
res_fill = imfill(res_bw22, 'holes');
figure, imagesc(res_fill), colormap gray, title('binary image with filled holes');


res_label_pre = bwlabeln(res_fill);

res_clustering_gfp = res_label_pre; 






