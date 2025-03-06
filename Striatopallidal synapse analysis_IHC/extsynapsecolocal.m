
function [res_cluster_synapse,dmat] = extsynapv1 (input1,input2,d_criterion2)

%input1 : presynapse
%input2: postsynapse
%d_criterion2: criterion for synapse (distance between pre & post)

res_presynapse = colocalfn (input_pre,res_fill_body);


%% Extracting indices between cluster 

% Nc1 = max(input2(:));
% % cluster_mat = zeros(100,2,Nc);
% for nc = 1:Nc1
%     
%     [row ,col] = find(input2==nc) ;
% temp = [row, col];
% 
% [t1 ,t2] = size(temp);
% cluster_mat_post(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
% end
% Nc2 = max(input1(:));
% % cluster_mat = zeros(100,2,Nc);
% for nc = 1:Nc2
%     
%     [row ,col] = find(input1==nc) ;
% temp = [row, col];
% 
% [t1 ,t2] = size(temp);
% cluster_mat_pre(1:t1,1:t2,nc) = temp; %cluster_mat: row, col, Number (label idx) of clusters
% end
% 
% %% Pairwise distance
% 
% 
% dmat = zeros(Nc1,Nc2);
% for c1 = 1:Nc1
%     temp1 = nonzeros(cluster_mat_post(:,:,c1));
%     [m1,m2] = size(temp1);
%     
%     cmat_nonzero1 = reshape(temp1,[m1/2,2]);
%     for c2 = 1:Nc2
%         
% %         if c1<c2
% %             dmat(c1,c2)=9999;
% %         else
%             
%             temp2 = nonzeros(cluster_mat_pre(:,:,c2));
%             [m1,m2] = size(temp2);
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
% %  dvec = dmat(:);
% %  ind = find(dvec <= d_criterion);
%  ind2 = find(dmat <= d_criterion2);
% 
% testmat = zeros(size(dmat));
% % D1mat = eye(size(testmat));
% testmat(ind2)=1;
% % testmat = testmat;%+D1mat;
% 
% 
% %% synapse clustering
% 
% res_clustering_temp = zeros(size(input1));
% [n11,n22] = size(testmat);
% 
% for n2 = 1:n22
%     tempvec = testmat(:,n2);
%     if find(tempvec)
%         temp_idx2 = find(tempvec);
%         for qi = 1: length(temp_idx2)
%             idxpost1 = find(input2==temp_idx2(qi));  
%             idxpre1 = find(input1==n2);
%             res_clustering_temp(idxpost1) = n2;    
%              res_clustering_temp(idxpre1) = n2;           
%         end
%     end
% end
% 
% % figure, imagesc(input2), colormap colorcube
% % title('post')
% % 
% % figure, imagesc(input1), colormap colorcube
% % title('pre')
% % 
% % figure, imagesc(res_clustering_temp), colormap colorcube
% % title('synapse')
% % 
%     labeln = 1;
% res_cluster_synapse = zeros(size(res_clustering_temp));
% for n3 = 1: max(res_clustering_temp(:))
%     
%   if find(res_clustering_temp==n3)
%       
%     indk = find(res_clustering_temp ==n3);
%     
%     res_cluster_synapse(indk) =labeln;
%     labeln = labeln+1;
%   else
%       
%   end
%     
% end
% 
% 
% % %% confirm
% % 
% % test6 = zeros(size(res_cluster_synapse));
% % 
% % for  ni= 1:max(res_cluster_synapse(:))
% %     test6 = zeros(size(res_cluster_synapse));
% %     
% %     idx5 =  find(res_cluster_synapse==ni);
% %     
% %     test6(idx5) = ni;
% %     figure(123), imagesc(test6), colormap colorcube
% %     name = sprintf('Lable number is %d', ni);
% %     title(name)
% %     % title(['Acquisition number: ',num2str(ni)])
% %     pause(1)
% %     
% %  
% % end
% end
