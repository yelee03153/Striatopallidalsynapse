
function [res_cluster_synapse] = extsynapseyelee (input1,input2)

%input1 : presynapse
%input2: postsynapse
%d_criterion2: criterion for synapse (distance between pre & post)

res_pre = imbinarize(input1, 0);

res_post = imbinarize(input2, 0);

res_synapse = colocalfn (res_pre,res_post);

res_cluster_synapse = bwlabeln(res_synapse);


end
