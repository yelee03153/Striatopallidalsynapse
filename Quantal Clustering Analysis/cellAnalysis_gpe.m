function [clusterdist_mat, MEAN, SEM] = cellAnalysis_gpe(condition, cell_label, cluster, clusternum)


cell_ = cell(length(condition),1);
for ii = 1:length(condition)
    cell_{ii} = cluster(cell_label == ii);
end
%%
clusterdist = cell(1,length(condition));
for ii = 1:length(condition)
    cell_cluster = [];
    for jj = 1:clusternum
        cell_cluster(jj,1) = nnz(cell_{ii}==jj);
    end
    clusterdist{ii} = cell_cluster;
end

clusterdist_mat = cell2mat(clusterdist);
MEAN = mean(clusterdist_mat,2);
STD = std(clusterdist_mat,0,2);
SEM = STD/sqrt(size(clusterdist_mat,2));
end
