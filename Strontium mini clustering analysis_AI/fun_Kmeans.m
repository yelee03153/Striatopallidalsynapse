function [K, clusterIndices] = fun_Kmeans(merged_data, dim, lim)
fh4 = @(X,K)(kmeans(X,K,'MaxIter',500));
eva = evalclusters(merged_data,fh4,"CalinskiHarabasz","KList",lim);
clear fh4
K = eva.OptimalK;
clusterIndices = eva.OptimalY;

% 군집 평가 기준 값 표시하기
figure
bar(eva.InspectedK,eva.CriterionValues);
xticks(eva.InspectedK);
xlabel("Number of clusters");
ylabel("Criterion values - 칼린스키-하라바츠(Calinski-Harabasz)");
legend("Optimal number of clusters is " + num2str(K))
title("Evaluation of Optimal Number of Clusters")
disp("Optimal number of clusters is " + num2str(K));

figure;
if dim == 2
    scatter(merged_data(:,1), merged_data(:,2),'k.'); 
    legend('data')
    xlabel('x1'); 
    ylabel('x2');
    hold on;
else
    scatter3(merged_data(:,1), merged_data(:,2), merged_data(:,3),'k.'); 
    legend('data')
    xlabel('x1'); 
    ylabel('x2');
    zlabel('x3');
    hold on;
end

% Labeling the clusters
if dim == 2
    for i = 1:K
        scatter(merged_data(clusterIndices ==i,1), merged_data(clusterIndices ==i,2), 'filled'); 
    end
else
    for i = 1:K
        scatter3(merged_data(clusterIndices ==i,1), merged_data(clusterIndices ==i,2), merged_data(clusterIndices ==i,3), 'filled'); 
    end
end
end
