function [cluster1, cluster2, cluster3, cluster4] = SEPARATION(clusterIndices, labelmap,K)

cluster1 = clusterIndices(labelmap==1);
cluster2 = clusterIndices(labelmap==2);
cluster3 = clusterIndices(labelmap==3);
cluster4 = clusterIndices(labelmap==4);

clusterdist1 = [];
clusterdist2 = [];
clusterdist3 = [];
clusterdist4 = [];
for ii = 1:K
    clusterdist1(ii,1) = nnz(cluster1==ii);
    clusterdist2(ii,1) = nnz(cluster2==ii);
    clusterdist3(ii,1) = nnz(cluster3==ii);
    clusterdist4(ii,1) = nnz(cluster4==ii);
end

x = linspace(1,K,K);
figure;
hold on;
plot(x, clusterdist1,'Color', 'k', 'Marker','o','MarkerSize',3);
plot(x, clusterdist2,'Color', 'r', 'Marker','o','MarkerSize',3);
plot(x, clusterdist3,'Color', 'g', 'Marker','o','MarkerSize',3);
plot(x, clusterdist4,'Color', 'b', 'Marker','o','MarkerSize',3);
hold off;
ticks = 1:K;
xticks(ticks)
xlim([0 K+1])
title("Cluster distribution")
xlabel("Cluster")
ylabel("Event number")
legend("DL", "VL", "DM", "VM")
end
