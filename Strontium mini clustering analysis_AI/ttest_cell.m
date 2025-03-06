function [h,p] = ttest_cell(group1, group2, MEAN1, MEAN2, SEM1, SEM2, K)
x = 1:K;
for ii = 1:K
    [h{ii},p{ii},~,~] = ttest2(group1(ii,:), group2(ii,:));
end
figure;
hold on;
errorbar(x,MEAN1,SEM1, "-s","MarkerSize",5,...
    "MarkerEdgeColor","blue","MarkerFaceColor",[0.65 0.85 0.90])
errorbar(x,MEAN2,SEM2, "-s","MarkerSize",5,...
    "MarkerEdgeColor","red","MarkerFaceColor",[0.65 0.85 0.90])

for ii = 1:K
    if h{ii} == 1
        plot(ii,max((MEAN1(ii)+SEM1(ii)),(MEAN2(ii)+SEM2(ii)))*1.3,'*',"MarkerSize",3, "Color", 'r')
    end
end
hold off;
ticks = 1:K;
xticks(ticks)
xlim([0, K+1])
title("Cluster distribution")
xlabel("Cluster")
ylabel("Frequency (/s)")
legend("Group 1", "Group 2")

end