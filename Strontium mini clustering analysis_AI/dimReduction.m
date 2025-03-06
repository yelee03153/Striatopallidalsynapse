function [transformedData] = dimReduction(Merged_data)
%% Dimensionality reduction

% 주성분 분석 - Excitatory
[coeffs, transformedData, ~, ~, explained] = pca(Merged_data);

% 결과 표시하기
figure
colorMap = colormap("lines");
barPlot = bar(explained);
hold on;
plot(cumsum(explained),"*-","Color",colorMap(3,:));
hold off;
barPlot(1).FaceColor = "flat";
barPlot(1).CData(1:2,:) = repmat(colorMap(2,:),[2 1]);
label = ["Selected components"; sprintf("Cumulative explained\nvariance")];
legend(label,"Location","best");
title("Scree Plot of Explained Variances")
xlabel("Principal component")
ylabel("Variance explained (%)")
clear barPlot colorMap label explained

figure
scatter(transformedData(:,1),transformedData(:,2), "Marker", ".");
title("Scatter Plot of Principal Components")
xlabel("Component 1")
ylabel("Component 2")

figure
for i = 1:size(coeffs, 2)
    varLabels{i} = ['Var', num2str(i)];
end
biplot(coeffs(:,[1 2]),"VarLabels",varLabels)
title("PCA Biplot")
xlabel("Component 1")
ylabel("Component 2")
clear coeffs explained i varLabels
transformedData = transformedData(:,1:2);
end