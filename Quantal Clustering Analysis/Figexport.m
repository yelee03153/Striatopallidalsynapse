
a = cellfun(@cell2mat,h9,'UniformOutput',false);
a = cellfun(@(x) nnz(x==1),a);
num(9) = nnz(a)

%%
figure;
scatter(transformedData_Excitatory(:,1), transformedData_Excitatory(:,2),'k.');
xlabel('x1');
ylabel('x2');
hold on;

% xidx1 = idx1{1}(1:2:length(idx1{1}));
% yidx1 = idx1{1}(2:2:length(idx1{1}));
% xidx2 = idx2{1}(1:2:length(idx2{1}));
% yidx2 = idx2{1}(2:2:length(idx2{1}));
% xidx3 = idx3{1}(1:2:length(idx3{1}));
% yidx3 = idx3{1}(2:2:length(idx3{1}));
% xidx4 = idx4{1}(1:2:length(idx4{1}));
% yidx4 = idx4{1}(2:2:length(idx4{1}));
% xidx42 = idx4{2}(1:2:length(idx4{2}));
% yidx42 = idx4{2}(2:2:length(idx4{2}));
% xidx5 = idx5{1}(1:2:length(idx5{1}));
% yidx5 = idx5{1}(2:2:length(idx5{1}));
% xidx6 = idx6{1}(1:2:length(idx6{1}));
% yidx6 = idx6{1}(2:2:length(idx6{1}));
% scatter(xidx4, yidx4, 'filled');
% scatter(xidx42, yidx42, 'filled');

for ii = 1:1000
    xidx = idx6{ii}(1:2:length(idx6{ii}));
    yidx = idx6{ii}(2:2:length(idx6{ii}));
    scatter(xidx, yidx, 'x');
    
end

%%
figure
hmap = [[0 num(1:3)]; [0 0 num(5) num(4)]; [0 0 0 num(6)]; [0 0 0 0]];
hmap = hmap.'+ hmap;
h2 = heatmap(hmap/1000*100);

XLabels = 1:4;
% Convert each number in the array into a string
CustomXLabels = categorical({'DL','VL','DM','VM'});
h2.XDisplayLabels = CustomXLabels;
h2.YDisplayLabels = CustomXLabels;
%%
figure; 
h = heatmap(num(7:9)/1000*100);
XLabels = 1:3;

CustomXLabels = categorical({'Diagonal','Lateral-Medial','Dorsal-Ventral'});
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = '';
