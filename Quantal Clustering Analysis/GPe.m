clear all
clc
close all
load('GPe_SR2+_mini_4.mat')
%% DL

condition1 = fields(RESULT.DL);
DL_amp = [];
DL_rt = [];
DL_dt = [];
DL_label = zeros(1,length(condition1));
for ii = 1:length(condition1)
    DL_amp = [DL_amp; RESULT.DL.(condition1{ii}).AMP];
    DL_rt = [DL_rt; RESULT.DL.(condition1{ii}).RT];
    DL_dt = [DL_dt; RESULT.DL.(condition1{ii}).DT];
    DL_label(ii) = height(RESULT.DL.(condition1{ii}));
end
DL_meged = [DL_amp DL_rt DL_dt];

%% VL

condition2 = fields(RESULT.VL);
VL_amp = [];
VL_rt = [];
VL_dt = [];
VL_label = zeros(1,length(condition2));
for ii = 1:length(condition2)
    VL_amp = [VL_amp; RESULT.VL.(condition2{ii}).AMP];
    VL_rt = [VL_rt; RESULT.VL.(condition2{ii}).RT];
    VL_dt = [VL_dt; RESULT.VL.(condition2{ii}).DT];
    VL_label(ii) = height(RESULT.VL.(condition2{ii}));
end
VL_meged = [VL_amp VL_rt VL_dt];

%% DM

condition3 = fields(RESULT.DM);
DM_amp = [];
DM_rt = [];
DM_dt = [];
DM_label = zeros(1,length(condition3));
for ii = 1:length(condition3)
    DM_amp = [DM_amp; RESULT.DM.(condition3{ii}).AMP];
    DM_rt = [DM_rt; RESULT.DM.(condition3{ii}).RT];
    DM_dt = [DM_dt; RESULT.DM.(condition3{ii}).DT];
    DM_label(ii) = height(RESULT.DM.(condition3{ii}));
end

DM_meged = [DM_amp DM_rt DM_dt];

%% VM

condition4 = fields(RESULT.VM);
VM_amp = [];
VM_rt = [];
VM_dt = [];
VM_label = zeros(1,length(condition4));
for ii = 1:length(condition4)
    VM_amp = [VM_amp; RESULT.VM.(condition4{ii}).AMP];
    VM_rt = [VM_rt; RESULT.VM.(condition4{ii}).RT];
    VM_dt = [VM_dt; RESULT.VM.(condition4{ii}).DT];
    VM_label(ii) = height(RESULT.VM.(condition4{ii}));
end
VM_meged = [VM_amp VM_rt VM_dt];

%% Merge

Allmerged = [DL_meged; VL_meged; DM_meged; VM_meged];

%% Label egion and cell #

region_label = [ones(1,length(DL_meged)), 2.*ones(1,length(VL_meged)), 3.*ones(1,length(DM_meged)), 4.*ones(1,length(VM_meged))];

cell_label_DL = [];
for ii = 1:length(DL_label)
    cell_label_DL = [cell_label_DL; repmat([ii], [DL_label(ii), 1])];
end
cell_label_VL = [];
for ii = 1:length(VL_label)
    cell_label_VL = [cell_label_VL; repmat([ii], [VL_label(ii), 1])];
end
cell_label_DM = [];
for ii = 1:length(DM_label)
    cell_label_DM = [cell_label_DM; repmat([ii], [DM_label(ii), 1])];
end
cell_label_VM = [];
for ii = 1:length(VM_label)
    cell_label_VM = [cell_label_VM; repmat([ii], [VM_label(ii), 1])];
end
%% Dimensionality reduction

% PCA
transformedData_Excitatory = dimReduction(Allmerged);

%% 2D kmeans
cutoff = 50;
rep = 1:1000;
idx1 = cell(1,length(rep));
idx2 = cell(1,length(rep));
idx3 = cell(1,length(rep));
idx4 = cell(1,length(rep));
idx5 = cell(1,length(rep));
idx6 = cell(1,length(rep));
idx7 = cell(1,length(rep));
idx8 = cell(1,length(rep));
idx9 = cell(1,length(rep));
for ii = rep
    Krange = 50;
    [K_2D_Excitatory, clusterIndices_2D_Excitatory] = fun_Kmeans(transformedData_Excitatory, 2, Krange);

    %% Region separation
    [DL_cluster,VL_cluster,DM_cluster,VM_cluster] = SEPARATION(clusterIndices_2D_Excitatory,region_label,K_2D_Excitatory);

    %% Cell comparison
    [cluster_num_DL, MEAN_DL, SEM_DL] = cellAnalysis_gpe(condition1, cell_label_DL, DL_cluster, K_2D_Excitatory);
    [cluster_num_VL, MEAN_VL, SEM_VL] = cellAnalysis_gpe(condition2, cell_label_VL, VL_cluster, K_2D_Excitatory);
    [cluster_num_DM, MEAN_DM, SEM_DM] = cellAnalysis_gpe(condition3, cell_label_DM, DM_cluster, K_2D_Excitatory);
    [cluster_num_VM, MEAN_VM, SEM_VM] = cellAnalysis_gpe(condition4, cell_label_VM, VM_cluster, K_2D_Excitatory);
    sum_DL = sum(cluster_num_DL,2);
    sum_VL = sum(cluster_num_VL,2);
    sum_DM = sum(cluster_num_DM,2);
    sum_VM = sum(cluster_num_VM,2);
    %% ttest
    [h1{ii},p1{ii}] = ttest_cell(cluster_num_DL, cluster_num_VL, MEAN_DL, MEAN_VL, ...
        SEM_DL, SEM_VL, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h1{ii}) == 1))
        for i = find(cell2mat(h1{ii})==1)
            if (sum_DL(i) + sum_VL(i)) > cutoff
                idx1{ii} = [idx1{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h2{ii},p2{ii}] = ttest_cell(cluster_num_DL, cluster_num_DM, MEAN_DL, MEAN_DM, ...
        SEM_DL, SEM_DM, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h2{ii}) == 1))
        for i = find(cell2mat(h2{ii})==1)
            if (sum_DL(i) + sum_DM(i)) > cutoff
                idx2{ii} = [idx2{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h3{ii},p3{ii}] = ttest_cell(cluster_num_DL, cluster_num_VM, MEAN_DL, MEAN_VM, ...
        SEM_DL, SEM_VM, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h3{ii}) == 1))
        for i = find(cell2mat(h3{ii})==1)
            if (sum_DL(i) + sum_VM(i)) > cutoff
                idx3{ii} = [idx3{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h4{ii},p4{ii}] = ttest_cell(cluster_num_VL, cluster_num_VM, MEAN_VL, MEAN_VM, ...
        SEM_VL, SEM_VM, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h4{ii}) == 1))
        for i = find(cell2mat(h4{ii})==1)
            if (sum_VL(i) + sum_VM(i)) > cutoff
                idx4{ii} = [idx4{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h5{ii},p5{ii}] = ttest_cell(cluster_num_VL, cluster_num_DM, MEAN_VL, MEAN_DM, ...
        SEM_VL, SEM_DM, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h5{ii}) == 1))
        for i = find(cell2mat(h5{ii})==1)
            if (sum_VL(i) + sum_DM(i)) > cutoff
                idx5{ii} = [idx5{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h6{ii},p6{ii}] = ttest_cell(cluster_num_DM, cluster_num_VM, MEAN_DM, MEAN_VM, ...
        SEM_DM, SEM_VM, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h6{ii}) == 1))
        for i = find(cell2mat(h6{ii})==1)
            if (sum_DM(i) + sum_VM(i)) > cutoff
                idx6{ii} = [idx6{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h7{ii},p7{ii}] = ttest_cell([cluster_num_DL cluster_num_VM], [cluster_num_VL cluster_num_DM], MEAN_DL+MEAN_VM, MEAN_VL+MEAN_DM, ...
        SEM_DL+SEM_VM,SEM_DM+SEM_VL, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h7{ii}) == 1))
        for i = find(cell2mat(h7{ii})==1)
            if (sum_DL(i) + sum_VL(i) + sum_DM(i) + sum_VM(i)) > cutoff*2
                idx7{ii} = [idx7{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h8{ii},p8{ii}] = ttest_cell([cluster_num_DL cluster_num_VL], [cluster_num_DM cluster_num_VM], MEAN_DL+MEAN_VL, MEAN_DM+MEAN_VM, ...
        SEM_DL+SEM_VL, SEM_DM+SEM_VM, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h8{ii}) == 1))
        for i = find(cell2mat(h8{ii})==1)
            if (sum_DL(i) + sum_VL(i) + sum_DM(i) + sum_VM(i)) > cutoff*2
                idx8{ii} = [idx8{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    %%
    [h9{ii},p9{ii}] = ttest_cell([cluster_num_DL cluster_num_DM], [cluster_num_VL cluster_num_VM], MEAN_DL+MEAN_DM, MEAN_VL+MEAN_VM, ...
        SEM_DL+SEM_DM, SEM_VL+SEM_VM, K_2D_Excitatory);
    if ~isempty(find(cell2mat(h3{ii}) == 1))
        for i = find(cell2mat(h3{ii})==1)
            if (sum_DL(i) + sum_VL(i) + sum_DM(i) + sum_VM(i)) > cutoff*2
                idx9{ii} = [idx9{ii} mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1)) mean(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2))];
            end
        end
    end
    close all
end
%%
idx_tot = {idx1, idx2, idx3, idx4, idx5, idx6, idx7, idx8, idx9};

%%
for ii = 1:length(idx_tot)
    a = cellfun(@length,idx_tot{ii});
    num(ii) = nnz(a);
end
%%
figure;
scatter(transformedData_Excitatory(:,1), transformedData_Excitatory(:,2),'k.');
xlabel('x1');
ylabel('x2');
hold on;

for ii = find(num>0)
    xidx = [];
    yidx = [];
    for jj = rep
        xidx = [xidx idx_tot{ii}{jj}(1:2:length(idx_tot{ii}{jj}))];
        yidx = [yidx idx_tot{ii}{jj}(2:2:length(idx_tot{ii}{jj}))];
    end
    scatter(xidx, yidx, '+');
end
labels = {'', 'DL vs VL', 'DL vs DM', 'DL vs VM', 'VL vs DM', 'VL vs VM', 'DM vs VM', 'Diagonal', 'LM', 'DV'};
legend(labels([1 find(num>0)+1]))
%%
figure
hmap = [[0 num(1:3)]; [0 0 num(5) num(4)]; [0 0 0 num(6)]; [0 0 0 0]];
hmap = hmap.'+ hmap;
h2 = heatmap(hmap/length(rep)*100);

XLabels = 1:4;
% Convert each number in the array into a string
CustomXLabels = categorical({'DL','VL','DM','VM'});
h2.XDisplayLabels = CustomXLabels;
h2.YDisplayLabels = CustomXLabels;
%%
figure; 
h = heatmap(num(7:9)/length(rep)*100);
XLabels = 1:3;

CustomXLabels = categorical({'Diagonal','Lateral-Medial','Dorsal-Ventral'});
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = '';

%{
%% ttest
[h1,p1] = ttest_cell(cluster_num_DL, cluster_num_VL, MEAN_DL, MEAN_VL, ...
    SEM_DL, SEM_VL, K_2D_Excitatory);
%%
[h2,p2] = ttest_cell(cluster_num_DL, cluster_num_DM, MEAN_DL, MEAN_DM, ...
    SEM_DL, SEM_DM, K_2D_Excitatory);
%%
[h3,p3] = ttest_cell(cluster_num_VL, cluster_num_VM, MEAN_VL, MEAN_VM, ...
    SEM_VL, SEM_VM, K_2D_Excitatory);
%%
[h8,p8] = ttest_cell(cluster_num_DM, cluster_num_VM, MEAN_DM, MEAN_VM, ...
    SEM_DM, SEM_VM, K_2D_Excitatory);
%%

[h4,p4] = ttest_cell([cluster_num_DL cluster_num_VM], [cluster_num_VL cluster_num_DM], MEAN_DL+MEAN_VM, MEAN_VL+MEAN_DM, ...
    SEM_DL+SEM_VM,SEM_DM+SEM_VL, K_2D_Excitatory);
%%
[h5,p5] = ttest_cell(cluster_num_DL, cluster_num_VM, MEAN_DL, MEAN_VM, ...
    SEM_DL, SEM_VM, K_2D_Excitatory);
%%
[h6,p6] = ttest_cell(cluster_num_VL, cluster_num_DM, MEAN_VL, MEAN_DM, ...
    SEM_VL, SEM_DM, K_2D_Excitatory);
%%
[h7,p7] = ttest_cell([cluster_num_DL cluster_num_VL], [cluster_num_DM cluster_num_VM], MEAN_DL+MEAN_VL, MEAN_DM+MEAN_VM, ...
    SEM_DL+SEM_VL, SEM_DM+SEM_VM, K_2D_Excitatory);
%%
[h9,p9] = ttest_cell([cluster_num_DL cluster_num_DM], [cluster_num_VL cluster_num_VM], MEAN_DL+MEAN_DM, MEAN_VL+MEAN_VM, ...
    SEM_DL+SEM_DM, SEM_VL+SEM_VM, K_2D_Excitatory);

%%
figure;
scatter(transformedData_Excitatory(:,1), transformedData_Excitatory(:,2),'k.');
xlabel('x1');
ylabel('x2');
hold on;

for i = find(cell2mat(h8)==1)
    scatter(transformedData_Excitatory(clusterIndices_2D_Excitatory==i,1), transformedData_Excitatory(clusterIndices_2D_Excitatory==i,2), 'filled');
end

%%

regions_number = zeros(4,K_2D_Excitatory);
for ii = 1:4
    for jj = 1:K_2D_Excitatory
        regions_number(ii,jj) = length(transformedData_Excitatory((region_label==ii) .* clusterIndices_2D_Excitatory' == jj,1));
    end
end

summed_number = sum(regions_number,1)
regions_number(:,summed_number < length(transformedData_Excitatory)*0.001) = [];
regions_number(regions_number == 0) = exp(-30);
%%
score12 = getdiff(regions_number(1,:),regions_number(2,:))
score13 = getdiff(regions_number(1,:),regions_number(3,:))
score14 = getdiff(regions_number(1,:),regions_number(4,:))
score23 = getdiff(regions_number(2,:),regions_number(3,:))
score24 = getdiff(regions_number(2,:),regions_number(4,:))
score34 = getdiff(regions_number(3,:),regions_number(4,:))

scoreDiagonal = getdiff(regions_number(1,:)+regions_number(4,:),regions_number(2,:)+regions_number(3,:))
scoreDV = getdiff(regions_number(1,:)+regions_number(3,:),regions_number(2,:)+regions_number(4,:))
scoreML = getdiff(regions_number(1,:)+regions_number(2,:),regions_number(3,:)+regions_number(4,:))
    
mat = zeros(4,4);

mat(1,2) = score12;
mat(1,3) = score13;
mat(1,4) = score14;
mat(2,3) = score23;
mat(2,4) = score23;
mat(3,4) = score34;

mat = mat + mat.';


figure(), 
h = heatmap(mat);
XLabels = 1:4;
% Convert each number in the array into a string
CustomXLabels = categorical({'DL','VL','DM','VM'});
h.XDisplayLabels = CustomXLabels;
h.YDisplayLabels = CustomXLabels;
%}
%% Local functions

function [Amp, RT, DT] = getparameters(struct)
    
    Amp = [];
    RT = [];
    DT = [];
    for ii = 1:length(struct)
        Amp = [Amp struct{ii}.AMP'];
        RT = [RT struct{ii}.RT'];
        DT = [DT struct{ii}.DT'];
    end
    
end

function [Amp, RT, DT] = removeNan(orgAmp, orgRT, orgDT)
    nanidx = isnan(orgDT);
    Amp = orgAmp(~nanidx);
    RT = orgRT(~nanidx);
    DT = orgDT(~nanidx);
end

%%
function [score] = getdiff(a,b)
    diff = abs(a-b)./(a+b);
    score = sum(diff);
end
