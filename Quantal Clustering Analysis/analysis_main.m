%% -------------------------------------------------------------------------
%   Cluster‑Region Distinctiveness Analysis
% -------------------------------------------------------------------------
%   • Load four brain regions (DL, VL, DM, VM)
%   • Clean 3‑D features (AMP, RT, DT) and z‑standardize
%   • Recursive k‑means sub‑clustering (Calinski‑Harabasz split rule)
%   • χ² (df = 3) test per cluster → Region‑by‑cluster Z‑matrix
%   • Region‑wise ΣZ² + pairwise ΔS permutation tests (Holm/FDR)
% -------------------------------------------------------------------------
%   Author : Minseok Jeong
%   Date   : 2025-06-13
% -------------------------------------------------------------------------

clear; clc; close all;

%% 0 | Load ----------------------------------------------------------------
load('GPe_SR2+_mini.mat');      % RESULT.DL / VL / DM / VM

%% 1 | Concatenate regions -------------------------------------------------
regions = {'DL','VL','DM','VM'};
X  = [];                         % raw feature matrix
RL = [];                         % region label (1–4)
for r = 1:numel(regions)
    condNames = fields(RESULT.(regions{r}));
    for k = 1:numel(condNames)
        T  = RESULT.(regions{r}).(condNames{k});
        X  = [X; [T.AMP, T.RT, T.DT]];
        RL = [RL; r*ones(height(T),1)];
    end
end

%% 2 | Cleaning ------------------------------------------------------------
X   = X(~any(isnan(X),2),:);       % drop rows with NaN
RL  = RL(~any(isnan(X),2));
flag = isoutlier(X,'percentiles',[0.01 99.99]);
X   = X(~any(flag,2),:);
RL  = RL(~any(flag,2));
fprintf('Rows kept after cleaning: %d\n',size(X,1));

%% 3 | Z‑score -------------------------------------------------------------
[Xz,mu,sigma] = zscore(X);

%% 3' | PCA ------------------------------------------------------------
[coeff, score, latent, tsquared, explained] = pca(X);

[coeff_z, score_z, latent_z, tsquared_z, explained_z] = pca(Xz);

% --- Explained Variance Bar Graph (Scree Plot) ---
figure;
bar(explained);
hold on;
cumulativeExplained = cumsum(explained);
plot(1:numel(explained), cumulativeExplained, ':o', 'LineWidth', 1, 'Color', 'r');
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot: Explained Variance by Principal Component (Raw data)');
grid on;
legend('Individual Variance', 'Cumulative Variance');

figure;
bar(explained_z);
hold on;
cumulativeExplained_z = cumsum(explained_z);
plot(1:numel(explained_z), cumulativeExplained_z, ':o', 'LineWidth', 1, 'Color', 'r');
xlabel('Principal Component');
ylabel('Variance Explained (%)');
title('Scree Plot: Explained Variance by Principal Component (Z-scored data)');
grid on;
legend('Individual Variance', 'Cumulative Variance');

% --- Extract first two principal components ---
X_2d = score(:,1:2);
% --- Visualization ---


figure;
scatter(X_2d(:,1), X_2d(:,2), 5, 'filled');
xlabel('Principal Component 1');
ylabel('Principal Component 2');
title('PCA 2D Projection');
grid on;

%% 4 | Recursive k‑means ---------------------------------------------------
finalLab = recursiveKMeansAuto(Xz,10,500,0.05);
fprintf('Final #clusters = %d\n',numel(unique(finalLab)));
tabulate(finalLab)

% Extended Data Fig.8a
figure; scatter3(Xz(:,1),Xz(:,2),Xz(:,3),6,finalLab,'filled');
title('Final recursive k‑means clusters');

%% 5 | χ² test per cluster -------------------------------------------------
C = numel(unique(finalLab)); R = 4;
N = numel(RL);
p_r = accumarray(RL,1,[R,1])/N;
p_min = min(p_r);
minCell = 20; 
minClustSize = ceil(minCell / p_min);
chiObs = NaN(C,1); Zmat = NaN(C,R);
for c = 1:C
    idx = finalLab==c; if sum(idx)<minClustSize, continue, end
    tbl = [accumarray(RL(idx),1,[R,1],@sum,0), accumarray(RL(~idx),1,[R,1],@sum,0)];
    if any(tbl(:)<minCell), continue, end
    expT = sum(tbl,2).*sum(tbl,1)/sum(tbl,'all');
    chiObs(c)=sum((tbl-expT).^2./(expT+eps),'all');
    Zmat(c,:)=(tbl(:,1)-expT(:,1))./sqrt(expT(:,1)+eps);
end

%% 6 | Region‑wise ΣZ² + permutation --------------------------------------
valid = ~isnan(Zmat(:,1));
Sobs  = nansum(Zmat(valid,:).^2,1);
B = 1e4; rng(4); n = numel(RL);
Sperm = zeros(B,R);
for b = 1:B
    perm = RL(randperm(n));
    tbl  = accumarray([finalLab,perm],1,[C,R]);
    expT = sum(tbl,2).*sum(tbl,1)/sum(tbl,'all');
    Z    = (tbl-expT)./sqrt(expT+eps);
    Sperm(b,:) = sum(Z(valid,:).^2,1);
end
pPerm = (sum(Sperm>=Sobs,1)+1)/(B+1);

fprintf('\nRegion   S_obs   p_perm\n');
for r = 1:R, fprintf('%6d  %7.1f   %.4f\n',r,Sobs(r),pPerm(r)); end



%% 7 | Pairwise ΔS permutation + Holm / BH --------------------------------
Pairs   = nchoosek(1:R,2);  m = size(Pairs,1); rawP = zeros(m,1);
for i = 1:m
    a=Pairs(i,1); b=Pairs(i,2);
    rawP(i) = (sum(abs(Sperm(:,a)-Sperm(:,b))>=abs(Sobs(a)-Sobs(b)))+1)/(B+1);
end
[pSrt,ix] = sort(rawP); holm=min(1,cummax(pSrt.*(m:-1:1)')); pHolm=zeros(m,1); pHolm(ix)=holm;
qBH = mafdr(rawP,'BHFDR',true);
Tpairs = table(Pairs(:,1),Pairs(:,2),rawP,pHolm,qBH,'VariableNames',{'R_A','R_B','pRaw','pHolm','qBH'});
disp(Tpairs);

%% 8 | Visualisations ------------------------------------------------------
% 8‑a ΔS heat‑map (row>col = red)
Delta = NaN(R); pMat=NaN(R);
for i=1:m
    a=Pairs(i,1); b=Pairs(i,2); d=Sobs(a)-Sobs(b);
    Delta(a,b)= d; Delta(b,a)=-d; pMat(a,b)=pHolm(i); pMat(b,a)=pHolm(i);
end
Delta(1:R+1:end)=0;

% Extended Data Fig.8c
figure('Name','DeltaS heat‑map'); imagesc(Delta); axis square; colormap(redbluecmap);
cb=colorbar; cb.Label.String='ΔS (S_r - S_s)';
set(gca,'XTick',1:R,'XTickLabel',{'DL','VL','DM','VM'},'YTick',1:R,'YTickLabel',{'DL','VL','DM','VM'});
title('ΔS between Regions (red: row greater)'); hold on;
[row,col]=find(pMat<0.05); plot(col,row,'k*','MarkerSize',8,'LineWidth',1.2);
rectangle('Position',[0.5 3.5 R 1],'EdgeColor','k','LineWidth',1.4);
rectangle('Position',[3.5 0.5 1 R],'EdgeColor','k','LineWidth',1.4);

% 8‑b Normalised ΔS bar (Region‑4 focus)
% Extended Data Fig.8b
figure('Name','VM ΔS bar');
bar(abs(Delta(4,[1 2 3]))/Sobs(4),'FaceColor',[0.4 0.6 1]);
set(gca,'XTick',1:3,'XTickLabel',{'VM‑DL','VM‑VL','VM‑DM'});
ylim([0, 1.0])
ylabel('|ΔS| / S_V_M'); title('Distinctiveness of VM vs others');

% 8‑c Region ΣZ² horizontal bar
% Extended Data Fig.8d
figure('Name','ΣZ² bar'); [Ssort,ord]=sort(Sobs,'descend');
barh(Ssort,'FaceColor',[.7 .7 .7]); hold on; barh(find(ord==4),Ssort(ord==4),'r');
set(gca,'YTick',1:R,'YTickLabel',compose('R%d',ord)); xlabel('ΣZ²'); title('Region‑wise distinctiveness');
text(Ssort+5,1:R,compose('p=%.4f',pPerm(ord)),'FontSize',8);

% 8‑d Permutation cloud vs observed ΣZ²
% Extended Data Fig.8e
figure('Name','Permutation cloud'); hold on;
for r=1:R, scatter(r*ones(B,1),Sperm(:,r),4,[.8 .8 .8],'filled'); end
plot(1:R,Sobs,'r*','MarkerSize',10,'LineWidth',1.1);
xlim([0.5 R+0.5]); ylabel('ΣZ²'); set(gca,'XTick',1:R,'XTickLabel',{'DL','VL','DM','VM'});
title('Permutation distribution of ΣZ² (red star = observed)');

%% 9 | Save figures and data ------------------------------------------------
figHandles = findall(groot, 'Type', 'Figure');
if isempty(figHandles)
    fprintf('No figures to save.\n');
else
    for i = 1:length(figHandles)
        figTitle = figHandles(i).CurrentAxes.Title.String;
        if ~isempty(figTitle)
            % Replace spaces with underscores for the filename
            filename = strrep(figTitle, ' ', '_');
            filename = strrep(filename, '(', '');
            filename = strrep(filename, ')', '');
            filename = strrep(filename, ':', '');
            % Save the figure as a PDF with 600 DPI
            print(figHandles(i), filename, '-dpdf', '-r600');
        end
    end
end

% Save resulting data
save('analysis_results.mat', 'Xz', 'Zmat', 'Sobs', 'Sperm', 'Delta', 'pPerm','pMat');

%% -------------------------------------------------------------------------
function [p,chi2,df] = chi2gof2D(tbl)
    expT=sum(tbl,2).*sum(tbl,1)/sum(tbl,'all');
    chi2=sum((tbl-expT).^2 ./ (expT+eps),'all');
    df=(size(tbl,1)-1)*(size(tbl,2)-1); p=1-chi2cdf(chi2,df);
end

