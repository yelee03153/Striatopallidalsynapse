regions = {'DL', 'VL', 'DM', 'VM'};

for jj = 1:4
    region = regions{jj};
    T = readtable('Sr2+ GPe mini data_YELEE.xlsx','Sheet',region,'VariableNamingRule', 'preserve');
    nanidx = find(isnan(T.Var1));
    newidx = zeros(1,length(nanidx)+1);
    newidx(1) = 0;
    newidx(2:length(nanidx)+1) = nanidx;
    newidx = [newidx height(T)+1];
    for ii = 1:length(newidx)-1
        AMP = T.Amplitude(newidx(ii)+1:newidx(ii+1)-1);
        RT = T.("Rise (ms)")(newidx(ii)+1:newidx(ii+1)-1);
        DT = T.("Decay (ms)")(newidx(ii)+1:newidx(ii+1)-1);
        tmpTable = table(AMP,RT,DT);
        RESULT.(region).(strcat('cell',num2str(ii))) = tmpTable;
    end
    
end
%%
save("GPe_SR2+_mini","RESULT");
