function labels = recursiveKMeansAuto(X,kMax,sizeTh,delta)
    N = size(X,1);
    labels = zeros(N,1);
    nextID = 1;
    
    function split(idxs)
        nCur = numel(idxs);
        if nCur <= sizeTh
            labels(idxs) = nextID;
            nextID = nextID + 1;
            return;
        end
        
        Xc = X(idxs,:);
        kList = 2:min(kMax,nCur-1);
        Ec = evalclusters(Xc,'kmeans','CalinskiHarabasz','KList',kList);
        kOpt = Ec.OptimalK;
        base = Ec.CriterionValues(1);
        imp  = (Ec.CriterionValues(kOpt-1)-base)/base;

        if (nCur > 1000 && imp > delta)
            subIdx = kmeans(Xc,kOpt,'Replicates',10,'Display','off');
            for s = 1:kOpt
                split(idxs(subIdx==s));
            end
        else
            labels(idxs) = nextID;
            nextID = nextID + 1;
        end
    end

    split(1:N);
end
