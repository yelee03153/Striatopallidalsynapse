function union = unionfn(gfp,rfp)
temp = imbinarize(gfp);
temp2 = imbinarize(rfp);
temp(temp2==1) =1;
union = bwlabeln(temp);
