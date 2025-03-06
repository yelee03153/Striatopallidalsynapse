function colocalbn = colocalfnbinarize(gfp,rfp)
temp = imbinarize(gfp);
temp2 = imbinarize(rfp);
temp(temp2~=1) =0;
colocalbn = bwlabeln(temp);

