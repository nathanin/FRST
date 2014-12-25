function bcgmask = findbcg(HE)

grayHE = rgb2gray(HE);  %rgb2gray !!!!!! 
bcgmask = grayHE >= 210;  % ad hoc thresholding to determine background areas

[L, ~] = bwlabel(bcgmask);
stat = regionprops(L,'Area');
trh = 50;
indx = ([stat.Area] >= trh); 
bcgmask = ismember(L, find(indx));
bcgmask = imclose(bcgmask,strel('disk',7));  % cleaned up background mask 

return