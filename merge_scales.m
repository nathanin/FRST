function [res, accepted, Snuclei] = merge_scales(L, incr, reg, H, S)

% add accounting for the number of regions in each layer of L, so during
% the for-loop execution it can skip it's own layer of L saving at least
% some time.
%
% Different way, not using overlaps:
% find distance between the centroid of the region in question (outer loop)
% and the centroid of the comparisons. if the distance is too large, skip
% finding the overlap? - find a radius of second regions that might be
% overlapping and only look at those? the centroid lookup form S will be
% faster than finding overlaps..? Use the first region centroid as a center
% frame of reference for the masks & comparison?

A = false(reg, reg);
f = zeros(reg,1);
accepted = zeros(1, reg);
% Construct adjacency matrix
for i=1:reg
    [imask, layer] = mfromL(L,i);
    ni = nnz(imask);
    f(i,1)= findfitness(S(i,1),imask,H);
    Atemp = zeros(1,reg);
    for j=incr(layer)+1:reg
%     for j=i+1:reg
        [jmask, ~] = mfromL(L,j);
        nj = nnz(jmask);
        
        Over = sum(imask(:) & jmask(:)) / min(ni,nj);
        
        if Over>0.2, Atemp(1,j) = true;     end 
    end
    A(i,:) = Atemp;
end

% DEBUG/CHECK FOR ACCURACY
% Find the best regions from the overlaps
% add more parameters than just solidity
% -----------------------------------------------------------------
% 1. Find largest s, mark as accepted.
% 2. Find the row of A that corresponds to the index of largest s
% 3. Mark all indices in that row which are true as rejected
% 4. Repeat for the next largest value of s.
% -----------------------------------------------------------------

iter=0;
while nnz(accepted)~=length(accepted)
     lrg = max(f); 
     lrg = find(f==lrg);
     lrg = lrg(1);
     
     vec = A(lrg,:);
     accepted(vec==true) = 1; %reject
     accepted(lrg) = 2; %accept
     
%      f1=f(1:lrg-1); f2=f(lrg+1:end);
%      f=[f1; f2];
     f(lrg) = 0;
     f(vec==true)=0;
     iter=iter+1;
     if iter>reg, break;    end
end
% iter
accepted = accepted-1;

% Reconstruct the accepted regions into a single mask.
res=zeros(size(L,1), size(L,2), 'uint16');
% accepted = unique(accepted);
accepted=find(accepted==1);
for i=1:length(accepted)
    [mask, ~] = mfromL(L,accepted(i));
    res(mask) = i;
end

% Get stats from the accepted regions
Snuclei = regionprops(res, H, 'Solidity', 'Area', 'Eccentricity', 'MeanIntensity');
return

function [f]= findfitness(S,mask,H)
% find a linear combination of fitness parameters:
%  Solidity, area, eccentricity, 'border saliency', mean intensity
% such that higher sum of parameters = more likely to be a nucleus
% Can weight these to favor 'stronger' factors

    sol = S.Solidity; %high solidity preferred
    l = border_saliency(mask,H)/255; %high saliency preferred
    e = 1-S.Eccentricity; %more circular areas preferred
    a = (S.Area - 142)/105; %larger areas preferred
    m = (112 - S.MeanIntensity)/255; %darker areas preferred
    f = abs(sol+l+m+a+e);
    
return

function [mask, layer] = mfromL(L,n)

bin = L==n;
layer = squeeze(max(max(bin)))==1;
mask = L(:,:,layer)==n;

return