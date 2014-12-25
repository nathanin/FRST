function [regions, LL, newL, STATout]=Post_process(Hdeconv, labels, r, STATS)
% 'Solidity', 'Area', 'Centroid', 'MajorAxisLength',...
%                                'MinorAxisLength', 'Orientation', 'WeightedCentroid',...
%                                'MeanIntensity', 'Eccentricity'


regions = []; nreg = 1;
newL = zeros(size(labels), 'uint16');
% STATout = struct('Area',[],'Centroid',[],'MajorAxisLength',[],'MinorAxisLength',[],...
%                  'Orientation',[],'Solidity',[],'WeightedCentroid',[]);
% masks = struct([]);
for j = 2:length(STATS)
    A = STATS(j).Area;
    S = STATS(j).Solidity;
    C = STATS(j).Centroid; Cx=C(1); Cy=C(2);
    WC = STATS(j).WeightedCentroid; WCx=WC(1); WCy=WC(2);
    d = sqrt((Cx-WCx)^2 + (Cy-WCy)^2);
    EC = STATS(j).Eccentricity;
%     I = STATS(j).MeanIntensity;
    
%     mask = labels==j;
%     l = border_saliency(mask, Hdeconv);
    
    if A<(r^2)*pi || A>((4*r)^2)*pi || S<0.815 || d<0.08 || EC>0.9
        labels(labels==j) = 0;
    else
        mask = labels==j;
        l = border_saliency(mask, Hdeconv);
        
        if l<20, continue;   end
        
        regions(nreg) = j;
        newL(mask) = nreg;
%         masks(nreg,1).mask = mask;
%         STATout(nreg) = STATS(j);
        nreg=nreg+1;
    end
end
STATout = regionprops(newL, Hdeconv, 'Solidity', 'Area', 'Centroid', 'MajorAxisLength',...
                      'MinorAxisLength', 'Orientation', 'WeightedCentroid', 'Eccentricity',...
                      'MeanIntensity');
LL = fit_ellipses(newL, STATout);

return

function l = border_saliency(mask, H)
%     mask=false(size(orim));
    invmask = ~mask;
    Bout = bwperim(invmask);
    Bout(1,:)=0; Bout(:,1)=0; Bout(end,:)=0; Bout(:,end)=0;
    Bin = bwperim(mask);
    Bin(1,:)=0; Bin(:,1)=0; Bin(end,:)=0; Bin(:,end)=0;
    Bout = bwperim(~Bout);
    Bout(1,:)=0; Bout(:,1)=0; Bout(end,:)=0; Bout(:,end)=0;
    Bout(Bin) = 0;

    inpx = H(Bin);
    outpx = H(Bout);
    mi = double(median(inpx));
    mo = double(median(outpx));
    l = abs(mi - mo);
return