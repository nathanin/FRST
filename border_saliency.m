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