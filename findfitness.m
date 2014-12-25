function f = findfitness(S,mask,H)
% find a linear combination of fitness parameters:
%  Solidity, area, eccentricity, 'border saliency'
% such that higher sum of parameters = more likely to be a nucleus
% Can weight these to favor 'stronger' factors
    sol = S.Solidity;
    l = border_saliency(mask,H)/255;
    e = (0.754 - S.Eccentricity)*10; %difference from a determined average
    a = (S.Area - 145)/60; %adjusted to the scale of the other parameters..
    f = abs(sol+l+e+a);
return