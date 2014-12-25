function LL = fit_ellipses(L, S)

% fit ellipses to a label matrix
% clear all; close all;
% load('/Users/nathaning/Documents/MATLAB/Cedars/FRST/nuclei_test1.mat');
% L = nuclei;
n=max(max(L));
% phi = linspace(0,2*pi,50);
% cosphi = cos(phi);
% sinphi = sin(phi);
% S = regionprops(L,'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Orientation');
LL = zeros(size(L));
for j=1:n
%     %``````````````` Parametric plot ellipse `````````````````````
    xbar = S(j).Centroid(1);
    ybar = S(j).Centroid(2);

    a = S(j).MajorAxisLength/2;
    b = S(j).MinorAxisLength/2;

    theta = pi*S(j).Orientation/180;
%     R = [ cos(theta)   sin(theta)
%          -sin(theta)   cos(theta)];
% 
%     xy = [a*cosphi; b*sinphi];
%     xy = R*xy;
% 
%     x = xy(1,:) + xbar;
%     y = xy(2,:) + ybar;
% 
%     ellipses(j).x = x;
%     ellipses(j).y = y;
%     ellipses(j).xbar = xbar;
%     ellipses(j).ybar = ybar;
 
%     %``````````````` mask `````````````````````
    % http://stackoverflow.com/questions/11079781/cropping-an-ellipse-from-an-image
    % http://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate
    
    [X,Y] = meshgrid(1:size(L,2), 1:size(L,1));
    ellipse_mask = ((X-xbar)*cos(theta) - (Y-ybar)*sin(theta)).^2/a^2 +...
                    ((X-xbar)*sin(theta) + (Y-ybar)*cos(theta)).^2/b^2 <= 1;
    
    % at this point, can check if there is already a similar region
    % segmented from another scale -- or accept/reject regions based on
    % some criteria
    LL(ellipse_mask) = j;
end
return