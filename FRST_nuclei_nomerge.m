% main script

% anisotropic diffusion filter H deconv image
clear; close all;

imgpath = '/Users/nathaning/Documents/MATLAB/Cedars/FRST/Eric ims';
load('/Users/nathaning/Documents/MATLAB/Cedars/FRST/Eric ims/files.mat');

d=datestr(now);
d = [strrep(d,':','-') ' Eric PCa ims r 3 5 11'];
outfolder = ['/Users/nathaning/Documents/MATLAB/Cedars/FRST/test results' filesep d];
% outfolder = '/Users/nathaning/Documents/MATLAB/Cedars/FRST/test results/40x';

if ~exist(outfolder, 'dir'), mkdir(outfolder);  end

tic;
% for NN= 1:length(tns)
for NN = [1 3 6]
% testname = ['test' num2str(NN)];
testname = tns{NN}(1:end-5);

I = imread([imgpath filesep testname '.tiff']);
I = colornormalization(I);
scale = 1;

H = Color_DeconvolutionPC(I);
H = H(:,:,1);
H_diff = anisodiff2D(H, 15, 1/7 , 30, 1);

ellips = struct([]);
Lmatrix = zeros(size(H,1), size(H,2), 1);
Lmat2 = zeros(size(H,1), size(H,2), 1, 'uint16');
n=1;

% colors = ['r', 'g', 'y', 'b', 'c']; i=1;
    for r = [3 5 9]; 
    SE = strel('disk', r);
    %open
    morph2 = imerode(H_diff,SE);
    morph = imreconstruct(morph2, H_diff);
    %close
    morph2 = imdilate(morph, SE);
    morph = imcomplement(imreconstruct(imcomplement(morph2),imcomplement(morph)));
    % morph(morph>210) = 255;

    SE = strel('disk', floor(r/2));
    %close again
    morph2 = imdilate(morph, SE);
    morph = imcomplement(imreconstruct(imcomplement(morph2),imcomplement(morph)));
    morph(morph>200) = 255;

    S = FRSTEX2(double(morph), [r:2:(r+5)], 0.05, 6.5, 1);
    [Gmag,~] = imgradient(double(morph));
    
    BW = imextendedmax(S,0.4);
    D = bwdist(BW);
    DL = watershed(D);
    bak = DL==0;
    BW = imdilate(BW,SE); %new
    Gmag(BW) = -Inf;
    Gmag(bak) = -Inf;
    L = watershed(Gmag);
    STATSfrst = regionprops(L, H, 'Solidity', 'Area', 'Centroid', 'MajorAxisLength',...
                               'MinorAxisLength', 'Orientation', 'WeightedCentroid',...
                               'MeanIntensity', 'Eccentricity');
    [regfrst, frstLL, Lfrst, STATSfrst] = Post_process(H, L, r, STATSfrst);
    
    
%     LL = fit_ellipses(Lfrst, STATSfrst);
    frstLL(frstLL~=0) = frstLL(frstLL~=0) + max(max(Lmatrix));
    Lmatrix = Lmatrix + frstLL;
    
    Lfrst(Lfrst~=0) = Lfrst(Lfrst~=0) + max(max(Lmat2));
    Lmat2 = Lmat2 + Lfrst;
    
    end

%% % ===========================================================================
% % Write images to file
LL = Lmatrix;
Hout = uint8(H_diff);
perims = bwperim(LL); perims(1,:)=0; perims(:,1)=0; perims(end,:)=0; perims(:,end)=0;
IoR = Hout; IoR(perims) = 255;
IoG = Hout; IoG(perims) = 0;
IoB = Hout; IoB(perims) = 0;
Io = cat(3,IoR, IoG, IoB);

Hout = I;
IoR = Hout(:,:,1); IoR(LL~=0) = 0;
IoG = Hout(:,:,2); IoG(LL~=0) = 255;
IoB = Hout(:,:,3); IoB(LL~=0) = 0;
Hout = cat(3,IoR, IoG, IoB);

% --H deconvolution with oulined elliptical regions
imwrite(Io, [outfolder filesep testname '_io.jpg'], 'jpg'); %

% --Original HE with filled in elliptical regions
% imwrite(Hout, [outfolder filesep testname '_fe.jpg'], 'jpg');

% --Label matrix converted to rgb - can remove label2rgb function for b/w
%   mask but the numbers are lost (all regions =255).
imwrite(LL, [outfolder filesep testname '_elps.png'], 'png');
% imwrite(double(Lmat2), [outfolder filesep testname '_bin.png'], 'png');

% fprintf(['\n' testname ' done' ])
% save([outfolder filesep testname '_LL.mat'], 'LL'); %LL is elliptical nuclei mask, 'nuclei' is raw regions


end

toc