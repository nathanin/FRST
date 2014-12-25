function S = FRSTEX2(im,radii,thresh,alpha,orient,A)
% FRST   compute Fast Radial Symmetry Transform (FRST).
%     Computes the Fast Radial Symmetry Transform (FRST) as described in 
%       
%       Loy & Zelinsky (2003), 
%       Fast Radial Symmetry for Detecting Points of Interest,
%       IEEE Transactions on Pattern Analysis and Machine Intelligence, 
%       August 2003.
%
%     This code returns the combined absolute magnitudes of the dark and 
%     light symmetry components which is well suited for generic interest 
%     point detection (this varies slightly from the combined output presented
%     in the papers where dark symmetry is shown as negative and light as 
%     positive values).
%
%     To restrict this code to only detect dark or light comment out
%     the 3 lines suffixed with "% DISABLE FOR DARK/LIGHT".  With these
%     lines disabled dark symmetry is computed for positive radii
%     parameter values and light symmetry for negative radii values.  (The
%     resulting dark and light symmetries can be combined to generate the
%     combined output presented in the papers.)
%
%     Parameters:
%       im     = input image (gray scale)
%       radii  = vector of radii at which to compute transform
%       thresh = optional gradient threshold (beta in ECCV and PAMI papers), 
%                default 0.05.
%       A      = optional set of 1D separable component of 2D blurring kernel, default
%                A{i} = gaussian size 2*radii with sd radii/2.
%       alpha  = optional radial strictness parameter, default 2.
%       orient = optional flag (1 or 0) for computing 'orientation symmetry', i.e. only
%                using gradient orientations and not magnitude, see PAMI paper,
%                default 0.
%
%     Examples:
%       S = FRST(im, [2 4 7 10]);
%       S = FRST(double(imread('g.jpg')), 5);
%
%
% release 1.0.
% Gareth Loy, KTH, Stockholm, 2005.
%
%--------------------------------------------------------------------------------
% The code was modified by Arkadiusz Gertych PhD, to improve speed,
% and rule out utilization of external mex code provided by G.Loy
% By A.Gertych PhD, Bioinformatics Laboratory at C-S, Los Angeles, 2010-2011.  
%--------------------------------------------------------------------------------

% set default values for undefined parameters
% Product transform smoothing mask %
% Filter mask size needs to be somehow modulated, 
% because for large foci the smoothing effect is to strong, and signals may fall below the threshold. 
% Modification by A.Gertych Oct-12-2011

% Filter settings. Remark: setting sigma >= radii(i) causes "boxing" effect
% on symmetry image. Keep sigma <0.5*radii(i), and change filter size to
% achive desired smoothing effect.

if nargin <= 5,
    parfor i = 1:length(radii),
        if radii(i) <= 9,
            A{i} = make_gauss_vec(round(2*radii(i)), 0.5*radii(i));
        else
%             fprintf('%s \n',['Stronger smoothing was applied for R=' num2str(radii(i))]);
            A{i} = make_gauss_vec(round(3*radii(i)), 0.5*radii(i));
        end
    end
end


if nargin<=4, orient = 0; end
if nargin<=3, alpha = 2; end
if nargin<=2, thresh = 0.05; end

% determine gradient magnitude and unit gradient
grad = gradient(im);
Mag = sqrt(grad.x.^2 + grad.y.^2);
unit_Grad.x = grad.x./(Mag + 1e-5);
unit_Grad.y = grad.y./(Mag + 1e-5);

% set edge mag to zero to avoid edge effects
Mag([1:2, end-2:end], :) = 0;
Mag(:, [1:2, end-2:end]) = 0;

% pad with zeros to accommodate voting outside the image
frame_r = max(abs(radii));
unit_Grad.x = pad(unit_Grad.x, frame_r, frame_r);
unit_Grad.y = pad(unit_Grad.y, frame_r, frame_r);
Mag = pad(Mag, frame_r, frame_r);

[Coords_x, Coords_y] = meshgrid(1:size(unit_Grad.x,2), 1:size(unit_Grad.x,1));

% only consider significant gradient elements
ind = find(Mag > thresh);   % pixels of significant gradient magnitide
num = length(ind);          % number of such pixels

affected_pix_x = zeros(length(radii)*num,1);  %temporary matrix
affected_pix_y = zeros(length(radii)*num,1);  %temporary matrix
affecting_ind = repmat(ind,length(radii),1);  %temporary matrix: replicates indices mtx with length(radii)rows and no columns

% if there are significant gradient elements present
if ~isempty(ind)
    S = 0;
    for i = 1:length(radii)

        % compute affected pixels and indices of affecting gradient
        % elements                
        %%%%%%% !!!!!!!!!!!!!!!!! add + and - combination
        for r = 1:length(radii)           
            affected_pix_x((r-1)*num+1:r*num) = round(Coords_x(ind) + unit_Grad.x(ind)*radii(r)); % x - coordinates of affected pixels
            affected_pix_y((r-1)*num+1:r*num) = round(Coords_y(ind) + unit_Grad.y(ind)*radii(r)); % y - coordinates of affected pixels
        end

        affected_ind = (affected_pix_x-1)*size(unit_Grad.x,1) + affected_pix_y;
        
        %the line below can replace the line above., not sure if this will be faster: AGertych 07/2011
        %affected_ind =  sub2ind(size(unit_Grad.x),affected_pix_y, affected_pix_x);
        
        
        % compute orientation image,
        hist_1D = histc(affected_ind, 1:size(unit_Grad.x,1)*size(unit_Grad.x,2)); %does a histogram of indices
        Oa_1D = abs(hist_1D).^(alpha);
        %indx = Oa_1D < 10;          Oa_1D(indx) = 0;
        Oa = reshape(Oa_1D, size(unit_Grad.x,1), size(unit_Grad.x,2));

        if orient~=1
            % compute magnitude image, and combine with orientation image
            
            % hist_1D_M = histc_weighted(affected_ind, Mag(affecting_ind(:)), [1:size(unit_Grad.x,1)*size(unit_Grad.x,2)]);
            
            % the line above will not work on 64bit OS and Matlab, because of the proprietary function histc_weighted
            % thus the weighted histogram calulation is replaced by the two lines
            % below. It is also faster than the original code. 
            % BY A.GERTYCH 2010-Nov-9
            % 
            
             padds = zeros((size(unit_Grad.x,1)*size(unit_Grad.x,2)) - max(affected_ind(:)),1,'double');
             hist_1D_M = [accumarray(affected_ind, Mag(affecting_ind(:))); padds];
            
             Ma = reshape(hist_1D_M, size(unit_Grad.x,1), size(unit_Grad.x,2));
             Oa = Oa./(10^alpha).*(Ma/10);   % kn = 10 is expreimentally adjusted
        end

        % remove frame of zeros
        Oa = unpad(Oa, frame_r, frame_r);

        % Use a separable kernel to blur result
        S = S + conv2(A{i}, A{i}', Oa, 'same');
        %S = S + Oa;
    end
else
    % return an image of zeros
    warning('IND is empty, returning an image of zeros');
    S = zeros(size(unit_Grad.x)-2*frame_r);
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE FUNCTIONS BELOW ARE GENERIC LIBRARY FUNCTIONS. THEY ARE INCLUDED IN 
% THIS FILE HERE AS THEY ARE REQUIRED BY THE FRST FUNCTION, HOWEVER, THEY 
% ARE USEFUL ON THEIR OWN AND AS SUCH CAN BE MOVED TO INDIVIDUAL FILES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function grad = gradient(I)
% grad = GRADIENT(I)
%
% determines the gradient of I returns the result as a vector field.
% 
% INPUT:  I = grayscale image
%
% OUTPUT: grad.x = x gradient image
%         grad.y = y gradient image
% 
% e.g.    grad = gradient(I);
% 

% Gareth Loy, KTH, Stockholm, 2003-2005.

if size(I,3) > 1
    error('gradient.m: unable to compute gradient of colour images.');
end

g = make_gauss_vec(5,1); 
h = [-1; -2; 0; 2; 1];
inv_step_edge_response = 0.4466; 

% g = [1; 1; 1];
% h = [-1; 0; 1];
% inv_step_edge_response = 0.1768;

Ix = conv2(I,h','same');            
Iy = conv2(I,h,'same');            
Ix = conv2(Ix,g,'same');            
Iy = conv2(Iy,g','same');            
% Ix = conv2(g,g,Ix,'same');            
% Iy = conv2(g,g,Iy,'same');            
grad.x = Ix*inv_step_edge_response;
grad.y = Iy*inv_step_edge_response;

return


function mx_new = pad(mx,x,y,val)
%
% Pads a mx with a frame x wide and y high consisting of the value val.
% Default value for 'val' is 0.
%
% E.g. to pad a matrix mx with a border x_rad and y_rad of zeros:
%  pad(mx,x_rad,y_rad,0);
%  pad(mx,x_rad,y_rad);
%

% Gareth Loy, KTH, Stockholm, 2003-2005.

if nargin < 4
    val = 0;
end
mx_new = repmat(val,[2*[y,x]+[size(mx,1),size(mx,2)], size(mx,3)]);
mx_new(y+1:end-y,x+1:end-x,:) = mx;


function mx = unpad(mx,x,y)
%
% Unpads matrix mx by removing the outer x-wide and y-high frame, i.e.
% "cropping" the matrix.
%
% E.g. to pad a matrix mx with a border x_rad and y_rad of zeros:
%  unpad(max,x_rad,y_rad);
%

% Gareth Loy, KTH, Stockholm, 2003-2005.

mx = mx(y+1:end-y, x+1:end-x, :);

return


function [gauss_vec,d_gauss_vec] = make_gauss_vec(n,sd)
% gauss_vec = make_gauss_vec(n,sd)
%
% returns a normalised 1D Gaussian mask with
%  - n elements, and 
%  - sd standard deviation.
% 
% [gauss_vec,d_gauss_vec] = make_gauss_vec(n,sd)
%  also returns the derivative of the gaussian 
%

% Gareth Loy, KTH, Stockholm, 2003-2005.

for i=1:n
    dist = i-(n+1)/2;
    gauss_vec(i) = exp(-0.5*(dist/sd)^2)  / (sqrt(2*pi)*sd);
end; %for
gauss_vec = gauss_vec/(sum(sum(gauss_vec)));

if nargout==2
    for i=1:n
        dist = i-(n+1)/2;
        d_gauss_vec(i) = -2*dist/(2*sd)/(sqrt(2*pi)*sd) * exp(-0.5*(dist/sd)^2);
    end
    d_gauss_vec = d_gauss_vec/(sum(sum(gauss_vec)));
end

return
