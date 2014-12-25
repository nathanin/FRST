function [newpixels] = Color_DeconvolutionPC(pixels)
 %img newpixels

    % load image
    %[name path]=uigetfile('*.tif; *jpg','Select TIF file:');
    %pixels=imread([path name]);
    
    % image size
    [height,width,~]=size(pixels);
    
    %HAem
%	MODx(1)= 0.644211; %//0.650;
%	MODy(1)= 0.716556; %//0.704;
%	MODz(1)= 0.266844; %//0.286;    
    
%    MODx(1)= 0.49015734;
%	MODy(1)= 0.76897085;
%	MODz(1)= 0.41040173;    
    
    %FastRed
 %   MODx(2)= 0.21393921;
%	MODy(2)= 0.85112669;
%	MODz(2)= 0.47794022;    
    
    %DAB
	MODx(3)= 0.268;
	MODy(3)= 0.570;
	MODz(3)= 0.776;    
    
    
    % GL Haem matrix
    MODx(1)= 0.644211; %0.650;
    MODy(1)= 0.716556; %0.704;
    MODz(1)= 0.266844; %0.286;
    % GL Eos matrix
    MODx(2)= 0.092789; %0.072;
    MODy(2)= 0.954111; %0.990;
    MODz(2)= 0.283111; %0.105;

    % Zero matrix
%     MODx(3)= 0.6359544;
%     MODy(3)= 0.0010;
%     MODz(3)= 0.7717266;
    

    % start
    cosx=zeros(3,1);cosy=zeros(3,1);cosz=zeros(3,1);
    len=zeros(3,1);
    for i=1:3
        %normalise vector length
        cosx(i)=0;
        cosy(i)=0;
        cosz(i)=0;
        len(i)=sqrt(MODx(i)*MODx(i) + MODy(i)*MODy(i) + MODz(i)*MODz(i));
        if (len(i) ~= 0.0)
            cosx(i)= MODx(i)/len(i);
            cosy(i)= MODy(i)/len(i);
            cosz(i)= MODz(i)/len(i);
        end
    end
    
    
    % translation matrix
    if (cosx(2)==0.0) %2nd colour is unspecified
        if (cosy(2)==0.0)
            if (cosz(2)==0.0)
                cosx(2)=cosz(1);
                cosy(2)=cosx(1);
                cosz(2)=cosy(1);
            end
        end
    end
    
    if (cosx(3)==0.0)  % 3rd colour is unspecified
        if (cosy(3)==0.0)
            if (cosz(3)==0.0)
                if ((cosx(1)*cosx(1) + cosx(2)*cosx(2))> 1)
                    cosx(3)=0.0;
                else
                    cosx(3)=sqrt(1.0-(cosx(1)*cosx(1))-(cosx(2)*cosx(2)));
                end
                
                if ((cosy(1)*cosy(1) + cosy(2)*cosy(2))> 1)
                    cosy(3)=0.0;
                else
                    cosy(3)=sqrt(1.0-(cosy(1)*cosy(1))-(cosy(2)*cosy(2)));
                end
                
                if ((cosz(1)*cosz(1) + cosz(2)*cosz(2))> 1)
                    cosz(3)=0.0;
                else
                    cosz(3)=sqrt(1.0-(cosz(1)*cosz(1))-(cosz(2)*cosz(2)));
                end
            end
        end
    end
    
    leng=sqrt(cosx(3)*cosx(3) + cosy(3)*cosy(3) + cosz(3)*cosz(3));
    
    cosx(3)= cosx(3)/leng;
    cosy(3)= cosy(3)/leng;
    cosz(3)= cosz(3)/leng;
    
    for i=1:3
        if (cosx(i) == 0.0), cosx(i) = 0.001; end
        if (cosy(i) == 0.0), cosy(i) = 0.001; end
        if (cosz(i) == 0.0), cosz(i) = 0.001; end
    end
    
    
    %matrix inversion
    A = cosy(2) - cosx(2) * cosy(1) / cosx(1);
    V = cosz(2) - cosx(2) * cosz(1) / cosx(1);
    C = cosz(3) - cosy(3) * V/A + cosx(3) * (V/A * cosy(1) / cosx(1) - cosz(1) / cosx(1));
    q(3) = (-cosx(3) / cosx(1) - cosx(3) / A * cosx(2) / cosx(1) * cosy(1) / cosx(1) + cosy(3) / A * cosx(2) / cosx(1)) / C;
    q(2) = -q(3) * V / A - cosx(2) / (cosx(1) * A);
    q(1) = 1.0 / cosx(1) - q(2) * cosy(1) / cosx(1) - q(3) * cosz(1) / cosx(1);
    q(6) = (-cosy(3) / A + cosx(3) / A * cosy(1) / cosx(1)) / C;
    q(5) = -q(6) * V / A + 1.0 / A;
    q(4) = -q(5) * cosy(1) / cosx(1) - q(6) * cosz(1) / cosx(1);
    q(9) = 1.0 / C;
    q(8) = -q(9) * V / A;
    q(7) = -q(8) * cosy(1) / cosx(1) - q(9) * cosz(1) / cosx(1);
    
    newpixels=zeros(height,width,3);
    
    for j=1:width
        for i=1:height
            % log transform the RGB data
            R=pixels(i,j,1);G=pixels(i,j,2);B=pixels(i,j,3);
            Rlog = -((255.0*log((double(R)+1)/255.0))/log(255));
            Glog = -((255.0*log((double(G)+1)/255.0))/log(255));
            Blog = -((255.0*log((double(B)+1)/255.0))/log(255));
            for index=1:3
                % rescale to match original paper values
                Rscaled = Rlog * q((index-1)*3+1);
                Gscaled = Glog * q((index-1)*3+2);
                Bscaled = Blog * q(index*3);
                output = exp(-((Rscaled + Gscaled + Bscaled) - 255.0) * log(255) / 255.0);
                if(output>255), output=255; end
                newpixels(i,j,index)=uint8(floor(output+.5));
            end
        end
    end
    
    newpixels=uint8(newpixels);
    %figure;imshow(newpixels(:,:,1));  %Haematoxylin image
    %figure;imshow(newpixels(:,:,2));  %Eosine image
    %figure;imshow(newpixels(:,:,3));  %Null image
    
    
    
    
end
