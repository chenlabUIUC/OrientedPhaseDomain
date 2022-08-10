% NAME:     MapStrainGradientMagnitude
% PURPOSE:  This code maps the magnitude of strain gradient in the cathode nanoparticle
% INPUT:
%           Strain matrix: 'strain_(g002).xlsx'
%           Mask of the cathode particle: 'mask.tif'
% Note:     For generating the mask of the cathode particle, please refer to 
%           the particle morphology in the 4D-STEM data. One can manually delineate the
%           particle shape and make a mask using, for example, imageJ. An example of 
%           mask "mask.tif" is provided in the folder.
% OUTPUT:
%           Strain gradient map (unit: % per pixel)
% HISTORY:  written by Lehan Yao and Wenxiang Chen, 2022

% Load strain matrix data and strain mask
filename = 'strain_(g002).xlsx';
B = xlsread(filename);
D = imread('mask.tif');
DD=D;
[XX,YY] = size(B);

% Set the background outside the nanoparticle as -5
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            B(x,y)=-5;   
        end
    end
end

% Calculate the gradient of strain, including vectors and magnitude
[px,py] = gradient(B);
pxy=sqrt(px.^2+py.^2);

% For the background outside the nanoparticle, set the strain gradient
% magnitude as -5 and gradient vectors as 0, respectively.
[XX,YY] = size(B);
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            pxy(x,y)=-5;   % Background outside the nanoparticle is set as -5
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
        if x==1 | y==1 | x==XX | y==YY
            pxy(x,y)=-5;
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
    end
end

% To avoid the sharp gradients at the boundary of the particle, set the
% strain gradient in the pixels close to the particle boundary (within 1 pixel) as 0.
XX=XX-1;
YY=YY-1;
for x=2:1:XX
    for y=2:1:YY
        if D(x-1,y)==0 | D(x+1,y)==0 | D(x,y-1)==0 | D(x,y+1)==0
            pxy(x,y)=-5;
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
    end
end

% To avoid the sharp gradients at the boundary of the particle, set the
% strain gradient in the pixels close to the particle boundary (within 2 pixels) as 0.
XX=XX-1;
YY=YY-1;
for x=3:1:XX
    for y=3:1:YY
        if D(x-2,y)==0 | D(x+2,y)==0 | D(x,y-2)==0 | D(x,y+2)==0
            pxy(x,y)=-5;
            px(x,y)=0;
            py(x,y)=0;
            DD(x,y)=0;
        end
    end
end

% Map the strain gradient magnitude
figure
P=imagesc(pxy);
caxis([-1 6])
load('cmapGradient.mat','mycmap');
colormap(mycmap)
axis image;

% Calculate the mean value of strain gradient magnitude
[XX,YY] = size(B);
n=1;
for x=1:1:XX
    for y=1:1:YY
        if DD(x,y)~=0
            result(n,1)=pxy(x,y);
            n=n+1;
        end
    end
end
mean(result)