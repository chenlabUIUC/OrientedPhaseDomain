% NAME:     MapStrainforXandY
% PURPOSE:  This code maps strain along x and y directions without background
% INPUT:
%           Strain matrix: 'strain_(g220).xlsx' and 'strain_(g002).xlsx'
%           Mask of the cathode particle: 'mask.tif'
% Note:     For generating the mask of the cathode particle, please refer to 
%           the particle morphology in the 4D-STEM data. One can manually delineate the
%           particle shape and make a mask using, for example, imageJ. An example of 
%           mask "mask.tif" is provided in the folder.
% OUTPUT:
%           Strain maps
% HISTORY:  written by Lehan Yao and Wenxiang Chen, 2022


% Load strain matrix data along x and y direction and the mask for the cathode particle
filename1 = 'strain_(g220).xlsx';
filename2 = 'strain_(g002).xlsx';
filename3= 'mask.tif';

B1 = xlsread(filename1);
B2 = xlsread(filename2);
D = imread('mask.tif');
[XX,YY] = size(B1);

% set the background value outside of the cathode particle as NaN
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            B1(x,y)=NaN;
            B2(x,y)=NaN;
        end
    end
end

%Map strain along x direction
load('cmapStrain.mat','mycmap');
number=size(mycmap,1)-1
colorRange = [-3, 16];
B1_RGB = (B1-colorRange(1))/(colorRange(2)-colorRange(1));
B1_RGB = uint8(B1_RGB*number+1); %%%%%62
B1_RGB = ind2rgb(B1_RGB,mycmap(2:end,:));
B1_RGB(~logical(cat(3,D,D,D))) = 1;
figure
set(gca,'position',[0.1 0.1 0.8 0.8])
imshow(B1_RGB)
colormap(mycmap)
caxis([-3,16])
colorbar()
title('Strain_x_x')


%Map strain along y direction
load('cmapStrain.mat','mycmap');
number=size(mycmap,1)-1
colorRange = [-3, 16];
B2_RGB = (B2-colorRange(1))/(colorRange(2)-colorRange(1));
B2_RGB = uint8(B2_RGB*number+1); %%%%%62
B2_RGB = ind2rgb(B2_RGB,mycmap(2:end,:));
B2_RGB(~logical(cat(3,D,D,D))) = 1;
figure
set(gca,'position',[0.1 0.1 0.8 0.8])
imshow(B2_RGB)
colormap(mycmap)
caxis([-3,16])
colorbar()
title('Strain_y_y')
