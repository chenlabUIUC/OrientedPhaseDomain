% NAME:  MapStrainforXandY
% PURPOSE:  This code maps strain along x and y directions without background
% INPUT:
%           Strain matrix: 'strain_(g220).xlsx' and 'strain_(g002).xlsx'
%           Mask: 'mask.tif'
% OUTPUT:
%           Strain maps
% HISTORY:  written by Wenxiang Chen, 2019

% Load strain matrix data along x and y direction and mask
filename1 = 'strain_(g220).xlsx';
filename2 = 'strain_(g002).xlsx';
filename3= 'mask.tif';

B1 = xlsread(filename1);
B2 = xlsread(filename2);
D = imread('mask.tif');
[XX,YY] = size(B1);

% set the background value as the lowest, -20
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            B1(x,y)=NaN;
            B2(x,y)=NaN;
        end
    end
end

% Set the low limit of the strain for the display only 
%{
for x=1:1:XX
    for y=1:1:YY
        if B1(x,y)<-2.6 & B1(x,y)~=-20
            B1(x,y)=-2.6; 
        end
        if B2(x,y)<-2.6 & B2(x,y)~=-20
            B2(x,y)=-2.6; 
        end
    end
end
%}

%Map strain along x direction
load('cmapStrain.mat','mycmap');
number=size(mycmap,1)-2
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

%Map strain along y direction
load('cmapStrain.mat','mycmap');
number=size(mycmap,1)-2
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
