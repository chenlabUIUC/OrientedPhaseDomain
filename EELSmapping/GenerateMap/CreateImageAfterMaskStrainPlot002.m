% NAME:  MapEELSData
% PURPOSE:  This code maps white line ratio without background
% INPUT:
%           White line ratio matrix: 'WLR matrix.xlsx'
%           Mask: 'mask.tif'
% OUTPUT:
%           White line ratio maps
% HISTORY:  written by Wenxiang Chen and Lehan Yao, 2022

% Load EELS white line intensity ratio (L3/L2) data. Load the mask for the cathode particle
% Note: For how to obtain L3/L2 data from raw EELS data, refer to the code "Hyperspy EELS Mn White Line Ratio.ipynb"
filename = 'WLR matrix.xlsx';
filename2= 'mask.tif';

B = xlsread(filename);
D = imread('mask.tif');
[XX,YY] = size(B);

% set the background value outside of the cathode particle as NaN
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            B(x,y)=NaN;
        end
    end
end

% Map white line ratio in the cathode particle
load('cmapEELS.mat','mycmap');
number=size(mycmap,1)-1
colorRange = [2, 4];
B1_RGB = (B-colorRange(1))/(colorRange(2)-colorRange(1));
B1_RGB = uint8(B1_RGB*number+1); %%%%%62
B1_RGB = ind2rgb(B1_RGB,mycmap(2:end,:));
B1_RGB(~logical(cat(3,D,D,D))) = 1;
figure
set(gca,'position',[0.1 0.1 0.8 0.8])
imshow(B1_RGB)
colormap(mycmap)
caxis([2,4])
colorbar()

% Optional: save the white line intensity ratio L3/L2 as excel file
%{
n=1;
result=[];
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)~=0
            result(n,1)=Backup(x,y);  % save the L3/L2 ratio
            n=n+1;
        end
    end
end
cc=mean(result); % save the average L3/L2 ratio in the second column
s=std(result);   % save the standard deviation of L3/L2 ratio in the second column
result(1,2)=cc;
result(2,2)=s;
filename3= 'L3L2ratio(noBG_column).xlsx';
xlswrite (filename3, result);
%}
