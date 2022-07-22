% NAME:  CreatPremask
% Purpose: This code is used to create a TIF image of strain map with grey scale 
% and with intensity between 0-1. The TIF image can be further used in ImageJ for
% manually delineating a mask.
% INPUT:
%           Strain matrix: 'strain_(g002).xlsx'
% OUTPUT:
%           'mask-Pre.tif'
% HISTORY:  written by Wenxiang Chen and Lehan Yao, 2019

% Load data
filename = 'strain_(g002).xlsx';

B = xlsread(filename);
[XX,YY] = size(B);

% In case that the original minimum/maxmimum value in B is too small/high, 
% it sets a low/high threshold here
for x=1:1:XX  
    for y=1:1:YY
        if B(x,y)<-4
            B(x,y)=-4; 
        end
        if B(x,y)>8
             B(x,y)=8; 
        end
    end
end

% Rescale the strain map to a range between 0 and 1.
D = B-min(min(B));
F = D/max(max(D)); 
imwrite(F, 'mask-Pre.tif');
