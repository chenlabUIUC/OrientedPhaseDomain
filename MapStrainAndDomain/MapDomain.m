% NAME:  MapDomain
% PURPOSE:  This code generates the domain map based on the fitting reuslts from tetragonality histogram
% INPUT:
%           Tetragonality matrix: 'ratiobetween(g002)and(g220).xlsx'
%           Mask: 'mask.tif'
% OUTPUT:
%           Domain map
% HISTORY:  written by Wenxiang Chen and Chang Qian, 2019

% Load data
filename = 'ratiobetween(g002)and(g220).xlsx';
filename2= 'mask.tif';

% Set the threshold from the fitting result here to distinguish the [100]t 
% and pesudo cubic phase
lowlimit = 0.65; 

% Set the threshold from the fitting result here to distinguish the [111]t 
% and pesudo cubic phase
highlimit = 0.724; 

count = 0;
B = xlsread(filename);
D = imread('mask.tif');
[XX,YY] = size(B);

% Generate the matrix of the domain map
for x=1:1:XX
    for y=1:1:YY
        if B(x,y)< lowlimit
            B(x,y)=2; % [100]t phase orientation
        elseif B(x,y)> highlimit
                B(x,y)=1; % [111]t phase orientation
        else B(x,y)=0; % Pseudo cubic phase [110]c
            
        end
    end
end

% The background region is set as a negative value
for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            B(x,y)=-1; 
        end
    end
end

% Plot the domain map
figure
P=imagesc(B);
caxis([-1,2])
load('cmapDomain.mat','mycmap');ax = gca; colormap(ax,mycmap);
