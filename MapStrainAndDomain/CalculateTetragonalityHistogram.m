% NAME:     CalculateTetragonalityHistogram
% PURPOSE:  This code calculates the histogram of tetragonality for a
%           cathode particle
% INPUT:
%           Tetragonality matrix: 'ratiobetween(g002)and(g220).xlsx'
%           Mask: 'mask.tif'
% Note:     For generating the mask of the cathode particle, please refer to 
%           the particle morphology in the 4D-STEM data. One can manually delineate the
%           particle shape and make a mask using, for example, imageJ. An example of 
%           mask "mask.tif" is provided in the folder.
% OUTPUT:
%           Histogram of tetragonality: 'ratiobetween(g002)and(g220)_hist.xlsx'
% HISTORY:  written by Wenxiang Chen and Chang Qian, 2019


% Load data
filename = 'ratiobetween(g002)and(g220).xlsx';
filename2= 'mask.tif';

% Output data
filename4= 'ratiobetween(g002)and(g220)_hist.xlsx';

B = xlsread(filename);
D = imread('mask.tif');
[XX,YY] = size(B);

for x=1:1:XX
    for y=1:1:YY
        if D(x,y)==0
            B(x,y)=-20; % set the background as the lowest value for the easy of removing the background
        end
    end
end

% calcualte the histogram of tetragonality
E=B(:); % convert the tetragonality matrix to a column
edges = [0.61:0.003:0.763]; % define the edges of the bin in the histogram
centers = [0.6115:0.003:0.7615]; % define the centers of the bin in the histogram
N = histcounts(E,edges); % calculate the frequency count

% save the histogram of tetragonality
centers=transpose(centers);
N=transpose(N);
results = [centers,N];
xlswrite (filename4, results);