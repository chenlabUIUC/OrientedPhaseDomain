% NAME:  FitTetragonalityHist
% PURPOSE:  This code opens the curve fitting tool for fitting the tetragonality histogram
% INPUT:
%           Tetragonality matrix: 'ratiobetween(g002)and(g220)_hist.xlsx'
% OUTPUT:
%           Fitting results
% HISTORY:  written by Wenxiang Chen and Chang Qian, 2019

% Load the data
filename = 'ratiobetween(g002)and(g220)_hist.xlsx';
B = xlsread(filename);

pos_x=B(:,1);
pos_y=B(:,2);

% An example of the fitting of the histogram is provided here:
[result, good] = createFit(pos_x, pos_y);
% Check the fitting results
result.b2  % threshold to distinguish the [100]t and pesudo cubic phase
result.b4  % threshold to distinguish the [111]t and pesudo cubic phase. 

% Note: One can perform fitting for cathode particles using the cftool as follows.
%       Fitting was performed with constraints, considering the theoretical 
%       values of tetragonality: R=0.613 for [100]t, R=0.707 for [110]c, R=0.756 
%       for [111]t. Please also refer to our manuscript.

% cftool(pos_x,pos_y);  % curve fitting tool