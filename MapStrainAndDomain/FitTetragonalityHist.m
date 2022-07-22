% NAME:  FitTetragonalityHist
% PURPOSE:  This code opens the curve fitting tool for fitting the tetragonality histogram
% INPUT:
%           Tetragonality matrix: 'ratiobetween(g002)and(g220)_hist.xlsx'
% HISTORY:  written by Wenxiang Chen and Chang Qian, 2019

% Load the data
filename = 'ratiobetween(g002)and(g220)_hist.xlsx';
B = xlsread(filename);

pos_x=B(:,1);
pos_y=B(:,2);

cftool(pos_x,pos_y);  % curve fitting tool

% The fitting was performed with certain constraints. Please refer to our
% manuscript for details.