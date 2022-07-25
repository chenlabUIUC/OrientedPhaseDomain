% NAME:  Calculate_strain_(g220)_(g002)
% PURPOSE:  This script is designed to calculate the local strain for
% cathode NPs along x_(g220) and y_(g002) directions based on reciprocal
% lattice vectors derived from 4D-STEM data
% INPUT:
%           g vector matrix processed from 4D-STEM data: '(g002).xlsx' and '(g220).xlsx'
% OUTPUT:
%           Strain matrix and tetragonality matrix: 'strain_(g002).xlsx', 'strain_(g220).xlsx', 'ratiobetween(g002)and(g220).xlsx'
% HISTORY:  written by Wenxiang Chen, 2019


% Input data
filename2= '(g002).xlsx'; % file with the absolute value of g002 in reciprocal space in the unit nm-1
filename3= '(g220).xlsx'; % file with the absolute value of g220 in reciprocal space in the unit nm-1
% Output data
filename4= 'strain_(g002).xlsx';
filename5= 'strain_(g220).xlsx';
filename6= 'ratiobetween(g002)and(g220).xlsx';

DA = xlsread(filename2); 
c1=0.84361;  % Input the lattice constant (nm) obtained from reference (pristine) particle for g002

gref002=2/c1;
strainA= (gref002-DA)./DA*100;
xlswrite (filename4, strainA); % save the strain along g002 direction (y direction)

DB = xlsread(filename3);
c2=0.82738;  % Input the lattice constant (nm) obtained from reference (pristine) particle for g220
gref220=sqrt(8)/c2;
strainB= (gref220-DB)./DB*100;
xlswrite (filename5, strainB); % save the strain along g220 direction (x direction)

ratio=DA./DB; % save the ratio between g002 and g220 as tetragonality
xlswrite (filename6, ratio);